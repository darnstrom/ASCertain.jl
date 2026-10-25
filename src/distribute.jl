function _normalize_worker_ids(workers::Vector{Int})
    return [pid for pid in unique(workers) if pid != Distributed.myid()]
end

function validate_distributed_settings(opts::CertSettings)
    isempty(opts.rm_callbacks) || error("Distributed certification requires empty rm_callbacks.")
    isempty(opts.add_callbacks) || error("Distributed certification requires empty add_callbacks.")
    isempty(opts.termination_callbacks) || error("Distributed certification requires empty termination_callbacks.")
    isempty(opts.pop_callbacks) || error("Distributed certification requires empty pop_callbacks.")
    isempty(opts.conditioned_callbacks) || error("Distributed certification requires empty conditioned_callbacks.")
    opts.overflow_handle === ASCertain.default_overflow_handle || error("Distributed certification requires the default overflow handler.")
    opts.output_limit == CertSettings().output_limit || error("Distributed certification requires the default output_limit.")
end

# Breath first search to generate a initial partition to distribute
function seed_regions(R0::Region, prob::CertProblem, P_theta, opts::CertSettings, target_count::Int64)
    ws = setup_workspace(P_theta,opts.max_constraints)
    settings(ws.DAQP_workspace,opts.daqp_settings)
    ws.warm_start = opts.warm_start_ldp
    pending, partition = Region[R0], Region[]
    iter_max, N_fin, lp_count = 0,0,0
    try
        while(length(pending) < target_count && !isempty(pending))
            region = popfirst!(pending)
            reset_workspace(ws)
            S = Region[]
            parametric_AS_iteration(prob,region,opts,ws,S)
            materialize_regions!(S,ws)
            append!(pending,S)
            append!(partition,ws.F)
            iter_max = max(iter_max,ws.iter_max)
            N_fin += ws.N_fin
            lp_count += ws.lp_count
        end
    finally
        DAQP.free_c_workspace(ws.DAQP_workspace)
    end
    return pending, partition, iter_max, N_fin, lp_count
end

function distributed_region_priority(region::Region)
    # TODO replace this with better heuristic
    # Maybe sampling a QP and using # of iterations?
    # Note that "less complex" should be the metric
    return (region.iter, -count(region.IS), -length(region.AS), -size(region.Ath,2))
end

function ensure_distributed_workers!(worker_ids::Vector{Int})
    @sync for pid in worker_ids # (in parallel)
        @async Distributed.remotecall_eval(Main, pid, :(import ASCertain))
    end
end

# One workspace per worker process, reused by all its tasks
const WORKER_WS = Ref{Union{Nothing,CertWorkspace}}(nothing)
# The problem that the worker certifies regions of: (prob, P_theta, opts)
const WORKER_PROBLEM = Ref{Any}(nothing)

function set_worker_problem!(prob::CertProblem, P_theta, opts::CertSettings)
    WORKER_PROBLEM[] = (prob, P_theta, opts)
    return nothing
end

function worker_workspace(P_theta, opts::CertSettings)
    ws = WORKER_WS[]
    max_radius = length(P_theta.ub)*(maximum(P_theta.ub)^2)/2 # (as in setup_workspace)
    if(isnothing(ws) || size(ws.Ath) != (length(P_theta.ub),opts.max_constraints) || ws.max_radius != max_radius)
        isnothing(ws) || DAQP.free_c_workspace(ws.DAQP_workspace)
        ws = WORKER_WS[] = setup_workspace(P_theta,opts.max_constraints)
    end
    ws.theta_lb .= P_theta.lb; ws.theta_ub .= P_theta.ub
    settings(ws.DAQP_workspace,opts.daqp_settings)
    ws.warm_key = (-1,-1)
    ws.F = Region[]; ws.ASs = falses(0,0); empty!(ws.ASs_state); empty!(ws.bin)
    ws.lp_count = 0; ws.n_witness = 0; ws.i_witness = 0; ws.N_fin = 0; ws.iter_max = 0
    return ws
end

# Certify the regions (on a worker). If this takes longer than time_budget, the remaining regions are
# returned (in explicit form) so that they can be redistributed.
function distributed_region_certify(regions::Vector{Region}, time_budget::Float64)
    prob, P_theta, opts = WORKER_PROBLEM[]
    ws = worker_workspace(P_theta,opts)
    deadline = time() + time_budget
    opts.conditioned_callbacks = [((S,prob,ws,opts) -> time() > deadline, (S,prob,ws,opts) -> true)]
    S = copy(regions)
    t = @elapsed part, iter_max, N_fin, _, _, _ = certify(S,prob,ws,opts)
    materialize_regions!(S,ws)
    return part, iter_max, N_fin, ws.lp_count, t, S
end

function distributed_certify(R0::Region, prob::CertProblem, P_theta, opts::CertSettings, worker_ids::Vector{Int})
    ts = Float64[]
    seed_target = max(length(worker_ids), opts.distributed_region_factor*length(worker_ids))

    t0 = @elapsed pending, part, iter_max, N_fin, lp_count = seed_regions(R0,prob,P_theta,opts,seed_target)
    push!(ts,t0)
    isempty(pending) && begin
        ASs, ASs_state = rebuild_stored_ASs(part,prob,opts)
        print_final_information(iter_max,lp_count,N_fin,size(ASs,2),opts)
        return part, iter_max, N_fin, ASs, Dict{Any,Any}(), ASs_state
    end

    ensure_distributed_workers!(worker_ids)
    worker_opts = deepcopy(opts)
    worker_opts.verbose = 0
    worker_opts.store_ASs = false
    @sync for pid in worker_ids
        @async Distributed.remotecall_wait(set_worker_problem!, pid, prob, P_theta, worker_opts)
    end

    # Dynamic scheduling: idle workers take the next job, and the regions that a job did not finish
    # within the time budget are queued as new jobs
    jobs = Channel{Union{Nothing,Vector{Region}}}(Inf)
    for region in sort(pending; by=distributed_region_priority)
        put!(jobs,[region])
    end
    n_jobs = length(pending) # Queued or running
    @sync for pid in worker_ids
        @async while true
            job = take!(jobs)
            isnothing(job) && break
            worker_part, worker_iter_max, worker_N_fin, worker_lp_count, t, remaining = try
                Distributed.remotecall_fetch(distributed_region_certify, pid, job, opts.distributed_time_budget)
            catch err
                close(jobs) # Do not leave the other workers waiting for jobs
                rethrow(err)
            end
            append!(part,worker_part)
            iter_max = max(iter_max,worker_iter_max)
            N_fin += worker_N_fin
            lp_count += worker_lp_count
            push!(ts,t)
            for region in remaining
                put!(jobs,[region])
            end
            n_jobs += length(remaining)-1
            if(n_jobs == 0) # Everything is certified -> stop the workers
                foreach(_ -> put!(jobs,nothing), worker_ids)
            end
        end
    end

    ASs, ASs_state = rebuild_stored_ASs(part,prob,opts)
    print_final_information(iter_max,lp_count,N_fin,size(ASs,2),opts)
    return part, iter_max, N_fin, ASs, Dict{Any,Any}(:ts=>ts), ASs_state
end
