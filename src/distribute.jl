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
    for pid in worker_ids
        Distributed.remotecall_eval(Main, pid, :(import ASCertain))
    end
end

function distributed_region_certify(prob::CertProblem, P_theta, region::Region, opts::CertSettings)
    ws = setup_workspace(P_theta,opts.max_constraints)
    settings(ws.DAQP_workspace,opts.daqp_settings)
    try
        t = @elapsed part, iter_max, N_fin, _, _, _ = certify(region,prob,ws,opts)
        return part, iter_max, N_fin, ws.lp_count, t
    finally
        DAQP.free_c_workspace(ws.DAQP_workspace)
    end
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
    ordered_regions = sort(copy(pending); by=distributed_region_priority)
    pool = Distributed.CachingPool(worker_ids)
    results = Distributed.pmap(region -> distributed_region_certify(prob,P_theta,region,worker_opts), pool, ordered_regions, batch_size=opts.distributed_batch_size)

    for (worker_part, worker_iter_max, worker_N_fin, worker_lp_count,t) in results
        append!(part,worker_part)
        iter_max = max(iter_max,worker_iter_max)
        N_fin += worker_N_fin
        lp_count += worker_lp_count
        push!(ts,t)
    end

    ASs, ASs_state = rebuild_stored_ASs(part,prob,opts)
    print_final_information(iter_max,lp_count,N_fin,size(ASs,2),opts)
    return part, iter_max, N_fin, ASs, Dict{Any,Any}(:ts=>ts), ASs_state
end
