## ASCertification main
function certify(mpQP,P_theta,AS::Vector{Int64}=Int64[];opts::CertSettings=CertSettings(),normalize=true, ws=nothing, workers::Vector{Int}=Int[])
    prob = setup_certproblem(mpQP;normalize);
    prob,P_theta,mpQP = ASCertain.normalize(prob,deepcopy(P_theta),deepcopy(mpQP));
    return certify(prob,P_theta,copy(AS),opts;ws,workers)
end

function certify(prob::CertProblem,P_theta,AS::Vector{Int64},opts::CertSettings;ws = nothing, workers::Vector{Int}=Int[])
    nth = length(P_theta.ub);

    # Pack initial region in the form A'θ ≤ b
    A = [Matrix{Float64}(I,nth,nth) Matrix{Float64}(-I,nth,nth) P_theta.A];
    b = [P_theta.ub;-P_theta.lb;P_theta.b];
    R0 = Region(AS,A,b,prob);
    if(opts.compute_flops)
        R0.kappa[:flops]= [0,0,0,0];
    end

    worker_ids = _normalize_worker_ids(workers)
    if(!isempty(worker_ids))
        isnothing(ws) || error("Distributed certification does not support a preallocated workspace.")
        validate_distributed_settings(opts)
        return distributed_certify(R0,prob,P_theta,opts,worker_ids)
    end

    owns_workspace = isnothing(ws);

    if(owns_workspace)
        ws = setup_workspace(P_theta,opts.max_constraints); # TODO larger?
    end

    settings(ws.DAQP_workspace,opts.daqp_settings)

    # Call the main certification function
    ret = certify(R0,prob,ws,opts);

    if(owns_workspace)
        DAQP.free_c_workspace(ws.DAQP_workspace);
    end
    return ret
end
function certify(S::Vector{Region}, prob::CertProblem, ws::CertWorkspace,opts::CertSettings)
    j = 0;
    while(!isempty(S) && ws.N_fin < opts.output_limit )
        # Run all callbacks
        terminate_cb = false
        for (condition,callback) in opts.conditioned_callbacks
            if(condition(S,prob,ws,opts))
                terminate_cb |= callback(S,prob,ws,opts)
            end
        end
        terminate_cb && return ws.F, ws.iter_max, ws.N_fin, ws.ASs, ws.bin, ws.ASs_state

        # Pop and process region
        region = pop!(S);
        if opts.verbose >= 2
            print("\r>> #$(j+=1)|Stack: $(length(S))|Fin: $(ws.N_fin)|Max: $(ws.iter_max)|   ");
        end
        parametric_AS_iteration(prob,region,opts,ws,S);
    end # Stack is empty

    if(!isempty(S) && ws.N_fin >= opts.output_limit) # Output limit reached
        return opts.overflow_handle(S,prob,ws,opts)
    end

    print_final_information(ws.iter_max,ws.lp_count,ws.N_fin,size(ws.ASs,2),opts)
    return ws.F, ws.iter_max, ws.N_fin, ws.ASs, ws.bin, ws.ASs_state
    #return ws.N_fin, ws.iter_max
end

function certify(region::Region, prob::CertProblem, ws::CertWorkspace,opts::CertSettings)
    ws.N_fin=0; ws.iter_max=0;
    S =[region];
    return certify(S,prob,ws,opts)
end

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

function materialize_regions!(regions::Vector{Region}, ws::CertWorkspace)
    for region in regions
        region.Ath = [ws.Ath[:,1:region.start_ind] region.Ath];
        region.bth = [ws.bth[1:region.start_ind]; region.bth];
        region.start_ind = 0;
    end
end

function print_final_information(iter_max::Int64, lp_count::Int64, N_fin::Int64, n_ASs::Int64, opts::CertSettings)
    if(opts.verbose>=1)
        ASs_string = (opts.store_ASs) ? "|ASs: $(n_ASs)" : ""
        println("\n======= Final information: =======");
        println("||Max:$(iter_max)|LPs: $(lp_count)"*ASs_string*"|Fin:$(N_fin)||");
        println("==================================");
    end
end

function rebuild_stored_ASs(part::Vector{Region}, prob::CertProblem, opts::CertSettings)
    opts.store_ASs || return falses(0,0), State[]
    n_constr = length(prob.bounds_table)
    ASs = falses(n_constr,0)
    ASs_state = State[]
    for region in part
        AS_bool = falses(n_constr)
        AS_bool[region.AS] .= true
        n_prev = size(ASs,2)
        ASs = update_ASs(ASs,AS_bool)
        if(size(ASs,2) > n_prev)
            push!(ASs_state,region.state)
        end
    end
    return ASs, ASs_state
end

function pop_seed_region!(pending::Vector{Region})
    best_ind = 1
    best_iter = pending[1].iter
    for i in 2:length(pending)
        if(pending[i].iter < best_iter)
            best_ind = i
            best_iter = pending[i].iter
        end
    end
    region = pending[best_ind]
    deleteat!(pending,best_ind)
    return region
end

function seed_regions(R0::Region, prob::CertProblem, P_theta, opts::CertSettings, target_count::Int64)
    ws = setup_workspace(P_theta,opts.max_constraints)
    settings(ws.DAQP_workspace,opts.daqp_settings)
    pending = Region[R0]
    part = Region[]
    iter_max = 0
    N_fin = 0
    lp_count = 0
    try
        while(length(pending) < target_count && !isempty(pending))
            region = pop_seed_region!(pending)
            reset_workspace(ws)
            S = Region[]
            parametric_AS_iteration(prob,region,opts,ws,S)
            materialize_regions!(S,ws)
            append!(pending,S)
            append!(part,ws.F)
            iter_max = max(iter_max,ws.iter_max)
            N_fin += ws.N_fin
            lp_count += ws.lp_count
        end
    finally
        DAQP.free_c_workspace(ws.DAQP_workspace)
    end
    return pending, part, iter_max, N_fin, lp_count
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
        part, iter_max, N_fin, _, _, _ = certify(region,prob,ws,opts)
        return part, iter_max, N_fin, ws.lp_count
    finally
        DAQP.free_c_workspace(ws.DAQP_workspace)
    end
end

function distributed_certify(R0::Region, prob::CertProblem, P_theta, opts::CertSettings, worker_ids::Vector{Int})
    seed_target = max(length(worker_ids), opts.distributed_region_factor*length(worker_ids))
    pending, part, iter_max, N_fin, lp_count = seed_regions(R0,prob,P_theta,opts,seed_target)
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
    results = Distributed.pmap(region -> distributed_region_certify(prob,P_theta,region,worker_opts), pool, ordered_regions)

    for (worker_part, worker_iter_max, worker_N_fin, worker_lp_count) in results
        append!(part,worker_part)
        iter_max = max(iter_max,worker_iter_max)
        N_fin += worker_N_fin
        lp_count += worker_lp_count
    end

    ASs, ASs_state = rebuild_stored_ASs(part,prob,opts)
    print_final_information(iter_max,lp_count,N_fin,size(ASs,2),opts)
    return part, iter_max, N_fin, ASs, Dict{Any,Any}(), ASs_state
end
## Update optimization model
function update_feas_model(reg::Region,ws::CertWorkspace)
    m = size(reg.Ath,2);
    for i in 1:m
        ws.Ath[:,reg.start_ind+i]=reg.Ath[:,i];
        ws.bth[reg.start_ind+i]=reg.bth[i];
    end
    reg.start_ind += m;
    ws.m = reg.start_ind;
end

## Prune candidates
# Determine candidates in M that yields a nonempty region
# when intersected with the reg.A th<= reg.b.
# The index of such feasible rows are retained in the vector cands
function prune_candidates(M::Matrix{Float64},ws::CertWorkspace,eps::Float64,eps_gap::Float64,cands::Vector{Int64},pos_cands::Vector{Int64})
    for i in size(M,2):-1:1
        ws.Ath[:,ws.m+1] = M[1:end-1,i];
        ws.bth[ws.m+1] = -M[end,i]-eps;
        # Normalize
        k=normalize_halfplane!(ws.Ath,ws.bth,ws.m+1;rhs_offset=eps_gap)-ws.m
        if(k<0)
            push!(pos_cands,cands[i]);
            deleteat!(cands,i);
        elseif(k>0 && ~isfeasible(ws.DAQP_workspace, ws.m+k, 0))
            push!(pos_cands,cands[i]);
            deleteat!(cands,i);
        end
    end
end

## Setup workspace
function setup_workspace(P_theta,m_max)::CertWorkspace
    nth = length(P_theta.ub);
    A = Array{Float64}(undef,nth,m_max);
    b = Array{Float64}(undef,m_max);
    blower = fill(-1e30,m_max);

    max_radius =  nth*(maximum(P_theta.ub)^2)/2; # The region is contained in a ball with this radius.
    # Create C workspace
    p=DAQP.setup_c_workspace(nth);
    ws = CertWorkspace(A,b,blower,zeros(Cint,m_max),Region[],0,0,0,
                       p,max_radius,falses(0,0),State[],0,Dict());

    DAQP.init_c_workspace_ldp(p,ws.Ath,ws.bth,ws.bth_lower,ws.sense_feasibility;max_radius)

    d_work = unsafe_load(Ptr{DAQP.Workspace}(p));
    unsafe_store!(Ptr{Cdouble}(d_work.settings+fieldoffset(DAQP.DAQPSettings,3)),1e-11);
    unsafe_store!(Ptr{Cdouble}(d_work.settings+fieldoffset(DAQP.DAQPSettings,4)),1e-6);

    return ws
end
## Reset workspace
function reset_workspace(ws::CertWorkspace)
    ws.F = Region[];
    ws.iter_max = 0;
    ws.N_fin = 0;
    ws.m=0;
    ws.ASs=falses(0,0)
    ws.lp_count=0;
end
## Get final region
function extract_regions(region::Region, ws::CertWorkspace;minrep_regions=false)
    if(minrep_regions)
        region.Ath,region.bth = minrep([ws.Ath[:,1:region.start_ind] region.Ath],
                                       [ws.bth[1:region.start_ind]; region.bth])
    else
        region.Ath= [ws.Ath[:,1:region.start_ind] region.Ath];
        region.bth= [ws.bth[1:region.start_ind]; region.bth];
    end
end
## Terminate a region
function terminate(region::Region,ws::CertWorkspace,opts::CertSettings,storage_level::Int8)

    # Run termination_callbacks
    for callback in opts.termination_callbacks
        callback(region,ws);
    end

    ws.iter_max = max(ws.iter_max,region.iter);
    ws.N_fin +=1;

    if(opts.store_ASs)
        AS_bool = falses(size(region.ASs,1)); AS_bool[region.AS] .=true;
        ws.ASs = update_ASs(ws.ASs,AS_bool);
        push!(ws.ASs_state,region.state)
    end
    storage_level==0 && return # Store nothing

    if storage_level < 2
        region.Lam, region.L, region.D  = zeros(0,0),zeros(0,0),zeros(0)
    end

    opts.store_regions && extract_regions(region,ws;minrep_regions=opts.minrep_regions)

    if(opts.compute_chebyball)
        if(opts.store_regions)
            c,r = center(region.Ath,region.bth)
        else
            c,r = center([ws.Ath[:,1:region.start_ind] region.Ath],
                         [ws.bth[1:region.start_ind]; region.bth])
        end
        r < opts.eps_cheby && return ## Do not store low-dimensional regions
        region.chebyball = (c,r)
    end

    push!(ws.F,region)
end

## Normalize half-plane (if the norm is to small, reduce the counter k)
function normalize_halfplane!(A,b,k; zero_tol = 1e-14, rhs_offset=0)
    ak = view(A, :, k)
    norm_factor = norm(ak);
    if(norm_factor < zero_tol)
        return (b[k]<0) ? -1 : k-1
    else
        inorm = inv(norm_factor)
        ak.*=inorm
        b[k]= inorm*b[k]-rhs_offset
        return k
    end
end
