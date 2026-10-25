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
    R0.crs = collect(Int32,1:length(opts.explicit))
    opts.use_region_box && (R0.box = [Float64.(P_theta.lb) Float64.(P_theta.ub)])
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
    ws.warm_start = opts.warm_start_ldp
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

## Update optimization model
function update_feas_model(reg::Region,ws::CertWorkspace)
    m = size(reg.Ath,2);
    for i in 1:m
        ws.Ath[:,reg.start_ind+i]=reg.Ath[:,i];
        ws.bth[reg.start_ind+i]=reg.bth[i];
    end
    reg.start_ind += m;
    ws.m = reg.start_ind;
    ws.n_regions += 1
    ws.warm_hint = reg.warm_WS
    update_box!(reg,ws)
end

## Bounding box of a region
# The box of a region is inherited from its parent and tightened by interval propagation over the
# constraints that have been added since (no LPs). Regions without a box use the global box.
function update_box!(reg::Region,ws::CertWorkspace)
    if(isempty(reg.box))
        ws.node_lb .= ws.theta_lb; ws.node_ub .= ws.theta_ub
        return
    end
    lb, ub = reg.box[:,1], reg.box[:,2]
    for sweep in 1:2, c in reg.box_m+1:ws.m
        tighten_box!(lb,ub,ws.Ath,ws.bth[c],c)
    end
    for l in eachindex(lb) # Keep the box nonempty (the region is nonempty up to tolerances)
        lb[l] > ub[l] && (lb[l] = ub[l] = (lb[l]+ub[l])/2)
    end
    reg.box, reg.box_m = [lb ub], ws.m
    ws.node_lb .= lb; ws.node_ub .= ub
end

# Tighten the box of the region ws.Ath[:,1:ws.m]'θ ≤ ws.bth[1:ws.m] with supporting halfspaces (2n_θ LDPs):
# the projection p of a point c outside the region gives the valid inequality (c-p)'θ ≤ (c-p)'p.
# The projections are also points in the region (witnesses).
function projection_box!(reg::Region,ws::CertWorkspace)
    nth, m = size(ws.Ath,1), ws.m
    lb, ub = reg.box[:,1], reg.box[:,2]
    b0 = ws.bth[1:m]
    c = zeros(nth); n = zeros(nth)
    for l in 1:nth, s in (-1.0,1.0)
        width = ub[l]-lb[l]
        width < 1e-9 && continue
        c .= (lb .+ ub) ./ 2
        c[l] = (s > 0 ? ub[l] : lb[l]) + s*width # Outside the box in direction s*e_l
        # Project c onto the region: θ = c + v with min ‖v‖ s.t. Ath'v ≤ bth - Ath'c
        @inbounds for j in 1:m
            acc = 0.0
            for k in 1:nth; acc += ws.Ath[k,j]*c[k]; end
            ws.bth[j] = b0[j]-acc
        end
        ws.lp_count += 1
        feas = isfeasible(ws.DAQP_workspace, m, 0)
        ws.bth[1:m] .= b0
        feas || break # (should not happen for a nonempty region)
        p = ws.daqp_u .+ c
        n .= c .- p
        nn = norm(n)
        nn < 1e-12 && continue
        n ./= nn
        tighten_box!(lb,ub,reshape(n,nth,1),dot(n,p)+1e-9,1)
        push_witness!(ws,p)
    end
    for l in 1:nth
        lb[l] > ub[l] && (lb[l] = ub[l] = (lb[l]+ub[l])/2)
    end
    reg.box = [lb ub]
    ws.node_lb .= lb; ws.node_ub .= ub
end

# Number of candidates i with a nonempty bound test for {θ : μᵢ(θ) < -eps} (as in prune_candidates and box_infeasible)
function count_box_undecided(μ::Matrix{Float64},ws::CertWorkspace,eps::Float64,eps_gap::Float64)
    tol = unsafe_load(ws.daqp_primal_tol)
    nth = size(ws.Ath,1)
    a = zeros(nth)
    n = 0
    for i in 1:size(μ,2)
        na = norm(view(μ,1:nth,i))
        na < 1e-14 && continue
        a .= view(μ,1:nth,i) ./ na
        box_infeasible(ws,a,(-μ[end,i]-eps)/na-eps_gap,tol) || (n += 1)
    end
    return n
end

# Is {θ in the bounding box of the region : a'θ ≤ b} empty? (with a margin for the primal tolerance)
function box_infeasible(ws::CertWorkspace,a::AbstractVector{Float64},b::Float64,tol::Float64)
    nth = length(a)
    min_ath, norm1 = 0.0, 0.0
    @inbounds for l in 1:nth
        al = a[l]
        min_ath += (al > 0) ? al*ws.node_lb[l] : al*ws.node_ub[l]
        norm1 += abs(al)
    end
    return min_ath-tol*norm1 > b+tol
end

# Tighten lb ≤ θ ≤ ub with A[:,c]'θ ≤ b
function tighten_box!(lb,ub,A,b,c)
    nth = length(lb)
    min_act = 0.0
    @inbounds for l in 1:nth
        a = A[l,c]
        min_act += (a > 0) ? a*lb[l] : a*ub[l]
    end
    @inbounds for l in 1:nth
        a = A[l,c]
        abs(a) < 1e-9 && continue
        rest = min_act - ((a > 0) ? a*lb[l] : a*ub[l])
        bound = (b-rest)/a
        if(a > 0)
            bound < ub[l] && (ub[l] = bound)
        else
            bound > lb[l] && (lb[l] = bound)
        end
    end
end

## Prune candidates
# Determine candidates in M that yields a nonempty region
# when intersected with the reg.A th<= reg.b.
# The index of such feasible rows are retained in the vector cands
function prune_candidates(M::Matrix{Float64},ws::CertWorkspace,eps::Float64,eps_gap::Float64,cands::Vector{Int64},pos_cands::Vector{Int64};use_witnesses=false)
    tol = unsafe_load(ws.daqp_primal_tol)
    nth = size(M,1)-1
    for i in size(M,2):-1:1
        # Bound test before normalizing (the same test as box_infeasible after normalization, scaled by the norm)
        if(use_witnesses)
            norm2, min_ath, norm1 = 0.0, 0.0, 0.0
            @inbounds for l in 1:nth
                a = M[l,i]
                norm2 += a*a; norm1 += abs(a)
                min_ath += (a > 0) ? a*ws.node_lb[l] : a*ws.node_ub[l]
            end
            na = sqrt(norm2)
            if(na >= 1e-14 && min_ath-tol*norm1 > (-M[end,i]-eps)-eps_gap*na+tol*na)
                push!(pos_cands,cands[i]);
                deleteat!(cands,i);
                continue
            end
        end
        ws.Ath[:,ws.m+1] .= view(M,1:size(M,1)-1,i);
        ws.bth[ws.m+1] = -M[end,i]-eps;
        # Normalize
        k=normalize_halfplane!(ws.Ath,ws.bth,ws.m+1;rhs_offset=eps_gap)-ws.m
        if(k<0)
            push!(pos_cands,cands[i]);
            deleteat!(cands,i);
        elseif(k>0 && find_feasible_point(ws,k,use_witnesses)==0)
            push!(pos_cands,cands[i]);
            deleteat!(cands,i);
        end
    end
end

## Feasibility check
# Check if {θ : ws.Ath[:,1:ws.m+k]'θ ≤ ws.bth[1:ws.m+k]} is nonempty.
# Returns 0 if it is empty, otherwise the column in ws.witnesses of a point in the set.
# (Without witnesses, 1 is returned if the set is nonempty)
# Before an LP is solved, known points and bounds on θ are used to decide feasibility.
function find_feasible_point(ws::CertWorkspace,k::Int64,use_witnesses::Bool)
    if(use_witnesses)
        j = check_witnesses(ws,k)
        j > 0 && return j
        box_infeasible(ws,k) && return 0
    end
    ws.warm_start && warm_start!(ws)
    ws.lp_count += 1
    isfeasible(ws.DAQP_workspace, ws.m+k, 0) || return 0
    return use_witnesses ? push_witness!(ws,ws.daqp_u) : 1
end

# Start DAQP from the optimal working set of the current region ws.Ath[:,1:ws.m]'θ ≤ ws.bth[1:ws.m]
# (the constraints that are tight at its minimum-norm point). The LDPs of a region only differ in the
# k latest constraints, so DAQP then typically only has to add those. The working set is computed by
# one LDP per region (its solution is also a witness), which in turn starts from the working set of the
# latest such LDP on the path to the region (the region only has a few more constraints).
function warm_start!(ws::CertWorkspace)
    key = (ws.n_regions,ws.m)
    if(ws.warm_key != key)
        ws.warm_key = key
        empty!(ws.warm_WS)
        ws.lp_count += 1
        activate_constraints!(ws,ws.warm_hint)
        isfeasible(ws.DAQP_workspace, ws.m, 0) || return # (only the k latest constraints are new)
        u = ws.daqp_u
        push_witness!(ws,u)
        nth = length(u)
        for j in 1:ws.m
            r = ws.bth[j]
            for l in 1:nth; r -= ws.Ath[l,j]*u[l]; end
            abs(r) < 1e-9 && push!(ws.warm_WS,j)
            length(ws.warm_WS) == nth && break
        end
        empty!(ws.warm_hint); append!(ws.warm_hint,ws.warm_WS) # (inherited by the regions spawned from this one)
    end
    activate_constraints!(ws,ws.warm_WS)
end

function activate_constraints!(ws::CertWorkspace,WS::Vector{Int64})
    isempty(WS) && return
    for j in WS
        ws.sense_feasibility[j] |= DAQP.ACTIVE
    end
    ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),ws.DAQP_workspace)
    ccall((:daqp_activate_constraints,DAQP.libdaqp),Cint,(Ptr{Cvoid},),ws.DAQP_workspace)
end

# Return index of a point in the witness pool that satisfies the k latest constraints
function check_witnesses(ws::CertWorkspace,k::Int64)
    nth = size(ws.Ath,1)
    @inbounds for j in 1:ws.n_witness
        valid = true
        for c in ws.m+1:ws.m+k
            ath = 0.0
            for l in 1:nth
                ath += ws.Ath[l,c]*ws.witnesses[l,j]
            end
            if(ath > ws.bth[c])
                valid = false
                break
            end
        end
        valid && return j
    end
    return 0
end

# Check if any of the k latest constraints is infeasible for all θ in the bounding box of the region
# (with a margin for the primal tolerance used by DAQP)
function box_infeasible(ws::CertWorkspace,k::Int64)
    tol = unsafe_load(ws.daqp_primal_tol)
    for c in ws.m+1:ws.m+k
        box_infeasible(ws,view(ws.Ath,:,c),ws.bth[c],tol) && return true
    end
    return false
end

function push_witness!(ws::CertWorkspace,θ::AbstractVector{Float64})
    n_max = size(ws.witnesses,2)
    ws.i_witness = mod1(ws.i_witness+1,n_max) # Overwrite oldest point if the pool is full
    ws.witnesses[:,ws.i_witness] .= θ
    ws.n_witness = min(ws.n_witness+1,n_max)
    return ws.i_witness
end

# Restart the pool with the point known to be in the region
function reset_witnesses!(ws::CertWorkspace,region::Region)
    ws.n_witness, ws.i_witness = 0, 0
    isempty(region.witness) || push_witness!(ws,region.witness)
end

witness_point(ws::CertWorkspace,j::Int64,use_witnesses::Bool) = use_witnesses ? ws.witnesses[:,j] : zeros(0)

## Setup workspace
function setup_workspace(P_theta,m_max)::CertWorkspace
    nth = length(P_theta.ub);
    A = Array{Float64}(undef,nth,m_max);
    b = Array{Float64}(undef,m_max);
    blower = fill(-1e30,m_max);

    max_radius =  nth*(maximum(P_theta.ub)^2)/2; # The region is contained in a ball with this radius.
    # Create C workspace
    p=DAQP.setup_c_workspace(nth);
    d_work = unsafe_load(Ptr{DAQP.Workspace}(p));
    daqp_u = unsafe_wrap(Vector{Cdouble}, d_work.u, nth, own=false)
    daqp_primal_tol = Ptr{Cdouble}(d_work.settings+fieldoffset(DAQP.DAQPSettings,1))
    ws = CertWorkspace(A,b,blower,zeros(Cint,m_max),Region[],0,0,0,
                       p,max_radius,falses(0,0),State[],0,Dict(),
                       Float64.(P_theta.lb),Float64.(P_theta.ub),Float64.(P_theta.lb),Float64.(P_theta.ub),zeros(nth,32),0,0,
                       daqp_u,daqp_primal_tol,false,0,(-1,-1),Int64[],Int64[]);

    DAQP.init_c_workspace_ldp(p,ws.Ath,ws.bth,ws.bth_lower,ws.sense_feasibility;max_radius)

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
        opts.store_regions || materialize_region!(region,ws)
        c,r = center(region.Ath,region.bth)
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
