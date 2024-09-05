## Parametric AS iteration 
function parametric_AS_iteration(prob::DualCertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,S::Vector{Region})
    for callback in opts.pop_callbacks
        # if callback returns true => skip iteration
        (callback(region,prob,ws,opts)) && return
    end
    # Check if iteration limit is reached
    if(region.iter >= opts.iter_limit)
        region.state=ITERLIM
        terminate(region,ws,opts,opts.storage_level) 
        return
    end

    # Perform usual iteration
    if(region.state==ADD)
        cert_add_constraint(prob,region,opts,ws,S);
    else
        cert_remove_constraint(prob,region,opts,ws,S);
    end
end

## Add constraint
function cert_add_constraint(prob::CertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,partition::Vector{Region})
    npart0 = length(partition)
    # No indices can be added -> Global optimum
    cands = region.IS[:];
    cands[prob.bounds_table[region.AS]].=false;# Both bounds cannot be active 
    cands[region.feas_cons].=false;# These constraints cannot be primal infeasible 
    ind_cands = findall(cands); 
    #ind_cands = setdiff(ind_cands,region.feas_cons);

    if(isempty(ind_cands))
        if(opts.prune_subsequences && region.n_siblings > 0)
            return nothing
        end
        region.state=OPTIMAL;
        terminate(region,ws,opts,opts.storage_level);
        return nothing 
    end

    # Update feasibility model 
    update_feas_model(region,ws)

    # Compute μ
    μ = compute_slack(region,prob,ind_cands)

    # Update flop count
    haskey(region.kappa,:flops) && flops_compute_slack(region,prob.n)

    # Execute callbacks
    for callback in opts.add_callbacks
        callback(region,ws,opts,μ);
    end
    # Set primal tolerance
    ϵp = opts.eps_primal;
    # Only i s.t. μᵢ(θ) < -ϵ can be added
    negative_cands = collect(1:length(ind_cands));
    pos_cands = Int64[]; 
    prune_candidates(μ,ws,ϵp-opts.delta_mu,opts.eps_gap,negative_cands,pos_cands);
    region.feas_cons[ind_cands[pos_cands]].=1;

    # Θʲ=∅ ∀j ⟹ Θ* = Θ  
    if(isempty(negative_cands))
        if(opts.prune_subsequences && region.n_siblings > 0)
            return nothing
        end
        region.state=OPTIMAL;
        region.Ath= zeros(prob.n_theta,0); 
        region.bth= zeros(0); 
        region.ASs=[region.ASs falses(size(region.ASs,1))];
        region.ASs[region.AS,end].=true;
        terminate(region,ws,opts,opts.storage_level);
        return nothing;
    end

    # Compute added regions
    n_cands = length(negative_cands);
    Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+n_cands)];
    bth_tmp = @view ws.bth[(ws.m+1):(ws.m+n_cands)];
    for i in negative_cands
        k=0;
        for j in negative_cands
            if(i==j)
                # μᵢ(θ) < -ϵ 
                Ath_tmp[:,k+1] .= μ[1:end-1,i];
                bth_tmp[k+1] = -μ[end,i]-ϵp
                bth_tmp[k+1] += opts.delta_mu; # Take into account round-off errors 
            else
                # μᵢ(θ) < μⱼ(θ)
                Ath_tmp[:,k+1] = μ[1:end-1,i]-μ[1:end-1,j];
                bth_tmp[k+1] = -μ[end,i]+μ[end,j];
                bth_tmp[k+1] += 2*opts.delta_mu  # Take into account round-off errors
            end
            # Normalize
            k=normalize_halfplane!(Ath_tmp,bth_tmp,k+1;rhs_offset=opts.eps_gap)
            (k==-1)&& break; # Trivially infeasible
        end
        # Check if Θᵢ ≂̸ ∅
        if(k>=0 && isfeasible(ws.DAQP_workspace, ws.m+k, 0))
            new_region=spawn_region(region,ind_cands[i],Ath_tmp[:,1:k],bth_tmp[1:k],Float64[],prob);
            push!(partition,new_region);
        end
    end

    if(opts.prune_subsequences && length(partition) > npart0)   
        return nothing
    end
    # Check if there are any parmameters in Θ which leads to μ(θ) ≥ - ϵ
    k=0;
    for neg_ind in negative_cands #for j ∈ AS  : j≂̸i
        ws.Ath[:,ws.m+k+1] = -μ[1:end-1,neg_ind];
        ws.bth[ws.m+k+1] = μ[end,neg_ind]+ ϵp
        ws.bth[ws.m+k+1] += opts.delta_mu; # Take into account round-off errors

        #Normalize
        k=normalize_halfplane!(ws.Ath,ws.bth,ws.m+k+1;rhs_offset=opts.eps_gap)-ws.m
        (k<=-1)&& return nothing; # Trivially infeasible
    end
    if(isfeasible(ws.DAQP_workspace, ws.m + k, 0))
        ws.m+=k #Update global model (since last)  
        region.state=OPTIMAL;
        region.start_ind =ws.m;
        region.ASs=[region.ASs falses(size(region.ASs,1))];
        region.ASs[region.AS,end].=true;
        terminate(region,ws,opts,opts.storage_level);
    end

    return nothing 
end

## Remove constraint
function cert_remove_constraint(prob::DualCertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,partition::Vector{Region})
    npart0 = length(partition)
    # Empty AS => Trivial CSP 
    if(isempty(region.AS)) 
        region.state=ADD;
        region.Lam= zeros(prob.n_theta+1,0);
        region.n_siblings=0;
        push!(partition,region);
        return nothing 
    end

    # Append constraints 
    update_feas_model(region,ws)

    n_active = length(region.AS);

    # TODO Singulariy only needs to be checked directly after addition...
    # TODO Pivoting...
    # Check singularity
    singular_ind = findfirst(region.D.<1e-14)
    singular = !isnothing(singular_ind);

    # Compute λᶜ (Constrained stationary point)
    if(!singular)
        λᶜ = -prob.d[:,region.AS]
        forward_L_para!(region.L,λᶜ);
        for j in 1:length(region.D) 
            λᶜ[:,j] ./=region.D[j];
        end
        backward_L_para!(region.L,λᶜ);
        #Update flops count
        haskey(region.kappa,:flops) && flops_solve_kkt(region)
    else
        #Update flops count
        haskey(region.kappa,:flops) && flops_singular_direction(region)
    end

    # Execute callbacks
    for callback in opts.rm_callbacks
        callback(region,ws,opts,λᶜ);
    end

    # === Partition Θ depending on constraints that are removed  === 
    if(region.add_ind == 0) # No constraint has been added to AS yet  (Case b2)
        # TODO  can we use region.feas_cons?
        #feas_inds=collect(1:n_active);
        feas_inds = Int64[]; # Only inequality constraints can be removed
        for (k,ind) in enumerate(region.AS)
            if(prob.senses[ind]&DAQP.IMMUTABLE==0)
                push!(feas_inds,k);
            end
        end
        # Only keep i such that λ*ᵢ(θ) < -ϵ for some θ∈Θ 
        prune_candidates(λᶜ[:,feas_inds],ws,opts.eps_dual,opts.eps_gap,feas_inds,Int64[]);
        # Candidates... 
        Δλ = -λᶜ+region.Lam;
        Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+length(feas_inds))];
        bth_tmp = @view ws.bth[(ws.m+1):(ws.m+length(feas_inds))];

        for (k_outer,i) in enumerate(feas_inds)
            k=0
            for j in filter(x->x!=i,feas_inds) #for j ∈ AS  : j≂̸i
                # α_i(θ) < α_j(θ)
                # TODO ===== Make sure that H*_k H cancels in the dual case.. =======  
                Ath_tmp[:,k+1] .= -region.Lam[end,j].*Δλ[1:end-1,i].+region.Lam[end,i].*Δλ[1:end-1,j];
                bth_tmp[k+1] = region.Lam[end,j]*Δλ[end,i]-region.Lam[end,i]*Δλ[end,j];
                # Normalize
                k=normalize_halfplane!(Ath_tmp,bth_tmp,k+1;rhs_offset=opts.eps_gap)
                (k==-1)&& break; # trivially infeasible
            end
            (k==-1) && continue # trivially infeasible 

            # λ*ᵢ(θ) < -ϵ
            Ath_tmp[:,k+1] .= λᶜ[1:end-1,i];
            bth_tmp[k+1] = -λᶜ[end,i]-opts.eps_dual;
            # Normalize
            norm_factor = norm(Ath_tmp[:,k+1],2);
            k=normalize_halfplane!(Ath_tmp,bth_tmp,k+1;rhs_offset=opts.eps_gap)
            (k==-1)&& continue; # trivially infeasible

            # Check if Θᵢ≂̸ ∅
            if(isfeasible(ws.DAQP_workspace, ws.m+k, 0))
                new_region = spawn_region(region,-i,Ath_tmp[:,1:k],bth_tmp[1:k],Float64[],prob);
                push!(partition,new_region);
            end
        end

    elseif(region.add_ind > 0) # A constraint had been added to AS (Case b1) 
        # Compute λ* and/or p̂ 
        if(!singular)
            p̂= -DAQP.compute_singulardirection(region.L,region.D,length(region.D))
        else
            p̂= -DAQP.compute_singulardirection(region.L,region.D,singular_ind)
        end

        # Update feasible candidates
        region.feas_cons[region.feas_cons] .= prob.MM[region.feas_cons,region.AS]*p̂.<0;
        feas_inds = Int64[]
        for i in 1:(n_active-1) # Never check last because p̂[end] < 0 after addition...
            if((prob.senses[region.AS[i]]&DAQP.IMMUTABLE == 0) && p̂[i]>opts.eps_zero)
                push!(feas_inds,i);
            end
        end

        if(singular)
            if(isempty(feas_inds)) # Primal infeasible problem 
                region.state=INFEASIBLE;
                region.Ath= zeros(prob.n_theta,0); 
                region.bth= zeros(0); 
                terminate(region,ws,opts,opts.storage_level);
                println("Infeasible region found")
                return nothing 
            else
                Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+length(feas_inds)-1)];
                bth_tmp = @view ws.bth[(ws.m+1):(ws.m+length(feas_inds)-1)];
            end
        else # Nonsingular reduced Hessian
            # Only keep i such that λ*ᵢ(θ) < -ϵ for some θ∈Θ 
            # TODO: eps_gap becomes normalized here...
            prune_candidates(λᶜ[:,feas_inds],ws,opts.eps_dual,opts.eps_gap,feas_inds,Int64[]);
            Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+length(feas_inds))];
            bth_tmp = @view ws.bth[(ws.m+1):(ws.m+length(feas_inds))];
            if(isempty(feas_inds))
                # Θʲ=∅ ∀j ⟹ Θ^CSP = Θ  
                region.Lam = λᶜ;
                region.Ath= zeros(prob.n_theta,0); 
                region.bth= zeros(0); 
                region.state=ADD;
                region.start_ind = ws.m;
                region.n_siblings=0;
                push!(partition,region);
                return nothing 
            end
        end

        # Partition Θ based on removals from AS
        for i in feas_inds
            k=0
            for j in (filter(x->x!=i,feas_inds)) #for j ∈ AS  : j≂̸i
                # α_i(θ) < α_j(θ)
                Ath_tmp[:,k+1] .= -p̂[i].*region.Lam[1:end-1,j].+p̂[j].*region.Lam[1:end-1,i];
                bth_tmp[k+1] = p̂[i]*region.Lam[end,j]-p̂[j]*region.Lam[end,i];
                # Normalize
                k=normalize_halfplane!(Ath_tmp,bth_tmp,k+1;rhs_offset=opts.eps_gap)
                (k==-1)&& break; # trivially infeasible
            end
            (k==-1) && continue #trivially infeasible

            if(!singular)
                # λ*ᵢ(θ) < -ϵ
                Ath_tmp[:,k+1] .= λᶜ[1:end-1,i];
                bth_tmp[k+1] = -λᶜ[end,i]-opts.eps_dual;
                # Normalize
                k=normalize_halfplane!(Ath_tmp,bth_tmp,k+1;rhs_offset=opts.eps_gap)
                (k==-1)&& continue; # trivially infeasible
            end

            # Check if Θᵢ≂̸ ∅
            if(isfeasible(ws.DAQP_workspace, ws.m+k,0))
                new_region = spawn_region(region,-i,Ath_tmp[:,1:k],bth_tmp[1:k],p̂,prob);
                push!(partition,new_region);
            end
        end
    end
    # Create Θ_CSP, i.e., θ∈Θ s.t. λ*(θ)≥-ϵ
    if(singular)  
        return nothing# CSP not defined for singular case
    end
    k=0
    for feas_ind in feas_inds #for j ∈ AS  : j≂̸i
        ws.Ath[:,ws.m+k+1] = -λᶜ[1:end-1,feas_ind];
        ws.bth[ws.m+k+1] = λᶜ[end,feas_ind]+opts.eps_dual;
        # Normalize
        k=normalize_halfplane!(ws.Ath,ws.bth,ws.m+k+1;rhs_offset=opts.eps_gap)-ws.m
        (k<=-1)&& return nothing; # trivially infeasible
    end
    if(isfeasible(ws.DAQP_workspace, ws.m+k, 0))
        ws.m+=k; #update model
        region.Lam = λᶜ;
        region.Ath= zeros(prob.n_theta,0); 
        region.bth= zeros(0); 
        region.state=ADD;
        region.start_ind = ws.m;
        region.n_siblings = length(partition)-npart0
        push!(partition,region);
    end
    return nothing 
end
## Spawn region
function spawn_region(region::Region, i::Int64, Ath::Matrix{Float64}, bth::Vector{Float64}, p̂, prob::DualCertProblem)
    n_active = length(region.AS);
    new_region=Region(region.IS[:], region.AS[:],
                      Ath,bth, REMOVE,region.iter+1,
                      region.start_ind,region.add_ind,region.reuse_ind,-1,
                      Array{Float64}(undef,size(region.Lam,1),n_active+sign(i)),
                      BitArray(undef,size(region.ASs).+(0,1)),
                      Array{Float64}(undef,0,0),
                      Array{Float64}(undef,0),
                      region.feas_cons[:],
                      (zeros(0),NaN),
                      deepcopy(region.kappa));
    if(i>0) # add
        region_add_constraint(i,region,new_region,prob)
    else # remove 
        region_remove_constraint(abs(i),region,new_region,prob,p̂)
    end
    # Update parent
    new_region.ASs[:, 1:end-1] = region.ASs;
    new_region.ASs[:,end]=.~region.IS;
    return new_region
end
## Update AS & LDL for region
function region_add_constraint(add_ind,src, dest, prob)
    m = prob.M[:,add_ind];
    dest.L,dest.D=DAQP.updateLDLadd(src.L,src.D,prob.M[:,src.AS]'*m,m'*m);
    haskey(src.kappa, :flops) && flops_add_constraint(src,size(prob.M,1));
    push!(dest.AS,add_ind);
    dest.IS[add_ind] = false;
    dest.Lam[:,1:end-1].=src.Lam;
    dest.Lam[:,end].=0;
    dest.add_ind=add_ind;
    dest.reuse_ind=length(dest.AS)-1;
end

function region_remove_constraint(rm_ind,src, dest, prob, p̂)
    n_active = length(src.AS)
    dest.L,dest.D=DAQP.updateLDLremove(src.L,src.D,rm_ind);
    haskey(src.kappa, :flops) && flops_remove_constraint(rm_ind,src,n_active)
    dest.IS[src.AS[rm_ind]]=true;
    deleteat!(dest.AS,rm_ind);
    if(!isempty(p̂))
        for (k,j) in enumerate(filter(x->x!=rm_ind,1:n_active))
            dest.Lam[:,k].= src.Lam[:,j].-(p̂[j]/p̂[rm_ind]).*src.Lam[:,rm_ind];
        end
    else
        dest.Lam.= src.Lam[:,1:end.!=rm_ind];
    end
    dest.reuse_ind=rm_ind-1;
end
## Compute slack 
function compute_slack(region,prob::DualCertProblem,ind_cands)
    μ = prob.d[:,ind_cands];
    if(~isempty(region.AS))
        μ+= region.Lam*prob.MM[region.AS,ind_cands]
    end
    return μ
end

## Normalize problem to a box -1 ≤ θ ≤ 1
function normalize(prob::Union{DualCertProblem,DualLPCertProblem},P_theta,mpQP)
    # TODO also normalize TH part of P_theta
    prob = deepcopy(prob);
    nth = length(P_theta.lb);
    center = (P_theta.lb+P_theta.ub)/2;
    norm_factors = (P_theta.ub-P_theta.lb)/2;

    prob.d[end:end,:] += center'*prob.d[1:end-1,:];
    mpQP.b[:] += mpQP.W*(center)
    mpQP.f[:] += mpQP.f_theta*(center)

    for i in 1:nth
        prob.d[i,:] *= norm_factors[i];
        mpQP.f_theta[:,i]*=norm_factors[i];
        mpQP.W[:,i]*=norm_factors[i];
    end
    if(isempty(P_theta.A))
        A = []
        b = []
    else
        b = P_theta.b - P_theta.A'*center;
        A = P_theta.A
        for i = 1:nth
            A[i,:] .*= norm_factors[i]
        end
    end

    P_theta =(A=A,
              b = b,
              lb = -ones(nth),
              ub = ones(nth));
    return prob, P_theta, mpQP
end
## Slice mpQP
function reduce_mpqp(mpQP, P_theta,ids;slice_vals=[])
    nth = size(mpQP.W,2)
    slice_ids= setdiff(1:nth,ids) #ids to remove 
    if(isempty(slice_vals)) # Default slice at 0
        slice_vals = zeros(length(slice_ids))
    end

    # Reduce the mpQP
    mpQP_new = (H = copy(mpQP.H),
                f = mpQP.f+mpQP.f_theta[:,slice_ids]*slice_vals,
                A = copy(mpQP.A),
                b = mpQP.b+mpQP.W[:,slice_ids]*slice_vals,
                f_theta = mpQP.f_theta[:,ids],
                W = mpQP.W[:,ids],
                H_theta = mpQP.H_theta[ids,ids],
                senses = copy(mpQP.senses),
                bounds_table = copy(mpQP.bounds_table))

    # Reduce the region 
    b_new = isempty(P_theta.b) ? P_theta.b :  P_theta.b + slice_vals'*P_theta.A[slice_ids,:]
    A_new = P_theta.A[ids,:];

    P_theta_new = (A=A_new,
                   b=b_new,
                   lb= P_theta.lb[ids], ub = P_theta.ub[ids])  
    return mpQP_new,P_theta_new
end
