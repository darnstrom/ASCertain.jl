## Parametric AS LP iteration 
function parametric_AS_iteration(prob::DualLPCertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,S::Vector{Region})
    # Update feasibility model 
    if(region.state==INFEASIBLE||region.state==UNBOUNDED)
        terminate(region,ws,opts,opts.storage_level);
        return
    end
    region.feas_cons[:].=false # since feas_cons is only valid for QPs 
    cert_add_constraint(prob,region,opts,ws,S)
    # No parameter dependence in objective => removal does not split the parameter space 
    # Hence, we do the removal region-wise in spawn_regions
end

## Compute slack 
function compute_slack(region,prob::DualLPCertProblem,ind_cands)
    x = prob.b[:,region.AS]/(prob.A[region.AS,:]')
    μ=prob.b[:,ind_cands]-x*(prob.A[ind_cands,:])'
    return μ 
end

## Spawn region
function spawn_region(region::Region, i::Int64, Ath::Matrix{Float64}, bth::Vector{Float64}, p̂::Array{Float64}, prob::DualLPCertProblem)
    new_region=Region(region.IS[:], region.AS[:],
                      Ath,bth, REMOVE,region.iter+1,region.start_ind,region.add_ind,region.reuse_ind,
                      Array{Float64}(undef,size(region.Lam,1),prob.n),
                      BitArray(undef,size(region.ASs).+(0,1)),
                      Array{Float64}(undef,0,0),
                      Array{Float64}(undef,0),
                      region.feas_cons[:],
                      deepcopy(region.kappa));

    # pivot to find new basis
    new_region.Lam[end,:],valid_pivot = pivot_add!(new_region.AS,new_region.IS,
                                                   region.Lam[:],i,prob.A)

    if(!valid_pivot) 
        # Dual unbounded => primal_infeasible
        new_region.state=INFEASIBLE;
    end
    # Update parent
    new_region.ASs[:, 1:end-1] = region.ASs;
    new_region.ASs[:,end]=.~region.IS;

    return new_region
end

## Region constructor
function Region(AS::Vector{Int64},A::Matrix{Float64},b::Vector{Float64},prob::DualLPCertProblem)
    m,n=size(prob.A)
    nth= size(A,1)
    if(isempty(AS))
        lam,AS,iter,flag_ph1=dphase1(prob.f,prob.A)
        if (flag_ph1!=0) 
            state = UNBOUNDED
            AS = Int64[]
            lam = zeros(0)
        else
            state = REMOVE 
        end
    else # Make sure that vertex is dual feasible
        lam = -prob.A[AS,:]'\prob.f;
        iter = 0;
        if(any(lam.<0))
            state = INFEASIBLE 
        else
            state = REMOVE
        end
    end
    IS = trues(m)
    IS[AS] .= false
    Lam = zeros(1,length(AS))
    Lam[1,:]=lam;
    # Create initial lambda (constant)
    return Region(IS,AS,A,b,state,iter,0,0,0,Lam,falses(m,0),zeros(0,0),zeros(0),falses(m),Dict())
end
