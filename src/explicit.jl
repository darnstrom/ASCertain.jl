## Pruning based on the explicit solution
# At every ADD state, the iterate u of the dual active-set method is the minimum-norm point of a
# relaxation that contains the optimal solution u*. Hence, u lies in the ball with diameter [0,u*],
# and constraint i can only be added at θ if d_i(θ) + s*_i(θ) < ‖u*(θ)‖, where s*_i is the optimal slack.
# Note: this is only valid for θ in the feasible parameter set (where u* exists).

# min c'[θ;1] s.t. Ath'θ ≤ bth (-Inf if the LP is not solved, which is conservative below)
function lp_min(c::AbstractVector{Float64}, Ath::Matrix{Float64}, bth::Vector{Float64})
    m = length(bth)
    x,_,exitflag,_ = DAQP.linprog(c[1:end-1],Ath,bth,fill(-1e30,m),zeros(Cint,m);A_rowmaj=true)
    return exitflag == 1 ? dot(c[1:end-1],x)+c[end] : -Inf
end

# Constraints that are never added for any θ in the given regions.
# Each region is a tuple (Ath,bth,U) with Ath'θ ≤ bth and u*(θ) = U*[θ;1], in the coordinates of prob.
# Only LPs are solved (no vertex enumeration): ‖u*‖ is bounded through the range of each of its components
function never_added_constraints(prob::DualCertProblem, regions; tol=1e-6)
    never = trues(size(prob.M,2))
    for (Ath,bth,U) in regions
        R2 = 0.0
        for k in 1:size(U,1)
            lo, hi = lp_min(U[k,:],Ath,bth), -lp_min(-U[k,:],Ath,bth)
            R2 += max(lo^2,hi^2)
        end
        R = sqrt(R2)
        for i in findall(never)
            g = 2*prob.d[:,i] - (prob.M[:,i]'*U)' # d_i + s*_i
            lp_min(g,Ath,bth) >= R + tol || (never[i] = false)
        end
    end
    return findall(never)
end

# Same as above, but with regions given in the original θ coordinates of mpQP as named tuples
# (Ath, bth, z) with Ath'θ ≤ bth and optimal solution x*(θ) = z'*[θ;1] (e.g., critical regions from ParametricDAQP).
function never_added_constraints(mpQP, P_theta, regions; tol=1e-6)
    prob, nregions = normalize_regions(mpQP, P_theta, regions)
    return never_added_constraints(prob,nregions;tol)
end

# Map regions (Ath, bth, z) in the original θ coordinates to tuples (Ath, bth, U) in the coordinates used by certify
function normalize_regions(mpQP, P_theta, regions)
    prob = setup_certproblem(mpQP)
    prob,_,_ = normalize(prob,deepcopy(P_theta),deepcopy(mpQP))
    center, scale = (P_theta.lb+P_theta.ub)/2, (P_theta.ub-P_theta.lb)/2
    nregions = map(regions) do r
        U = prob.R.U*r.z' + prob.V # u*(θ) = U*[θ;1] (original θ)
        Un = [U[:,1:end-1].*scale' U[:,end]+U[:,1:end-1]*center]
        (r.Ath.*scale, r.bth-r.Ath'*center, Un)
    end
    return prob, nregions
end

## Pruning per region in the certification
# Critical regions with the constraints that are never added in each of them, to be passed as
# CertSettings(explicit=...). At an ADD state, a candidate is removed if it is never added in any of the
# critical regions that intersect the current region. (Only valid if Θ ⊆ feasible parameter set.)
function explicit_regions(mpQP, P_theta, regions; tol=1e-6)
    prob, nregions = normalize_regions(mpQP, P_theta, regions)
    m = size(prob.M,2)
    return map(nregions) do r
        never = falses(m)
        never[never_added_constraints(prob,[r];tol)] .= true
        Ath, bth = copy(r[1]), copy(r[2])
        for j in 1:length(bth) # Normalize (as the constraints in the certification)
            nj = norm(Ath[:,j]); Ath[:,j] ./= nj; bth[j] /= nj
        end
        lb = [lp_min([Float64.(1:size(Ath,1).==l);0.0],Ath,bth) for l in 1:size(Ath,1)]
        ub = [-lp_min([-Float64.(1:size(Ath,1).==l);0.0],Ath,bth) for l in 1:size(Ath,1)]
        ExplicitRegion(Ath,bth,never,lb,ub)
    end
end

# Remove candidates in ind_cands that are never added in any of the critical regions that might intersect
# the region (given by ws.Ath[:,1:ws.m]'θ ≤ ws.bth[1:ws.m]). No LPs are solved: region.crs is refined by
# removing critical regions whose bounding box violates one of the constraints of the region.
function explicit_prune!(ind_cands::Vector{Int64},region::Region,ws::CertWorkspace,opts::CertSettings;tol=1e-7)
    nth = size(ws.Ath,1)
    keep = trues(length(region.crs))
    for (j,k) in enumerate(region.crs)
        cr = opts.explicit[k]
        for c in region.crs_m+1:ws.m
            min_ath = 0.0
            for l in 1:nth
                a = ws.Ath[l,c]
                min_ath += (a > 0) ? a*cr.lb[l] : a*cr.ub[l]
            end
            if min_ath > ws.bth[c] + tol
                keep[j] = false
                break
            end
        end
    end
    all(keep) || (region.crs = region.crs[keep])
    region.crs_m = ws.m
    isempty(region.crs) && return ind_cands # Outside all critical regions => no pruning
    filter!(i -> !all(k -> opts.explicit[k].never[i], region.crs), ind_cands)
    return ind_cands
end

