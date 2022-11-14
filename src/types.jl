
# Multi-parametric quadratic program
mutable struct MPQP
    H::Matrix{Float64}
    f::Matrix{Float64}
    f_theta::Matrix{Float64}
    H_theta::Matrix{Float64}

    A::Matrix{Float64}
    b::Matrix{Float64}
    W::Matrix{Float64}

    bounds_table::Vector{Int64}
    senses::Vector{Cint}
    MPQP()=new()
    MPQP(H,f,f_theta,H_theta,A,b,W,bounds_table,senses) =new(H,f,f_theta,H_theta,A,b,W,bounds_table,senses) 
end

abstract type CertProblem{T<:Real} end
abstract type AbstractRegion end
# DualCertProblem 
struct DualCertProblem{T} <: CertProblem{T}
    MM::Matrix{T}
    M::Matrix{T}
    d::Matrix{Float64}
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
    senses::Vector{Cint}
end

# DualLPCertProblem 
struct DualLPCertProblem{T} <: CertProblem{T}
    f::Vector{T}
    A::Matrix{T}
    b::Matrix{Float64}
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
    senses::Vector{Cint}
end

function DualCertProblem(mpQP,bounds_table::Vector{Int64};normalize=true)
    R = cholesky((mpQP.H+mpQP.H')/2)
    M = mpQP.A/R.U
    V = (R.L)\[mpQP.f_theta mpQP.f]
    d = [mpQP.W mpQP.b] + M*V
    d = d'[:,:]# Col major...
    if(normalize)
        # Normalize
        norm_factor = 0
        for i in 1:size(M,1)
            norm_factor = norm(M[i,:],2) 
            M[i,:]./=norm_factor
            d[:,i]./=norm_factor
        end
    end

    MM = M*M'
    n_theta = size(d,1)-1
    if(isempty(bounds_table))
        bounds_table = 1:length(mpQP.b)
    end
    return DualCertProblem(MM,M, d, n_theta, length(mpQP.f), bounds_table, mpQP.senses)
end
function DualCertProblem(mpQP;normalize=true)
    return DualCertProblem(mpQP,mpQP.bounds_table;normalize)
end

mutable struct Region <:AbstractRegion
    IS::BitVector
    AS::Vector{Int64}
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    state::State
    iter::Int64
    start_ind::Int64
    add_ind:: Int64
    reuse_ind::Int64
    Lam::Matrix{Float64}
    ASs::BitMatrix
    L::Matrix{Real}
    D::Vector{Real}
    feas_cons::BitVector
    kappa::Dict{Symbol,Any}
end
function Region(AS::Vector{Int64},A::Matrix{Float64},b::Vector{Float64},prob::DualCertProblem)
    n_constr=size(prob.M,1)
    nth= size(A,1)
    # Create initial L and D
    L = zeros(eltype(prob.M),0,0)
    D = zeros(eltype(prob.M),0)
    for (k,ind) in enumerate(AS)
        m = prob.M[ind,:]
        L,D=DAQP.updateLDLadd(L,D,prob.M[AS[1:k-1],:]*m,m'*m)
    end
    IS = trues(n_constr)
    IS[AS] .= false
    # Create initial lambda (constant)
    Lam = zeros(nth+1,length(AS))
    Lam[end,:] .=1
    return Region(IS,AS,A,b,REMOVE,1,0,0,0,Lam,falses(n_constr,0),L,D,falses(n_constr),Dict())
end


Base.@kwdef mutable struct CertSettings 
    eps_primal::Float64 = 1e-6
    eps_dual::Float64 = 0
    eps_zero::Float64 = 1e-10
    eps_gap::Float64  = 1e-6
    verbose::Int8 = 2
    iter_limit::Int64  = 1e3
    storage_level::Int8 = 1
    compute_flops::Bool = false
    max_constraints::Int64 = 1000
    delta_lam::Float64 = 0
    delta_mu::Float64 = 0
    delta_alpha::Float64 =0
    rm_callbacks::Vector{Function} = Function[] 
    add_callbacks::Vector{Function} = Function[]
    termination_callbacks::Vector{Function} = Function[]
    pop_callbacks::Vector{Function} = Function[]
end

mutable struct CertWorkspace
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    bth_lower::Vector{Float64}
    sense_feasibility::Vector{Cint}
    F::Vector{Region}
    iter_max::Int64
    N_fin::Int64
    m::Int64
    DAQP_workspace::Ptr{Cvoid}
    max_radius::Float64
    ASs::BitMatrix
    lp_count::Int64
    bin::Dict{Any,Any}
end


struct QP
    H::Matrix{Float64}
    f::Vector{Float64}
    A::Matrix{Float64}
    b::Vector{Float64}
    senses::Vector{Cint}
end

function QP(mpQP, theta::Vector{Float64})
    return QP(mpQP.H,mpQP.f[:,1]+mpQP.f_theta*theta,mpQP.A,mpQP.b[:,1]+mpQP.W*theta,mpQP.senses) # TODO, refine f and b as vectors in LMPC...
end
