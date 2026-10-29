module ASCertain

using DAQPBase, Distributed, LinearAlgebra, PolyDAQP
const DAQP = DAQPBase

export certify, merged_certify,
       Region,
       DualCertProblem,DualLPCertProblem,
       MPQP,QP,
       CertSettings
@enum State REMOVE ADD OPTIMAL INFEASIBLE ITERLIM UNBOUNDED
include("types.jl");
include("ascert.jl");
include("distribute.jl");


include("dualQPCert.jl");
include("dualLPCert.jl");

include("flops.jl")
include("utils.jl");
include("plot.jl");

# Implementation of solvers
include("solvers/dsimplex.jl");

# Compile the certification during precompilation (also for the workers in distributed certification)
using PrecompileTools: @setup_workload, @compile_workload
@setup_workload begin
    # Small MPC problem: double integrator, horizon 3, |u| ≤ 1, |x| ≤ 2, θ = x0 ∈ [-1,1]² (feasible)
    Ad, Bd, N = [1.0 1; 0 1], [0.5, 1.0], 3
    Φ, Γ = zeros(2N,2), zeros(2N,N)
    for k in 1:N
        Φ[2k-1:2k,:] = Ad^k
        for j in 1:k; Γ[2k-1:2k,j] = Ad^(k-j)*Bd; end
    end
    G, E = [Matrix(1.0I,N,N); Γ], [zeros(N,2); Φ]
    m = size(G,1)
    bu = [ones(N); 2*ones(2N)]
    mpQP = MPQP(Γ'Γ+I, zeros(N,1), Γ'Φ, zeros(2,2), [G;-G], reshape([bu;bu],:,1), [-E;E],
                [collect(m+1:2m);collect(1:m)], zeros(Cint,2m))
    P_theta = (A=zeros(2,0), b=zeros(0), lb=-ones(2), ub=ones(2))
    @compile_workload begin
        certify(mpQP,P_theta;opts=CertSettings(verbose=0))
        certify(mpQP,P_theta;opts=CertSettings(verbose=0,storage_level=0))
        # What the workers run (see distribute.jl)
        prob,Pn,_ = normalize(setup_certproblem(mpQP),deepcopy(P_theta),deepcopy(mpQP))
        R0 = Region(Int64[],[Matrix(1.0I,2,2) Matrix(-1.0I,2,2)],[Pn.ub;-Pn.lb],prob)
        R0.box = [Pn.lb Pn.ub]
        set_worker_problem!(prob,Pn,CertSettings(verbose=0,storage_level=0))
        distributed_region_certify([R0],Inf)
        # Do not keep any workspace (C memory) in the precompiled module
        DAQP.free_c_workspace(WORKER_WS[].DAQP_workspace)
        WORKER_WS[] = nothing
        WORKER_PROBLEM[] = nothing
    end
end
end
