module ASCertain

using DAQP, LinearAlgebra

export certify,
       Region,
       DualCertProblem,DualLPCertProblem,
       MPQP,QP,
       CertSettings
@enum State REMOVE ADD OPTIMAL INFEASIBLE ITERLIM UNBOUNDED
include("types.jl");
include("ascert.jl");


include("dualQPCert.jl");
include("dualLPCert.jl");

include("flops.jl")
include("utils.jl");

# Implementation of solvers
include("solvers/dsimplex.jl");
end
