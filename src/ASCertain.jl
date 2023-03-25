module ASCertain

using DAQP, LinearAlgebra, PolyDAQP

export certify, merged_certify,
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
include("plot.jl");

# Implementation of solvers
include("solvers/dsimplex.jl");
end
