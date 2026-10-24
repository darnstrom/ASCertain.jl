module ASCertain

using DAQPBase, Distributed, LinearAlgebra, PolyDAQP
const DAQP = DAQPBase

export certify, merged_certify, never_added_constraints, explicit_regions,
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
include("explicit.jl");

include("flops.jl")
include("utils.jl");
include("plot.jl");

# Implementation of solvers
include("solvers/dsimplex.jl");
end
