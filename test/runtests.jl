using ASCertain
using DAQP
using Test
using LinearAlgebra
#include(joinpath(dirname(@__FILE__), "utils.jl"))


# Generate example mpQP
n,m,nth = 5,5,4
mpQP,P_theta = ASCertain.generate_mpQP(n,m,nth) 
opts = CertSettings();
prob = DualCertProblem(mpQP); 
prob,P_theta,mpQP = ASCertain.normalize(prob,P_theta,mpQP);
opts.verbose=0
opts.storage_level=2

@testset "Cold start" begin
  # Run certificatoin
  AS0 = Int64[];
  (part,iter_max) = certify(prob,P_theta,AS0,opts);
  println("Part: $(length(part))")
  # Test for random samples 
  N = 1000
  ths = 2*rand(nth,N).-1;

  containment_inds = Int64[];
  diff_iters = zeros(Int64,N);
  for n = 1:N
	th = ths[:,n]
	inds = ASCertain.pointlocation(th,part,eps_gap=opts.eps_gap);
	push!(containment_inds,length(inds))
	if(length(inds)>0)
	  x,lam,AS,J,iter = DAQP.daqp_jl(QP(mpQP,th),deepcopy(AS0));
	  diff_iters[n] = part[inds[1]].iter-iter;
	end
  end
  println(sum(containment_inds.==0))
  @test ~any(containment_inds.==0) # No holes
  @test ~any(containment_inds.>1) # No overlap
  @test sum(abs.(diff_iters))==0 # Agreement with MC simulations
end

@testset "Warm start" begin
  # Run certificatoin
  AS0 = Int64[1,2];
  (part,iter_max) = certify(prob,P_theta,AS0,opts);
  # Test for random samples 
  N = 1000
  ths = 2*rand(nth,N).-1;

  containment_inds = Int64[];
  diff_iters = zeros(Int64,N);
  for n = 1:N
	th = ths[:,n]
	inds = ASCertain.pointlocation(th,part,eps_gap=opts.eps_gap);
	push!(containment_inds,length(inds))
	if(length(inds)>0)
	  x,lam,AS,J,iter = DAQP.daqp_jl(QP(mpQP,th),deepcopy(AS0));
	  diff_iters[n] = part[inds[1]].iter-iter;
	end
  end
  @test ~any(containment_inds.==0) # No holes
  @test ~any(containment_inds.>1) # No overlap
  @test sum(abs.(diff_iters))==0 # Agreement with MC simulations
end
