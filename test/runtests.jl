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

@testset "LPcert" begin
  n,m,nth = 5,20,4;
  f = randn(n);
  A = randn(m,n);
  b =  rand(nth+1,m);
  bounds_table=collect(1:m);
  senses = zeros(Cint,m);

  N=1000;
  ths = 2*rand(nth,N).-1;

  prob = ASCertain.DualLPCertProblem(f,A,b,nth,n,bounds_table,senses)
  P_theta = (A = zeros(nth,0), b=zeros(0), ub=ones(nth),lb=-ones(nth)) 
  part,max_iter = certify(prob,P_theta,Int64[],opts);

  containment_inds = Int64[];
  diff_iters = zeros(Int64,N);
  for n = 1:N
	th = ths[:,n]
	inds = ASCertain.pointlocation(th,part,eps_gap=opts.eps_gap);
	push!(containment_inds,length(inds))
	if(length(inds)>0)
	  x,lam,~,~,iter = ASCertain.dsimplex(f,A,b[1:end-1,:]'*th+b[end,:],Int64[]);
	  diff_iters[n] = part[inds[1]].iter-iter;
	end
  end
  @test ~any(containment_inds.==0) # No holes
  @test ~any(containment_inds.>1) # No overlap
  @test sum(abs.(diff_iters))==0 # Agreement with MC simulations

  # Unbounded
  n,m,nth = 5,1,4; # m<n => unboudned LP
  f = randn(n);
  A = randn(m,n);
  b =  rand(nth+1,m);
  bounds_table=collect(1:m);
  senses = zeros(Cint,m);

  prob = ASCertain.DualLPCertProblem(f,A,b,nth,n,bounds_table,senses)
  P_theta = (A = zeros(nth,0), b=zeros(0), ub=ones(nth),lb=-ones(nth)) 
  part,max_iter = certify(prob,P_theta,Int64[],opts);
  @test length(part)==1 && part[1].state == ASCertain.UNBOUNDED
end
