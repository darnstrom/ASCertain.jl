using ASCertain
using DAQP
using PolyDAQP
using Test
using LinearAlgebra

function is_finished(r::Region)
    return(r.state ∈ [ASCertain.OPTIMAL, 
                      ASCertain.INFEASIBLE, 
                      ASCertain.ITERLIM, 
                      ASCertain.UNBOUNDED]) 
end

# Generate example mpQP
n,m,nth = 5,5,4
mpQP,P_theta = ASCertain.generate_mpQP(n,m,nth) 

opts = CertSettings();
opts.verbose=0
opts.compute_flops=true
opts.storage_level=1
opts.compute_chebyball=true
opts.store_ASs=true

@testset "Cold start" begin
    # Run certificatoin
    AS0 = Int64[];
    (part,iter_max,N_fin,ASs, bin) = certify(mpQP,P_theta,AS0;opts);
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

    # Check that analysis is correct in Chebyshev centers 
    diff_iters = Int[]
    for p in part
        θ,r = p.chebyball
        θ ∉ Polyhedron(p.Ath,p.bth) && continue # Might be hard to find a point in narrow regions
        x,lam,AS,J,iter = DAQP.daqp_jl(QP(mpQP,θ),deepcopy(AS0));
        push!(diff_iters,p.iter-iter)
    end
    @test sum(abs.(diff_iters))==0

    # Test some utils
    display(part)
    display(part[end])
    @test size(ASs,2) >= size(ASCertain.get_unique_ASs(part),2)

    # Test some other options 
    opts.minrep_regions = true
    opts.prune_subsequences = true
    AS0 = Int64[];
    (part_new,iter_max,N_fin,ASs, bin) = certify(mpQP,P_theta,AS0;opts);
    @test length(part_new) < length(part)
    Am,bm = PolyDAQP.minrep(part_new[end].Ath,part_new[end].bth)
    @test length(bm) == length(part_new[end].bth)
    opts.minrep_regions = false
    opts.prune_subsequences = false

end

@testset "Warm start" begin
    # Run certificatoin
    AS0 = Int64[1,2];
    (part,iter_max) = certify(mpQP,P_theta,AS0;opts);
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
    mpLP,P_theta = ASCertain.generate_mpQP(n,m,nth;double_sided=false) 
    mpLP.H[:,:] .= 0 # Make LP
    part,max_iter = certify(mpLP,P_theta;opts,normalize=false);

    N=1000;
    ths = 2*rand(nth,N).-1;
    containment_inds = Int64[];
    diff_iters = zeros(Int64,N);
    for n = 1:N
        th = ths[:,n]
        inds = ASCertain.pointlocation(th,part,eps_gap=opts.eps_gap);
        push!(containment_inds,length(inds))
        if(length(inds)>0)
            x,lam,exitflag,AS,iter = ASCertain.dsimplex(mpLP.f[:],mpLP.A,mpLP.b+mpLP.W*th,Int64[]);
            diff_iters[n] = part[inds[1]].iter-iter;
            if(diff_iters[n] != 0) 
                println(diff_iters[n])
                readline()
            end
        end
    end
    @test ~any(containment_inds.==0) # No holes
    @test ~any(containment_inds.>1) # No overlap
    @test sum(abs.(diff_iters))==0 # Agreement with MC simulations

    # Unbounded
    n,m,nth = 5,1,4; # m < n => unboudned LP
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

@testset "LP merged" begin
    n,m,nth = 5,20,4;
    f = randn(n);
    A = randn(m,n);
    b =  rand(nth+1,m);
    bounds_table=collect(1:m);
    senses = zeros(Cint,m);

    prob = ASCertain.DualLPCertProblem(f,A,b,nth,n,bounds_table,senses)
    P_theta = (A = zeros(nth,0), b=zeros(0), ub=ones(nth),lb=-ones(nth)) 
    exp_sol = merged_certify(prob,P_theta,Int64[],opts);
end

@testset "Callbacks" begin
    # Test pop callback by always returning true (leading to skipping an iteration)
    test_pop_callback=(r,p,w,o) -> true
    cb_opts =  CertSettings()
    cb_opts.verbose=0
    cb_opts.storage_level=2
    push!(cb_opts.pop_callbacks,test_pop_callback)
    (part,iter_max,N_fin,ASs,bin) = certify(mpQP,P_theta;opts=cb_opts);
    @test length(part) == 0

    # Reset pop_callbacks
    cb_opts.pop_callbacks = Function[]
    # Test termination callback by settings regions state to ADD for all terminated regions
    test_termination_callback=(reg,ws) -> (reg.state=ASCertain.ADD)
    push!(cb_opts.termination_callbacks,test_termination_callback)
    (part,iter_max,N_fin,ASs,bin) = certify(mpQP,P_theta;opts=cb_opts);
    @test length(part) > 0
    @test all([p.state == ASCertain.ADD for p in part])
end

@testset "Default Overflow" begin
    (part,iter_max,N_fin,ASs, bin) = certify(mpQP,P_theta;opts);
    of_opts = CertSettings();
    of_opts.verbose=2
    of_opts.storage_level=2
    of_opts.output_limit = ceil(N_fin/2)
    (part_of,iter_max_of,N_fin_of,ASs_of,bin_of) = certify(mpQP,P_theta;opts=of_opts);
    @test(iter_max_of == iter_max)
    @test(N_fin_of < N_fin)
end

@testset "Try simple callback" begin
    opts = CertSettings();
    condition = (S,prob,ws,opts) -> (length(S) > 10)
    function test_callback(S,prob,ws,opts)
        println("Inside test_callback | stack length: $(length(S)) | final stack length: $(length(ws.F))")
        return true  # true  means that certify will terminate 
    end
    push!(opts.conditioned_callbacks, (condition, test_callback))
    certify(mpQP,P_theta;opts);
end

@testset "Step Overflow" begin
    n,m,nth = 5,5,4
    mpQP,P_theta = ASCertain.generate_mpQP(n,m,nth)  
    opts_baseline = CertSettings()
    opts.storage_level = 2;
    (part,iter_max,N_fin,ASs, bin) = certify(mpQP,P_theta;opts=opts_baseline); # Baseline

    function my_overflow_handle(S::Vector{Region}, prob, ws, opts)
        # Implicit representation -> Explicit representation
        for region in S 
            region.Ath= [ws.Ath[:,1:region.start_ind] region.Ath]; 
            region.bth= [ws.bth[1:region.start_ind]; region.bth];
            region.start_ind = 0;
        end
        return [S;ws.F],NaN,NaN 
    end

    of_opts = CertSettings();
    of_opts.verbose=0
    of_opts.storage_level=2
    of_opts.output_limit = 1 
    of_opts.overflow_handle = my_overflow_handle


    prob = ASCertain.setup_certproblem(mpQP;normalize=true); 
    prob,P_theta,mpQP = ASCertain.normalize(prob,deepcopy(P_theta),deepcopy(mpQP));
    ws = ASCertain.setup_workspace(P_theta,1000);

    A = [Matrix{Float64}(I,nth,nth) Matrix{Float64}(-I,nth,nth) P_theta.A];
    b = [P_theta.ub;-P_theta.lb;P_theta.b];
    R0 = Region(Int64[],A,b,prob);
    S = [R0]

    part_of_final = Region[]
    while !isempty(S)
        ASCertain.reset_workspace(ws)
        part_of,_ = certify(S,prob,ws,of_opts);
        S = Region[]
        while !isempty(part_of)
            p = pop!(part_of)
            if(is_finished(p))
                push!(part_of_final,p)
            else
                push!(S,p)
            end
        end
    end
    @test(length(part_of_final) == length(part))
end

@testset "Compoute explicit solution LP" begin
    n,m,nth = 5,5,4
    mpQP,P_theta = ASCertain.generate_mpQP(n,m,nth) 
    dprob = ASCertain.setup_certproblem(mpQP)
    opts = CertSettings();
    opts.storage_level=2
    opts.verbose=0
    AS0 = Int64[];
    (part,iter_max,N_fin,ASs, bin) = certify(mpQP,P_theta,AS0;opts);
    ASCertain.explicit_solution(part[1],dprob)

    n,m,nth = 5,20,4;
    mpLP,P_theta = ASCertain.generate_mpQP(n,m,nth;double_sided=false);
    mpLP.H[:,:] .= 0 # Make LP
    part,max_iter = certify(mpLP,P_theta;opts,normalize=false);
    dprob = ASCertain.setup_certproblem(mpLP);
    ASCertain.explicit_solution(part[1],dprob)
end

@testset "Compute explicit solution" begin
    n,m,nth = 5,5,4
    mpQP,P_theta = ASCertain.generate_mpQP(n,m,nth) 
    dprob = ASCertain.setup_certproblem(mpQP)
    opts = CertSettings();
    opts.storage_level=2
    opts.verbose=0
    AS0 = Int64[];
    (part,iter_max,N_fin,ASs, bin) = certify(mpQP,P_theta,AS0;opts);
    ASCertain.explicit_solution(part[1],dprob)

    n,m,nth = 5,20,4;
    mpLP,P_theta = ASCertain.generate_mpQP(n,m,nth;double_sided=false);
    mpLP.H[:,:] .= 0 # Make LP
    part,max_iter = certify(mpLP,P_theta;opts,normalize=false);
    dprob = ASCertain.setup_certproblem(mpLP);
    ASCertain.explicit_solution(part[1],dprob)
end
