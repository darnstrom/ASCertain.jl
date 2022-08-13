
## Update set of active sets 
function update_ASs(ASs::BitMatrix, AS::BitVector)
    (m,n) = size(ASs);
    if(n==0)
        return AS[:,:]
    end
    AS_found = 0;
    for j in 1:n
        AS_found = 1;
        for i in 1:m
            if (AS[i] != ASs[i,j])
                AS_found=0;
                break
            end
        end
        if(AS_found==1)
            break # AS already exists in ASs 
            #TODO return the index also for later use... 
        end
    end
    if(AS_found==0)
        return [ASs AS]
    else
        return ASs 
    end
end
function get_unique_ASs(part::Vector{Region})
    N = length(part)
    if(N==0)
        println("The partition is empty");
        return
    end
    ASs_unique = part[1].ASs[:,end:end];
    for i = 2:N
        ASs_unique = update_ASs(ASs_unique,part[i].ASs[:,end])
    end
    return ASs_unique
end

## Print ASs in a readable way
function print_ASs(ASs::BitMatrix)
    AS = Int64[];
    inds = collect(1:size(ASs,1))
    for i in 1:size(ASs,2)-1
        println("$AS ")
        if(sum(ASs[:,i])<sum(ASs[:,i+1])) # Addition	
            add_ind = inds[.!(ASs[:,i].⊻ .!ASs[:,i+1])]
            push!(AS,add_ind[1])
            printstyled("+$(lpad(add_ind[1],2," ")) "; color = :green)
        else # Removal
            rm_ind = inds[.!(ASs[:,i].⊻ .!ASs[:,i+1])]
            deleteat!(AS,findfirst(AS.==rm_ind[1]))
            printstyled("-$(lpad(rm_ind[1],2," ")) "; color = :red)
        end
    end
    println("$AS")
end
## Check containment in partition 
function pointlocation(th::Vector{Float64}, partition::Vector{Region};eps_gap=0.0)
    contained_in= Int64[]
    for (ind,region) in enumerate(partition)
        violation = minimum(region.bth-region.Ath'*th)
        if(violation>=-eps_gap)
            push!(contained_in,ind)
        end
    end
    return contained_in
end

## Parametric forward/backward substitution 
function forward_L_para(L,b)
    # Solve L x = b
    n = size(b,2);
    x = deepcopy(b);
    l = 0.0;
    for i in 1:n
        for j in 1:(i-1)
            l = L[i,j]
            x[:,i] -= l*x[:,j];
        end
    end
    return x
end

# Row instead of column vector
function backward_L_para!(L,x)
    # Solve L'x = b
    n = size(x,2);
    l=0.0;
    for i = n:-1:1
        for j = i+1:n
            l = L[j,i];
            x[:,i] -= l*x[:,j];
        end
    end
end

## Generate random mpQP
function generate_mpQP(n,m,nth)
    M = randn(n,n)
    H = M*M'
    f = randn(n,1)
    f_theta = randn(n,nth)
    A = randn(m,n)
    b = [rand(m,1);rand(m,1)]
    F0 = randn(n,nth); # The point F0*th will be primal feasible
    W =[A;-A]*(-F0);
    bounds_table = [collect(m+1:2m);collect(1:m)]
    senses = zeros(Cint,2m)
    mpQP = MPQP(H,f,f_theta,zeros(0,0),
                [A;-A],b,W,bounds_table,senses)

    P_theta = (A = zeros(nth,0), b=zeros(0), ub=ones(nth),lb=-ones(nth),F0=F0) 

    return mpQP,P_theta
end
## Remove redundant  
# Ar,br = remove_redundant(A,b)
# remove constraints for the polyhedron P = {x : A' x ≤ b} 
# such that {x : Ar' x ≤ br}  = {x : A' x ≤ b} 
function remove_redundant(A,b;sense=[],max_radius=1e30)

    # Setup DAQP workspace 
    nth,m = size(A)  
    ms = length(b)-m;
    p=DAQP.setup_c_workspace(nth);
    blower = fill(-1e30,m);
    if(isempty(sense))
        sense = zeros(Cint,m)
    end
    DAQP.init_c_workspace_ldp(p,A,b,blower,sense;max_radius)
    unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,3)),m) # set m 
    unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,4)),ms) # set ms 

    # Start finding redundant constraints
    is_redundant = -ones(Cint,m); 
    for i = 1:m
        (is_redundant[i]!= -1) && continue; # Decided from previous iteration 

        ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);

        sense[i] = 5; # Force ith constraint to equality
        test =ccall((:add_constraint,DAQP.libdaqp),Cint,(Ptr{Cvoid},Cint,Float64),p,i-1,1.0)

        # Check if system is infeasible (infeasible => reudandant) 
        exitflag =ccall((:daqp_ldp,DAQP.libdaqp), Int32, (Ptr{Cvoid},),p);
        ws = unsafe_load(Ptr{DAQP.Workspace}(p))
        AS = unsafe_wrap(Vector{Cint}, ws.WS, ws.n_active, own=false)
        if(exitflag==-1)
            is_redundant[i]=1
        else
            is_redundant[i]=0
            sense[i] &=~4; # Should be modifiable in later iterations 
            if(exitflag==1) # All activate constraints must also be nonredundant 
                is_redundant[AS.+1] .=0; 
            end
        end
        sense[AS.+1].&=~1; # Deactivate TODO: make sure equality constraint are not deactivated 
    end
    # Free DAQP workspace
    DAQP.free_c_workspace(p);
    nonred_ids = findall(is_redundant.==0)
    return A[:,nonred_ids],b[nonred_ids] 
end
