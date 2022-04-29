using Serialization 

## I/O 
# Load from MAT
function loadMPQP(s::String)
  vars = deserialize(s);
  m=length(vars["mpQP"]["b"][:,1]);
  mpQP = MPQP(vars["mpQP"]["H"],
			  vars["mpQP"]["f"],
			  vars["mpQP"]["f_theta"],
			 zeros(0,0),
			  vars["mpQP"]["A"],
			  vars["mpQP"]["b"],
			  vars["mpQP"]["W"],
			  collect(1:m),
			  zeros(Cint,m)
			  );

  Ath = vars["P_theta"]["A"]'[:,:];
  bth = vars["P_theta"]["b"]
  if(isempty(Ath)) # Make sure the nth dimension is correct.
	Ath=zeros(size(mpQP.W,2),0);
	bth = zeros(0);
  end
  P_theta =	(A = Ath,
			 b = bth,
			 lb = vars["P_theta"]["lb"][:,1],
			 ub = vars["P_theta"]["ub"][:,1]);
  return mpQP,P_theta
end

function exportMPQP(mpQP::MPQP, P_theta, s::String)
  mpQP_d = Dict("H"=>mpQP.H, 
		   "f"=>mpQP.f, 
		   "f_theta"=>mpQP.f_theta, 
		   "H_theta"=>mpQP.H_theta,
		   "A"=>mpQP.A,
		   "b"=>mpQP.b,
		   "W"=>mpQP.W)
  P_theta_d = Dict("A"=>P_theta.A,
			   "b"=>P_theta.b,
			   "lb"=>P_theta.lb,
			   "ub"=>P_theta.ub)

  d = Dict("mpQP"=>mpQP_d, "P_theta"=>P_theta_d)
  open(s, "w") do f
	serialize(f, d)
  end
  return
end

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
