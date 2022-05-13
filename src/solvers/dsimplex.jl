#Dual simplex method for solving min f'x s.t. A <=b 
# AS is assumed to be a dual feasible basis 
# to perform phase 1, use AS = [];
# exitflags:
# -2 - Initial point not dual feasible
#  2 - optimal solution
#  3 - primal infeasibible
# x,lam,exitflag,AS,info = dsimplex(f,A,b,AS) 
function dsimplex(f,A,b,AS) 
  iter = 0
  exitflag=0
  x=nothing
  if(isempty(AS)) # No basic starting point => perform phase 1
	lam,AS,iter_ph1,flag_ph1 = dphase1(f,A)
	iter += iter_ph1;
	if(flag_ph1!=0) # Unbounded problem
	  return x,lam,flag_ph1,AS,iter
	end
  else
	lam = -A[AS,:]'\f; # lam should be nonnegative to be feasible 
  end

  if(any(lam.<0))
	exitflag=-2
	return NaN,lam,exitflag,AS,iter; # Starting basis not dual feasible
  end
  # Setup inactive-set
  IS = trues(length(b)) 
  IS[AS].=false
  while(exitflag ==0)
	# Check primal feasibility
	x = A[AS,:]\b[AS];
	s = b[IS]-A[IS,:]*x;
	val,add_id = findmin(s);
	if(val>-1e-6) 
	  exitflag = 2 
	else
	  iter += 1 
	  # Find constraint to remove 
	  add_id = findall(IS)[add_id] # transform to global ID
	  lam,valid_pivot=pivot_add!(AS,IS,lam,add_id,A);
	  if(!valid_pivot)
		exitflag=3 #Infeasible
	  end
	end
  end
  return x,lam,exitflag,AS,iter
end

function pivot_add!(AS,IS,lam,add_id,A)
  dlam = -(A[AS,:]')\A[add_id,:];
  block_ids = findall(dlam .< 0);
  if(isempty(block_ids)) 
	return lam,false # infeasible 
  else
	alphas = -lam[block_ids]./dlam[block_ids]; 
	α,rm_id = findmin(alphas);
	rm_id = block_ids[rm_id]; # map id in block_ids to corresponding id in AS; 
	# Add constraint
	push!(AS,add_id); 
	IS[add_id] = false; 
	# Remove constraint
	IS[AS[rm_id]]=true; 
	deleteat!(AS,rm_id); 
	# Update lambda
	lam += α*dlam;
	lam[rm_id:end-1] = lam[rm_id+1:end]; # remove deactivated lambda
	lam[end] = α; # add activated lambda
	return lam,true
  end
end

function dphase1(f,A)
  m,n = size(A);
  Aph1= [A;-diagm(sign.(f))];
  bph1= [zeros(m);ones(n)];
  AS = collect((m+1):(m+n))
  ~,lam,eflag_ph1,AS,iter= dsimplex(f,Aph1,bph1,collect((m+1):(m+n)))
  flag = all(AS.<=m) ? 0 : -3;
  return lam,AS,iter,flag
end
