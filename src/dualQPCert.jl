## Parametric AS iteration 
function parametric_AS_iteration(prob::DualCertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,S::Vector{Region})
  # Check if iteration limit is reached
  if(region.iter >= opts.iter_limit)
	region.state=ITERLIM
	if(any(isnan.(region.Ath))||any(isnan.(ws.Ath[:,1:region.start_ind]))) 
	  println("Error in ITERLIM!")
	end
	terminate(region,ws,opts,opts.storage_level) 
	return
  end

  # Perform usual iteration
  if(region.state==ADD)
	cert_add_constraint(prob,region,opts,ws,S);
  else
	cert_remove_constraint(prob,region,opts,ws,S);
  end
end

## Add constraint
function cert_add_constraint(prob::DualCertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,partition::Vector{Region})
  
  # No indices can be added -> Global optimum
  cands = region.IS[:];
  cands[prob.bounds_table[region.AS]].=false;# Both bounds cannot be active 
  cands[region.feasible_constrs].=false;# These constraints cannot be primal infeasible 
  ind_cands = findall(cands); 
  #ind_cands = setdiff(ind_cands,region.feasible_constrs);

  if(isempty(ind_cands))
	region.state=OPTIMAL;
	if(any(isnan.(region.Ath))||any(isnan.(ws.Ath[:,1:region.start_ind]))) 
	  println("Error in Optimal (early)!")
	end
	terminate(region,ws,opts,opts.storage_level);
	return nothing 
  end
  
  # Update feasibility model 
  update_feas_model(region,ws)
  
  # Compute μ
  μ = compute_slack(region,prob,ind_cands)

  # Update flop count
  flops_compute_slack(region,size(prob.M,2))
  
  # Execute callbacks
  for callback in opts.add_callbacks
	callback(region,ws,opts,μ);
  end
  # Set primal tolerance
  ϵp = opts.eps_primal;
  # Only i s.t. μᵢ(θ) < -ϵ can be added
  negative_cands = collect(1:length(ind_cands));
  pos_cands = Int64[]; 
  prune_candidates(μ,ws,ϵp-opts.delta_mu,opts.eps_gap,negative_cands,pos_cands);
  region.feasible_constrs[ind_cands[pos_cands]].=1;

  # Θʲ=∅ ∀j ⟹ Θ* = Θ  
  if(isempty(negative_cands))
    region.state=OPTIMAL;
	region.Ath= zeros(prob.n_theta,0); 
	region.bth= zeros(0); 
    region.ASs=[region.ASs falses(size(region.ASs,1))];
    region.ASs[region.AS,end].=true;
	if(any(isnan.(region.Ath))||any(isnan.(ws.Ath[:,1:region.start_ind]))) 
	  println("Error in Early!")
	end
	terminate(region,ws,opts,opts.storage_level);
    return nothing;
  end

  # Compute added regions
  n_cands = length(negative_cands);
  norm_factor=0.0;
  Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+n_cands)];
  bth_tmp = @view ws.bth[(ws.m+1):(ws.m+n_cands)];
  for (k_outer,i) in enumerate(negative_cands)
	k=0;
	for j in negative_cands
	  if(i==j)
		# μᵢ(θ) < -ϵ 
		Ath_tmp[:,k+1] .= μ[1:end-1,i];
		bth_tmp[k+1] = -μ[end,i]-ϵp
		bth_tmp[k+1] += opts.delta_mu; # Take into account round-off errors 
		# Normalize
	  else
		# μᵢ(θ) < μⱼ(θ)
		Ath_tmp[:,k+1] = μ[1:end-1,i]-μ[1:end-1,j];
		bth_tmp[k+1] = -μ[end,i]+μ[end,j]-opts.eps_zero;
		bth_tmp[k+1] += 2*opts.delta_mu  # Take into account round-off errors
	  end
	  norm_factor = norm(Ath_tmp[:,k+1],2);
	  if(norm_factor < opts.eps_zero)
		if(bth_tmp[k+1]<0) 
		  k=-1; # mark infeasible with k=-1
		  break 
		end
	  else
		k+=1;
		Ath_tmp[:,k]./=norm_factor
		bth_tmp[k]/=norm_factor
		bth_tmp[k]-=opts.eps_gap
	  end
	end
	# Check if Θᵢ ≂̸ ∅
	if(k>=0 && isfeasible(ws,k))
	  new_region=spawn_region(region,ind_cands[i],Ath_tmp[:,:],bth_tmp[:],Float64[],prob);
	  push!(partition,new_region);
	end
  end
  
  # Check if there are any parmameters in Θ which leads to μ(θ) ≥ - ϵ
  norm_factor = 0.0;
  k=0;
  for neg_ind in negative_cands #for j ∈ AS  : j≂̸i
	norm_factor = norm(μ[1:end-1,neg_ind],2);
	if(norm_factor<1e-14)
	  (μ[end,neg_ind]< -(ϵp+opts.delta_mu))&& return nothing # trivally infeasible  
	else
	  k+=1
	  ws.Ath[:,ws.m+k] = -μ[1:end-1,neg_ind];
	  ws.bth[ws.m+k] = μ[end,neg_ind]+ ϵp
	  ws.bth[ws.m+k] += opts.delta_mu; # Take into account round-off errors
	  # Normalize
	  ws.Ath[:,ws.m+k]./=norm_factor;
	  ws.bth[ws.m+k]/=norm_factor;
	  ws.bth[ws.m+k]-=opts.eps_gap
	end
  end
  if(isfeasible(ws,k))
	ws.m+=k #Update global model (since last)  
	region.state=OPTIMAL;
	region.start_ind =ws.m;
	region.ASs=[region.ASs falses(size(region.ASs,1))];
	region.ASs[region.AS,end].=true;
	if(any(isnan.(region.Ath))||any(isnan.(ws.Ath[:,1:region.start_ind]))) 
	  println("Error in Optimal (true)!")
	end
	terminate(region,ws,opts,opts.storage_level);
  end

  return nothing 
end

## Remove constraint
function cert_remove_constraint(prob::DualCertProblem,region::Region,opts::CertSettings,ws::CertWorkspace,partition::Vector{Region})
  
  # Empty AS => Trivial CSP 
  if(isempty(region.AS)) 
    region.state=ADD;
    region.Lam= zeros(prob.n_theta+1,0);
    push!(partition,region);
    return nothing 
  end

  # Append constraints 
  update_feas_model(region,ws)
  
  n_active = length(region.AS);
  
  # TODO Singulariy only needs to be checked directly after addition...
  # TODO Pivoting...
  # Check singularity
  singular_ind = findfirst(region.D.<1e-14)
  singular = !isnothing(singular_ind);

  # Compute λ* 
  if(!singular)
	LamStar = -forward_L_para(region.L,prob.d[:,region.AS]);
	for j in 1:length(region.D) 
	  LamStar[:,j] /=region.D[j];
	end
	backward_L_para!(region.L,LamStar);
	#Update flops count
	flops_solve_kkt(region)
  else
	#Update flops count
	flops_singular_direction(region)
  end

  # Execute callbacks
  for callback in opts.rm_callbacks
	callback(region,ws,opts,LamStar);
  end

  # === Partition Θ depending on constraints that are removed  === 
  if(region.add_ind == 0) # No constraint has been added to AS yet  (Case b2)
	# TODO  can we use region.feasible_constrs?
	#feas_inds=collect(1:n_active);
	feas_inds = Int64[]; # Only inequality constraints can be removed
	for (k,ind) in enumerate(region.AS)
	  if(prob.senses[ind]&DAQP.IMMUTABLE==0)
		push!(feas_inds,k);
	  end
	end
	# Only keep i such that λ*ᵢ(θ) < -ϵ for some θ∈Θ 
	prune_candidates(LamStar[:,feas_inds],ws,opts.eps_dual,opts.eps_gap,feas_inds,Int64[]);
	# Candidates... 
	DeltaLam = -LamStar+region.Lam;
	Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+length(feas_inds))];
	bth_tmp = @view ws.bth[(ws.m+1):(ws.m+length(feas_inds))];

	for (k_outer,i) in enumerate(feas_inds)
	  k=0
	  for j in filter(x->x!=i,feas_inds) #for j ∈ AS  : j≂̸i
		# α_i(θ) < α_j(θ)
		# TODO ===== Make sure that H*_k H cancels in the dual case.. =======  
		Ath_tmp[:,k+1] .= -region.Lam[end,j].*DeltaLam[1:end-1,i].+region.Lam[end,i].*DeltaLam[1:end-1,j];
		bth_tmp[k+1] = region.Lam[end,j]*DeltaLam[end,i]-region.Lam[end,i]*DeltaLam[end,j]-opts.eps_zero;
		# Normalize
		norm_factor = norm(Ath_tmp[:,k+1],2);
		if(norm_factor < 1e-14)
		  if(bth_tmp[k+1]<0) 
			k=-1; # mark infeasible with k=-1
			break 
		  end
		else
		  k+=1;
		  Ath_tmp[:,k]./=norm_factor;
		  bth_tmp[k]/=norm_factor;
		  bth_tmp[k]-=opts.eps_gap;
		end
	  end
	  (k==-1) && continue # trivially infeasible 

	  # λ*ᵢ(θ) < -ϵ
	  Ath_tmp[:,k+1] .= LamStar[1:end-1,i];
	  bth_tmp[k+1] = -LamStar[end,i]-opts.eps_dual;
	  # Normalize
	  norm_factor = norm(Ath_tmp[:,k+1],2);
	  if(norm_factor<1e-14)
		(bth_tmp[k+1]<0) && continue  #trivially infeasible
	  else	
		k+=1
		Ath_tmp[:,k]./=norm_factor;
		bth_tmp[k]/=norm_factor;
		bth_tmp[k]-=opts.eps_gap;
	  end

	  # Check if Θᵢ≂̸ ∅
	  if(isfeasible(ws,k))
		new_region = spawn_region(region,-i,Ath_tmp[:,:],bth_tmp[:],Float64[],prob);
	    push!(partition,new_region);
	  end
	end

  elseif(region.add_ind > 0) # A constraint had been added to AS (Case b1) 
	# Compute λ* and/or p̂ 
	if(!singular)
	  p̂= -DAQP.compute_singulardirection(region.L,region.D,length(region.D))
	else
	  p̂= -DAQP.compute_singulardirection(region.L,region.D,singular_ind)
	end

	# Update feasible candidates
	region.feasible_constrs[region.feasible_constrs] .= prob.MM[region.feasible_constrs,region.AS]*p̂.<0;
	feas_inds = Int64[]
	for i in 1:(n_active-1) # Never check last because p̂[end] < 0 after addition...
	  if((prob.senses[region.AS[i]]&DAQP.IMMUTABLE == 0) && p̂[i]>opts.eps_zero)
		push!(feas_inds,i);
	  end
	end

	if(singular)
	  if(isempty(feas_inds)) # Primal infeasible problem 
		region.state=INFEASIBLE;
	    region.Ath= zeros(prob.n_theta,0); 
	    region.bth= zeros(0); 
		if(any(isnan.(region.Ath))||any(isnan.(ws.Ath[:,1:region.start_ind]))) 
		  println("Error in infeasible")
		end
		terminate(region,ws,opts,opts.storage_level);
		println("Infeasible region found")
		return nothing 
	  else
		Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+length(feas_inds)-1)];
		bth_tmp = @view ws.bth[(ws.m+1):(ws.m+length(feas_inds)-1)];
	  end
	else # Nonsingular reduced Hessian
	  # Only keep i such that λ*ᵢ(θ) < -ϵ for some θ∈Θ 
	  # TODO: eps_gap becomes normalized here...
	  prune_candidates(LamStar[:,feas_inds],ws,opts.eps_dual,opts.eps_gap,feas_inds,Int64[]);
	  Ath_tmp = @view ws.Ath[:,(ws.m+1):(ws.m+length(feas_inds))];
	  bth_tmp = @view ws.bth[(ws.m+1):(ws.m+length(feas_inds))];
	  if(isempty(feas_inds))
	    # Θʲ=∅ ∀j ⟹ Θ^CSP = Θ  
	    region.Lam = LamStar;
	    region.Ath= zeros(prob.n_theta,0); 
	    region.bth= zeros(0); 
	    region.state=ADD;
	    region.start_ind = ws.m;
	    push!(partition,region);
	    return nothing 
	  end
	end
	
	# Partition Θ based on removals from AS
	norm_factor = 0.0;
	for (k_outer,i) in enumerate(feas_inds)
	  k=0
	  for j in (filter(x->x!=i,feas_inds)) #for j ∈ AS  : j≂̸i
		# α_i(θ) < α_j(θ)
		Ath_tmp[:,k+1] .= -p̂[i].*region.Lam[1:end-1,j].+p̂[j].*region.Lam[1:end-1,i];
		bth_tmp[k+1] = p̂[i]*region.Lam[end,j]-p̂[j]*region.Lam[end,i]-opts.eps_zero;
		# Normalize
		norm_factor = norm(Ath_tmp[:,k+1],2);
		if(norm_factor<1e-14)
		  if(bth_tmp[k+1]<0) 
			k=-1; # mark infeasible with k=-1
			break 
		  end
		else
		  k+=1
		  Ath_tmp[:,k]./=norm_factor;
		  bth_tmp[k]/=norm_factor;
		  bth_tmp[k]-=opts.eps_gap;
		end
	  end
	  (k==-1) && continue #trivially infeasible

	  if(!singular)
		# λ*ᵢ(θ) < -ϵ
		Ath_tmp[:,k+1] .= LamStar[1:end-1,i];
		bth_tmp[k+1] = -LamStar[end,i]-opts.eps_dual;
		# Normalize
		norm_factor = norm(Ath_tmp[:,k+1],2);
		if(norm_factor<1e-14)
		  (bth_tmp[k+1]<0) && continue  #trivially infeasible
		else
		  k+=1
		  Ath_tmp[:,k]./=norm_factor;
		  bth_tmp[k]/=norm_factor;
		  bth_tmp[k]-=opts.eps_gap;
		end
	  end

	  # Check if Θᵢ≂̸ ∅
	  if(isfeasible(ws,k))
		new_region = spawn_region(region,-i,Ath_tmp[:,:],bth_tmp[:],p̂,prob);
	    push!(partition,new_region);
	  end
	end
  end
  # Create Θ_CSP, i.e., θ∈Θ s.t. λ*(θ)≥-ϵ
  if(singular)  
	return nothing# CSP not defined for singular case
  end
  norm_factor = 0.0;
  k=0
  for feas_ind in feas_inds #for j ∈ AS  : j≂̸i
	ws.Ath[:,ws.m+k+1] = -LamStar[1:end-1,feas_ind];
	ws.bth[ws.m+k+1] = LamStar[end,feas_ind]+opts.eps_dual;
	# Normalize
	norm_factor = norm(ws.Ath[:,ws.m+k+1],2);
	if(norm_factor<1e-14)
	  (ws.bth[ws.m+k+1]<0) && return nothing #trivially infeasible
	else
	  k+=1
	  ws.Ath[:,ws.m+k]./=norm_factor;
	  ws.bth[ws.m+k]/=norm_factor;
	  ws.bth[ws.m+k]-=opts.eps_gap;
	end	
  end
  if(isfeasible(ws,k))
	ws.m+=k; #update model
	region.Lam = LamStar;
	region.Ath= zeros(prob.n_theta,0); 
	region.bth= zeros(0); 
	region.state=ADD;
	region.start_ind = ws.m;
	push!(partition,region);
  end
  return nothing 
end
## Spawn region
function spawn_region(region::Region, i::Int64, Ath::Matrix{Float64}, bth::Vector{Float64}, p̂, prob::DualCertProblem)
  n_active = length(region.AS);
  new_region=Region(region.IS[:], region.AS[:],
					Ath,bth, REMOVE,region.iter+1,region.start_ind,region.add_ind,region.reuse_ind,
					Array{Float64}(undef,size(region.Lam,1),n_active+sign(i)),
					BitArray(undef,size(region.ASs).+(0,1)),
					Array{Float64}(undef,0,0),
					Array{Float64}(undef,0),
					region.feasible_constrs[:],
					deepcopy(region.kappa));
  if(i>0) # add
	m = prob.M[i,:];
	new_region.L,new_region.D=DAQP.updateLDLadd(region.L,region.D,prob.M[region.AS,:]*m,m'*m);
	flops_add_constraint(region,size(prob.M,2));
	push!(new_region.AS,i);
	new_region.IS[i] = false;
	new_region.Lam[:,1:end-1].=region.Lam;
	new_region.Lam[:,end].=0;
	new_region.add_ind=i;
	new_region.reuse_ind=length(new_region.AS)-1;
	

  else # remove 
	i = abs(i)
	new_region.L,new_region.D=DAQP.updateLDLremove(region.L,region.D,i);
	flops_remove_constraint(i,region,length(bth))
	new_region.IS[region.AS[i]]=true;
	deleteat!(new_region.AS,i);
	if(!isempty(p̂))
	  for (k,j) in enumerate(filter(x->x!=i,1:n_active)) #Remove constraint from AS
		new_region.Lam[:,k].= region.Lam[:,j].-(p̂[j]/p̂[i]).*region.Lam[:,i];
	  end
	else
	  new_region.Lam.= region.Lam[:,1:end.!=i];
	end
	new_region.reuse_ind=i-1;

  end
  # Update parent
  new_region.ASs[:, 1:end-1] = region.ASs;
  new_region.ASs[:,end]=.~region.IS;
  return new_region
end
## Compute slack 
function compute_slack(region,prob::DualCertProblem,ind_cands)
  μ = prob.d[:,ind_cands];
  if(~isempty(region.AS))
	μ+= region.Lam*prob.MM[region.AS,ind_cands]
  end
  return μ
end

## Normalize problem to a box -1 ≤ θ ≤ 1
function normalize(prob::DualCertProblem,P_theta,mpQP)
  # TODO also normalize TH part of P_theta
  prob = deepcopy(prob);
  nth = length(P_theta.lb);
  center = (P_theta.lb+P_theta.ub)/2;
  norm_factors = (P_theta.ub-P_theta.lb)/2;
  for i in 1:nth
	prob.d[i,:] *= norm_factors[i];
	mpQP.f_theta[:,i]*=norm_factors[i];
	mpQP.W[:,i]*=norm_factors[i];
  end
  prob.d[end:end,:] -= center'*prob.d[1:end-1,:];
  mpQP.b -= mpQP.W*(center)
  mpQP.f -= mpQP.f_theta*(center)
  P_theta =(A=P_theta.A,
			b = P_theta.b,
			lb = -ones(nth),
			ub = ones(nth));
  return prob, P_theta, mpQP
end
## Slice mpQP
function reduce_problem(mpQP, P_theta,ids;slice_vals=[])
  nth = size(mpQP.W,2)
  slice_ids= setdiff(1:nth,ids) #ids to remove 
  if(isempty(slice_vals)) # Default slice at 0
	slice_vals = zeros(length(slice_ids))
  end
  
  # Reduce the mpQP
  mpQP.f += mpQP.f_theta[:,slice_ids]*slice_vals
  mpQP.b += mpQP.W[:,slice_ids]*slice_vals
  mpQP.f_theta = mpQP.f_theta[:,ids] 
  mpQP.W= mpQP.W[:,ids] 
  # TODO also acount for H_theta?
  
  # Reduce the region 
  b_new = isempty(P_theta.b) ? P_theta.b :  P_theta.b + slice_vals'*P_theta.A[slice_ids,:]
  A_new = P_theta.A[ids,:];

  P_theta_new = (A=A_new,
				 b=b_new,
				 lb= P_theta.lb[ids], ub = P_theta.ub[ids])  
  return mpQP,P_theta_new
end
