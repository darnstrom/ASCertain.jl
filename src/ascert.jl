## ASCertification main
function certify(prob::CertProblem,P_theta,AS::Vector{Int64},opts::CertSettings)
  nth = length(P_theta.ub);
  # Pack initial region in the form A'θ ≤ b
  A = [Matrix{Float64}(I,nth,nth) Matrix{Float64}(-I,nth,nth) P_theta.A];
  b = [P_theta.ub;-P_theta.lb;P_theta.b];
  R0 = Region(AS,A,b,prob);
  R0.kappa[:flops]= [0,0,0,0];

  ws = setup_workspace(P_theta,opts.max_constraints); # TODO larger?
  ret = certify(R0,prob,ws,opts);
  DAQP.free_c_workspace(ws.DAQP_workspace);
  return ret
end
function certify(region::Region, prob::CertProblem, ws::CertWorkspace,opts::CertSettings) 
  ws.N_fin=0; ws.iter_max=0;j=0; 
  S =[region];
  while(!isempty(S))
	region = pop!(S);
	(opts.verbose>=2)&&print("\r>> #$(j+=1)|Stack: $(length(S))|Fin: $(ws.N_fin)|Max: $(ws.iter_max)|   ");
	parametric_AS_iteration(prob,region,opts,ws,S);
  end # Stack is empty 
  if(opts.verbose>=1)
	println("\n======= Final information: =======");
	println("||Max:$(ws.iter_max)|LPs: $(ws.lp_count)|ASs: $(size(ws.ASs,2))|Fin:$(ws.N_fin)||");
	println("==================================");
  end
  return ws.F, ws.iter_max, ws.N_fin, ws.ASs, ws.bin
  #return ws.N_fin, ws.iter_max
end
## Update optimization model
function update_feas_model(reg::Region,ws::CertWorkspace)
  m = size(reg.Ath,2);
  for i in 1:m
	ws.Ath[:,reg.start_ind+i]=reg.Ath[:,i];
	ws.bth[reg.start_ind+i]=reg.bth[i];
  end
  reg.start_ind += m; 
  ws.m = reg.start_ind;  
end
## Feasibility check
function isfeasible(ws::CertWorkspace ,m_new::Int64)::Bool
  
  ws.sense_feasibility[1:ws.m+m_new].=0; # TODO Remove after new version of DAQP_jll

  p = ws.DAQP_workspace;
  # Update A and bupper/blower 
  unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,3)),ws.m+m_new)
  unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,4)),0) # TODO: change to 2*nth

  exitflag =ccall((:daqp_ldp,DAQP.libdaqp), Int32, (Ptr{Cvoid},),p);
  ws.lp_count+=1;
  feas = (exitflag==1);

  # Reset the workspace
  ccall((:deactive_constraints,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
  ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);

  return feas
end

## Prune candidates 
# Determine candidates in M that yields a nonempty region 
# when intersected with the reg.A th<= reg.b. 
# The index of such feasible rows are retained in the vector cands
function prune_candidates(M::Matrix{Float64},ws::CertWorkspace,eps::Float64,eps_gap::Float64,cands::Vector{Int64},pos_cands::Vector{Int64})
  norm_factor=0.0;
  for i in size(M,2):-1:1
	ws.Ath[:,ws.m+1] = M[1:end-1,i];
	ws.bth[ws.m+1] = -M[end,i]-eps;
	# Normalize
	norm_factor = norm(ws.Ath[:,ws.m+1],2);
	ws.Ath[:,ws.m+1]./=norm_factor;
	ws.bth[ws.m+1]/=norm_factor;
	ws.bth[ws.m+1]-=eps_gap;
	if(~isfeasible(ws,1))
	  push!(pos_cands,cands[i]);
	  deleteat!(cands,i);
	end
  end
end

## Setup workspace
function setup_workspace(P_theta,m_max)::CertWorkspace
  nth = length(P_theta.ub);
  A = Array{Float64}(undef,nth,m_max);
  b = Array{Float64}(undef,m_max);
  blower = fill(-1e30,m_max);

  max_radius =  nth*(maximum(P_theta.ub)^2); # The region is contained in a ball with this radius.
  # Create C workspace
  p=DAQP.setup_c_workspace(nth);
  ws = CertWorkspace(A,b,blower,zeros(Cint,m_max),Region[],0,0,0,p,max_radius,falses(0,0),0,Dict());
  
  DAQP.init_c_workspace_ldp(p,ws.Ath,ws.bth,ws.bth_lower,ws.sense_feasibility;max_radius)

  return ws 
end
## Reset workspace
function reset_workspace(ws::CertWorkspace)
  ws.F = Region[];
  ws.iter_max = 0;
  ws.N_fin = 0;
  ws.m=0;
  ws.ASs=falses(0,0)
  ws.lp_count=0;
end
## Get final region
function extract_regions(region::Region, ws::CertWorkspace)
  region.Ath= [ws.Ath[:,1:region.start_ind] region.Ath]; 
  region.bth= [ws.bth[1:region.start_ind]; region.bth];
end
## Terminate a region
function terminate(region::Region,ws::CertWorkspace,opts::CertSettings,storage_level::Int8) 

  # Run termination_callbacks
  for callback in opts.termination_callbacks
	callback(region,ws);
  end

  ws.iter_max = max(ws.iter_max,region.iter);
  ws.N_fin +=1;
  (storage_level==0)&& return # Store nothing
  (storage_level>=2)&& extract_regions(region,ws) # Extract regions
   
  AS_bool = falses(size(region.ASs,1))
  AS_bool[region.AS] .=true;
  ws.ASs = update_ASs(ws.ASs,AS_bool);
  push!(ws.F,region)
end
