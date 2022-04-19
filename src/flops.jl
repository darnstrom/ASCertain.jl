function flops_add_constraint(region,n) 
  # Update LDL after addition
  nA = length(region.AS)
  
  Fadd = nA*n +n + nA*(nA-1)/2 
  Fmult= nA*(n+2)+n+nA*(nA-1)/2 
  Fdiv = nA 
  Fcomp = 0
  region.kappa[:flops].+=[Fadd,Fmult,Fdiv,Fcomp];
  return
end

function flops_remove_constraint(i,region,nB) 
  nA = length(region.AS)
  j = nA-i

  # Update LDL after removal 
  Fadd = 2*nA+j*(j+2) 
  Fmult= j*(j+5)+nA 
  Fdiv = j+nB 
  Fcomp = nB-1 
  region.kappa[:flops].+=[Fadd,Fmult,Fdiv,Fcomp];
  return
end

function flops_solve_kkt(region) 
  # Compute λ
  nA = length(region.AS)
  reuse_id = region.reuse_ind 

  Fadd = nA*(nA-1)-reuse_id*(reuse_id-1)/2
  Fmult= nA*(nA-1)-reuse_id*(reuse_id-1)/2
  Fdiv = nA-reuse_id
  Fcomp = nA #(From checking λ≥0)
  region.kappa[:flops].+=[Fadd,Fmult,Fdiv,Fcomp];
  return
end

function flops_singular_direction(region) 
  # Compute p̂
  singular_id = 0 #TODO
  
  Fadd = (singular_id-1)*(singular_id-2)/2 
  Fmult= (singular_id-1)*(singular_id-2)/2 
  Fdiv = 0 
  Fcomp = 0 
  region.kappa[:flops].+=[Fadd,Fmult,Fdiv,Fcomp];
  return
end

function flops_compute_slack(region,n) 
  # Compute μ 
  # TODO: take into account upper/lower
  nA = length(region.AS)
  nI = length(region.IS)

  Fadd = n*(nA-1)+n*nI 
  Fmult= n*(nA+nI)
  Fdiv = 0 
  Fcomp = nA 
  region.kappa[:flops].+=[Fadd,Fmult,Fdiv,Fcomp];
  return
end
