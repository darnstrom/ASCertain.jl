# ASCertain
[![CI](https://github.com/darnstrom/ASCertain.jl/workflows/CI/badge.svg)](https://github.com/darnstrom/ASCertain.jl/actions)
[![Code coverage](http://codecov.io/gh/darnstrom/ASCertain.jl/graphs/badge.svg)](http://codecov.io/github/darnstrom/ASCertain.jl)

An implementation of the parametric complexity certification method presented in the article [A unifying complexity certification framework for active-set methods for convex quadratic programming](https://ieeexplore.ieee.org/abstract/document/9461599). The package is specifically adapted for certifying the complexity of the dual active-set QP solver [daqp](https://github.com/darnstrom/daqp).

## Basic example  
The following code applies the certification method on the mpQP labeled "contrived mpQP" in the paper [A unifying complexity certification framework for active-set methods for convex quadratic programming](https://ieeexplore.ieee.org/abstract/document/9461599). 
```julia
using ASCertain 
# Setup problem data
H = [0.97 0.19 0.15;0.19 0.98 0.05; 0.15 0.05 0.99];
f = zeros(3,1);
f_theta = [11.3 -44.3; -3.66 -11.9; -32.6 7.81];
H_theta = zeros(0,0);
A = [2.5*0.15 2.5*0.88 2.5*0.17;0.49 0.57 0.22; 0.77 0.46 0.41];
b = zeros(3,1);
b[1:3] = [4.1;3.7;4.3];
W = [0.19 -0.89; 0.64 -1.54; -0.59 -1.01];

# Create mpQP
mpQP=ASCertain.MPQP(H,f,f_theta,H_theta,A,b,W)

# Create region of interest
P_theta = (A=zeros(2,0), b=zeros(0),lb=zeros(2),ub=1.5*ones(2))

# Run certification with default settings and an empty working set 
opts = CertSettings();
opts.storage_level = 2; # Store all regions
AS = Int64[]
@time (part,iter_max) = certify(mpQP,P_theta,AS;opts);
```

## Distributed certification

`certify` can also distribute independent certification subtrees over Julia workers:

```julia
using Distributed
addprocs(4; exeflags="--project=$(Base.active_project())")
@everywhere using ASCertain

opts = CertSettings(verbose=0)
(part, iter_max) = certify(mpQP, P_theta, AS; opts, workers=workers())
```

The distributed path first expands the search tree on the main process until there are enough pending regions, then sends those regions to workers. During this seeding stage, ASCertain expands the pending region with the fewest certification iterations first, which keeps the distributed subtrees more even without changing the certified partition. Each worker uses its own `CertWorkspace`, so the final partition is the same as for the sequential algorithm; only the execution is parallelized. The seeding granularity is controlled by `opts.distributed_region_factor`.

Current limitations of the distributed mode:

1. It requires the default overflow handling and output limit.
2. It does not support callbacks (`rm_callbacks`, `add_callbacks`, `termination_callbacks`, `pop_callbacks`, `conditioned_callbacks`).
3. Workers should be started with the same project environment, which is typically what you want on an HPC cluster as well.

Other parallelization options:

1. **Initial geometric partition of** `P_theta`: easy to implement and maps well to cluster job arrays, but it changes the region boundaries unless the resulting partition is merged afterward, and that merge can be expensive.
2. **Dynamic work stealing over certification regions**: can improve load balancing when subtree sizes differ a lot, but it needs a distributed queue, more synchronization, and extra care to keep the returned ordering deterministic.
3. **Thread-level parallelism inside a single region step**: could speed up candidate pruning or feasibility checks without shipping regions between workers, but the current implementation relies on a mutable DAQP workspace, so safe threading would require larger internal refactoring.
4. **Hybrid seeding plus distributed workers**: seed a moderate number of regions centrally and then let workers pull more work dynamically. This can give the best scalability on large trees, but it is more complex to maintain and harder to make reproducible.
