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

The distributed path first expands the search tree breadth first on the main process until there are about `opts.distributed_region_factor` pending regions per worker. The pending regions are then queued as jobs, with the regions that have used the fewest iterations first, and each idle worker takes the next job. A job that runs longer than `opts.distributed_time_budget` seconds returns its remaining regions, which are queued as new jobs, so that large subtrees are split over several workers. Each worker uses its own `CertWorkspace`, so the final partition is the same as for the sequential algorithm; only the execution is parallelized.

Current limitations of the distributed mode:

1. It requires the default overflow handling and output limit.
2. It does not support callbacks (`rm_callbacks`, `add_callbacks`, `termination_callbacks`, `pop_callbacks`, `conditioned_callbacks`).
3. Workers should be started with the same project environment, which is typically what you want on an HPC cluster as well.

## Settings for efficiency
The following settings in `CertSettings` reduce the time spent on the feasibility problems (LPs) that are solved during the certification. They are all enabled by default and do not change the certified partition.

* `use_witnesses`: Known points in a region (e.g., the solutions of earlier feasibility problems) are used to decide feasibility without solving an LP.
* `use_region_box`: Each region is bounded by a box (propagated from its parent), which is used to decide infeasibility without solving an LP.
* `warm_start_ldp`: The feasibility problems of a region are warm started from the optimal working set of the region.

In addition, `box_tightening_threshold` (default `0`, i.e., disabled) tightens the box of a region with 2n_θ feasibility problems if at least this many candidate constraints are not decided by the box. This can pay off for problems with many constraints.
