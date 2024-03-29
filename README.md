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
