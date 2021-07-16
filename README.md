# On the stability of robust dynamical low-rank approximations for hyperbolic problems

This code framework can be used to reproduce all numerical results of the paper "On the stability of robust dynamical low-rank approximations for hyperbolic problems". The code is written in the programming language Julia (Version 1.6.0).

The code uses the following Julia packages: PyPlot, NPZ, ProgressMeter, LinearAlgebra, GSL, FastTransforms, FastGaussQuadrature. Installing these packages can be done via the Julia package manager (https://docs.julialang.org/en/v1/stdlib/Pkg/).

To run the code, go into the Julia environment and type include("runAll.jl"). Since the cfl study in Figures 3 and 7 are time consuming, one can comment out the generation of these plots in lines 8 and 11.
