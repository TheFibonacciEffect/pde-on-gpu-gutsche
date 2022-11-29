# Porous Convection

## Introduction
Modelling of multi-physical processes poses challanges which can be addressed using HPC. This couse introduces the psudo-transient method to targets clusters featuring nodes with multiple GPUs which reward fine-grained parallisation of homogeneous tasks.

The aim of this procet is to to extend a previously implemented pseudo-transient porous convection solver in two dimensions written for the use on a single CPU architecture, to three dimensions and shared memory applications. This is achieved through several intermediate steps; 1) enhancing the code for the use on CPU/GPU alike making use of the ParallelStencil.jl and CUDA.jl package, 2) extending the XPU version of the code to 3D, 3) enabling the code to be used on multi-XPU architectures using MPI.jl and ImplizitGlobalGrid.jl. 

## Physics
The underlying equations are:


## Numerical methods

## Results

[![Build Status](https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche/actions/workflows/CI.yml/badge.svg)](https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche/actions/workflows/CI.yml)

<!-- [![Build Status](https://github.com/omlins/ParallelStencil.jl/workflows/CI/badge.svg)](https://github.com/omlins/ParallelStencil.jl/actions) -->

### Porous convection 2D
This is the temperature distribution and flux when running the porous convection 2D code on a 511 1023 with 4000 timesteps
![Fig1](docs/PorousConvection2D.gif)


### Porous convection 3D
this is the same phenomenon in 3D. Here is the temperature distribution after 2000 timesteps with `nx,ny,nz = 255,127,127`.

## Discussion/Conclusion
