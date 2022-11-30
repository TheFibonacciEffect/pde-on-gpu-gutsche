# Porous Convection

## Introduction
Modelling of multi-physical processes poses challanges which can be addressed using HPC. This couse introduces the psudo-transient method to targets clusters featuring nodes with multiple GPUs which reward fine-grained parallisation of homogeneous tasks.

The aim of this procet is to to extend a previously implemented pseudo-transient porous convection solver in two dimensions written for the use on a single CPU architecture, to three dimensions and shared memory applications. This is achieved through several intermediate steps; 1) enhancing the code for the use on CPU/GPU alike making use of the ParallelStencil.jl and CUDA.jl package, 2) extending the XPU version of the code to 3D, 3) enabling the code to be used on multi-XPU architectures using MPI.jl and ImplizitGlobalGrid.jl. 

## Physics
The underlying equations are:
$$ \begin{alignat}{1}
& \boldsymbol{q}_D =-\frac{k}{\eta}\left(\nabla p-\rho_0 \alpha \boldsymbol{g} T\right) \\
& \nabla \cdot \boldsymbol{q}_{\boldsymbol{D}}=0 \\
& \boldsymbol{q}_{\boldsymbol{T}}=-\frac{\lambda}{\rho_0 c_p} \nabla T \\
& \frac{\partial T}{\partial t}+\frac{1}{\varphi} \boldsymbol{q}_{\boldsymbol{D}} \cdot \nabla T+\nabla \cdot \boldsymbol{q}_{\boldsymbol{T}}=0
\end{alignat}$$

where $q_{D}$ is the Darcy flux, $k$ is the permeability, $\eta$ id the fluid viscosity, $p$ is the pressure, $\rho_{0}$ is the density, $\alpha$ is the termal expansion coeficient, $T$ is the temperature, $q_{T}$ is the conductive heat flux, $c_{p}$ is the specific heat capacity, $t$ it the physical time and $\varphi$ is the porosity.

## Numerical methods

## Results

[![Build Status](https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche/actions/workflows/CI.yml/badge.svg)](https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche/actions/workflows/CI.yml)

[![Literate Status](https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche/actions/workflows/Literate.yml/badge.svg)](https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche/actions/workflows/Literate.yml)

<!-- [![Build Status](https://github.com/omlins/ParallelStencil.jl/workflows/CI/badge.svg)](https://github.com/omlins/ParallelStencil.jl/actions) -->

### Porous convection 2D
This is the temperature distribution and flux when running the porous convection 2D code on a 511 1023 with 4000 timesteps
![Fig1](docs/PorousConvection2D.gif)


### Porous convection 3D
this is the same phenomenon in 3D. Here is the temperature distribution after 2000 timesteps with `nx,ny,nz = 255,127,127`.

## Discussion/Conclusion
