const USE_GPU = true
using ParallelStencil, ImplicitGlobalGrid
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end
using Plots, Printf, MPI, MAT

include("PorousConvection_3D_multixpu.jl")

porous_convection_3D(;nz=63,do_visu=true)