# Path: lecture8/src/ex2/t5/weak_scaling.jl
const USE_GPU = true
using ParallelStencil, ImplicitGlobalGrid
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots, Printf, MPI, MAT, CUDA
include("../t3/l8_diffusion_2D_perf_multixpu.jl")

t = diffusion_2D(; do_visu=false,do_save=false,ttot = 1e-5)

# write to file
fname = "weak_scaling.dat"
fid = open(fname,"a") #append
write(fid, string(length(CUDA.devices()),",",t,"\n")
close(fid)

