print("run_l8_diffusion_2D_perf_multixpu.jl")
print("beginning Execution")
t0 = time()
function print_time(line)
    println("on line $line Elapsed time: ", time() - t0, " seconds")
end
print_time(@__LINE__)
# juliap -O3 --check-bounds=no --math-mode=fast diffusion_2D_perf_xpu.jl
const USE_GPU = true
using ParallelStencil, ImplicitGlobalGrid
using ParallelStencil.FiniteDifferences2D
print_time(@__LINE__)
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
print_time(@__LINE__)
using Plots, Printf, MPI, MAT

include("l8_diffusion_2D_perf_multixpu.jl")

porous_convection_3D(;nz=63,do_viz=false)
