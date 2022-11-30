# Path: lecture8/src/ex2/t5/weak_scaling.jl
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
using Plots, Printf, MPI, MAT
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
print_time(@__LINE__)

include("../t3/l8_diffusion_2D_perf_multixpu.jl")

t,me = diffusion_2D(; do_visu=false,do_save=false,ttot = 1)

# write to file
if me==0
    fname = "weak_scaling.dat"
    print("writing to file $(fname)...")
    fid = open(fname,"a") #append
    write(fid, string(ARGS[1],",",t,"\n"))
    close(fid)
end
