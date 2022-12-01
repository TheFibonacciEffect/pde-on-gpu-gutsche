# juliap -O3 --check-bounds=no --math-mode=fast diffusion_2D_perf_xpu.jl
const USE_GPU = true
using ParallelStencil, ImplicitGlobalGrid
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots, Printf, MPI, MAT

include("../t3/l8_diffusion_2D_perf_multixpu.jl")

times = []
is = []
js = []
for i,j in ([2,16,16],[2,4,16])
    t_toc, me = diffusion_2D(; do_visu=false,do_save=true)
    if me == 0
        push!(times, t_toc)
        #TODO append to file
    end
end

# TODO: Plot