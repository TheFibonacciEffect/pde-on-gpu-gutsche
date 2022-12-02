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


t_toc, me = diffusion_2D(; do_visu=false,do_save=true,hidecom=false)
me == 0 && open("diffusion_2D_perf_multixpu.txt", "a") do f
    println(f, "no hidecomm, t_toc = $(t_toc)")
end

for (i,j) in ([2,16,16],[2,4,16])
    t_toc, me = diffusion_2D(; do_visu=false,do_save=true)
    if me == 0
        push!(times, t_toc)
        push!(is, i)
        push!(js, j)
        # append to file
        open("diffusion_2D_perf_multixpu.txt", "a") do f
            println(f, "i = $(i), j = $(j), t_toc = $(t_toc)")
        end
    end
end

# Plot
plot(is, times, label="time hidecomm", xlabel="i", ylabel="time (s)", title="time hidecomunication",xticks=[no-hidecomm, (2,2), (8,2), (16,4), (16,16)])

# Save
savefig("time_communtication.png")