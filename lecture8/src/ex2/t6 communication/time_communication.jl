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


# Physics
Lx, Ly  = 10.0, 10.0
D       = 1.0
# Numerics
nx, ny  = 64, 64 # number of grid points
nout    = 20
# Derived numerics
me, dims = init_global_grid(nx, ny, 1)  # Initialization more...

ENV["GKSwstype"]="nul"
include("l8_diffusion_2D_perf_multixpu.jl")

times = []
is = []
js = []


t_toc, me = diffusion_2D(Lx, Ly,D,nx, ny,nout,me, dims; do_visu=false,do_save=true,hidecom=false)
me == 0 && open("diffusion_2D_perf_multixpu.txt", "a") do f
    println(f, "no hidecomm, t_toc = $(t_toc)")
    push!(times,t_toc)
    push!(is, 0)
    push!(js, 0)
end

# ([2,16,16],[2,4,16])
for (i,j) in ([2,2],[8,2],[16,4],[16,16])
    t_toc, _ = diffusion_2D(Lx, Ly,D,nx, ny,nout,me, dims; do_visu=false,do_save=true,hc_x=i, hc_y=j)
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
finalize_global_grid()

if me==0
    println(is,times)
    # Plot
    p = plot(is, times./times[1],markershape=:circle, label="time hidecomm", xlabel="j", ylabel="time (s)/t_no_hidecom", title="time hidecomunication",xticks=is)

    println("plot:")
    println(p)
    # Save
    savefig(p,"time_communtication.png")
end