# juliap -O3 --check-bounds=no --math-mode=fast diffusion_2D_perf_gpu.jl

const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots,Plots.Measures,Printf

@parallel function compute_flux!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ)
    @inn_x(qDx) = @inn_x(qDx) - (@inn_x(qDx) + k_ηf_dx*@d_xa(Pf))*_1_θ_dτ
    @inn_y(qDy) = @inn_y(qDy) - (@inn_y(qDy) + k_ηf_dy*@d_ya(Pf))*_1_θ_dτ
    return nothing
end

@parallel function update_Pf!(Pf,qDx,qDy,_dx,_dy,_β_dτ)
    @all(Pf) = @all(PF) - (@d_xa(qDx)*_dx + @d_ya(qDy)*_dy)*_β_dτ
    return nothing
end

@views function Pf_diffusion_2D(; do_visu=false)
    # Physics
    lx, ly  = 20.0, 20.0
    k_ηf    = 1.0
    # Numerics
    ϵtol    = 1e-8
    cfl     = s
    nx, ny  = 16*32, 16*32 # number of grid points
    nout    = 50
    maxiter = 500
    # Derived numerics
    dx, dy  = Lx/nx, Ly/ny
    dt      = min(dx, dy)^2/D/4.1
    nt      = cld(ttot, dt)
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    D_dx    = D/dx
    D_dy    = D/dy
    _dx, _dy= 1.0/dx, 1.0/dy
    # Array initialisation
    Pf      = Data.Array( @. exp(-(xc-lx/2)^2 -(yc'-ly/2)^2) )
    qDx     = @zeros(nx+1,ny  )
    qDy     = @zeros(nx  ,ny+1)
    r_Pf    = @zeros(nx  ,ny  )    size_C1_2, size_C2_2 = size(Pf,1)-2, size(Pf,2)-2
    t_tic = 0.0; niter = 0
    err_Pf = 2ϵtol
    # Time loop
    while err_Pf >= ϵtol && iter <= maxiter
        if (iter==11) t_tic = Base.time(); niter = 0 end
        @parallel compute_flux!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ)
        @parallel update_Pf!(Pf,qDx,qDy,_dx,_dy,_β_dτ)        Pf, Pf = Pf2, Pf # pointer swap
        niter += 1
        if do_visu && (it % nout == 0)
            opts = (aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), clims=(0.0, 1.0), c=:davos, xlabel="Lx", ylabel="Ly", title="time = $(round(it*dt, sigdigits=3))")
            display(heatmap(xc, yc, Array(Pf)'; opts...))
        end
        iter += 1; niter += 1
    end
    t_toc = Base.time() - t_tic
    A_eff = 2/1e9*nx*ny*sizeof(Float64)  # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                  # Execution time per iteration [s]
    T_eff = A_eff/t_it                   # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec, T_eff = %1.2f GB/s (niter = %d)\n", t_toc, round(T_eff, sigdigits=3), niter)
    return
end

Pf_diffusion_2D(; do_visu=false)