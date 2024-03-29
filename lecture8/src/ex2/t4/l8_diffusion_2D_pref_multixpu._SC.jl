# juliap -O3 --check-bounds=no --math-mode=fast diffusion_2D_perf_xpu.jl
print("strong scaling")
print("l8_diffusion_2D_pref_multixpu._SC.jl")
const USE_GPU = false
using ParallelStencil, ImplicitGlobalGrid
using ParallelStencil.FiniteDifferences2D

@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end

using Plots, Printf, MPI, MAT


# macros to avoid array allocation
macro qx(ix,iy)  esc(:( -D_dx*(C[$ix+1,$iy+1] - C[$ix,$iy+1]) )) end
macro qy(ix,iy)  esc(:( -D_dy*(C[$ix+1,$iy+1] - C[$ix+1,$iy]) )) end

@parallel_indices (ix,iy) function compute!(C2, C, D_dx, D_dy, dt, _dx, _dy, size_C1_2, size_C2_2)
    if (ix<=size_C1_2 && iy<=size_C2_2)
        C2[ix+1,iy+1] = C[ix+1,iy+1] - dt*( (@qx(ix+1,iy) - @qx(ix,iy))*_dx + (@qy(ix,iy+1) - @qy(ix,iy))*_dy )
    end
    return
end

@views function diffusion_2D(nx, ny; do_visu=false,do_save=false,do_init_MPI=true)
    # Physics
    Lx, Ly  = 10.0, 10.0
    D       = 1.0
    ttot    = 1e-4
    # Numerics
    nout    = 20
    # Derived numerics
    me, dims = init_global_grid(nx, ny, 1, init_MPI=do_init_MPI)  # Initialization more...
    dx, dy  = Lx/nx_g(), Ly/ny_g()
    dt      = min(dx, dy)^2/D/4.1
    nt      = cld(ttot, dt)
    D_dx    = D/dx
    D_dy    = D/dy
    _dx, _dy= 1.0/dx, 1.0/dy
    # Array initialisation
    C       = @zeros(nx,ny)
    C      .= Data.Array([exp(-(x_g(ix,dx,C)+dx/2 -Lx/2)^2 -(y_g(iy,dy,C)+dy/2 -Ly/2)^2) for ix=1:size(C,1), iy=1:size(C,2)])   
    C2      = copy(C)
    size_C1_2, size_C2_2 = size(C,1)-2, size(C,2)-2
    t_tic = 0.0; niter = 0
    # Visualisation preparation
    if do_visu || do_save
        nx_v, ny_v = (nx-2)*dims[1], (ny-2)*dims[2]
        C_v   = zeros(nx_v, ny_v) # global array for visu and output
        C_inn = zeros(nx-2, ny-2) # no halo local array for visu
        if (nx_v*ny_v*sizeof(Data.Number) > 0.8*Sys.free_memory()) error("Not enough memory for visualization.") end
        # TODO This dir does not exist, but it doesn't matter because visualisation is not called. If I had more time I would fix this.
        do_visu && if (me==0) ENV["GKSwstype"]="nul"; if isdir("../docs/viz2D_mxpu_out")==false mkdir("../docs/viz2D_mxpu_out") end; loadpath = "../docs/viz2D_mxpu_out/"; anim = Animation(loadpath,String[]); println("Animation directory: $(anim.dir)") end
        xi_g, yi_g = LinRange(dx+dx/2, Lx-dx-dx/2, nx_v), LinRange(dy+dy/2, Ly-dy-dy/2, ny_v) # inner points only
    end
        # Time loop
    for it = 1:nt
        if (it==11) t_tic = Base.time(); niter = 0 end
        @hide_communication (8, 2) begin #with @hide_communication since it was the task description
            @parallel compute!(C2, C,_dx, _dy, D_dx, D_dy, dt, size_C1_2, size_C2_2)
            C, C2 = C2, C # pointer swap
            update_halo!(C)
        end
        niter += 1
        if do_visu && (it % nout == 0)
            C_inn .= Array(C)[2:end-1,2:end-1]; gather!(C_inn, C_v)
            if (me==0)
                opts = (aspect_ratio=1, xlims=(xi_g[1], xi_g[end]), ylims=(yi_g[1], yi_g[end]), clims=(0.0, 1.0), c=:turbo, xlabel="Lx", ylabel="Ly", title="time = $(round(it*dt, sigdigits=3))")
                heatmap(xi_g, yi_g, Array(C_v)'; opts...); frame(anim)
            end
        end
    end
    finalize_global_grid()
    # Create animation
    if (do_visu && me==0) gif(anim, "../docs/diffusion_2D_mxpu.gif", fps = 5)  end
    # Benchmarking
    t_toc = Base.time() - t_tic
    A_eff = 2/1e9*nx*ny*sizeof(Float64)  # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                  # Execution time per iteration [s]
    T_eff = A_eff/t_it                   # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec, T_eff = %1.2f GB/s (niter = %d)\n", t_toc, round(T_eff, sigdigits=3), niter)
    if (do_save && me==0)
        if isdir("../../../docs/l8ex2t3")==false mkdir("../../../docs/l8ex2t3") end
        file = matopen("../../../docs/l8ex2t3/mpigpu_out.mat", "w"); write(file, "C", Array(C_v)); close(file) 
    end
    return T_eff
end

@parallel function compute_flux!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ)
    @inn_x(qDx) = @inn_x(qDx) - (@inn_x(qDx) + k_ηf_dx*@d_xa(Pf))*_1_θ_dτ
    @inn_y(qDy) = @inn_y(qDy) - (@inn_y(qDy) + k_ηf_dy*@d_ya(Pf))*_1_θ_dτ
    return nothing
end

@parallel function update_Pf!(Pf,qDx,qDy,_dx,_dy,_β_dτ)
    @all(Pf) = @all(Pf) - (@d_xa(qDx)*_dx + @d_ya(qDy)*_dy)*_β_dτ
    return nothing
end


function Pf_diffusion_2D(nx=16*32,ny=16*32;do_check=false)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    threads = (32,4)
    ϵtol    = 1e-8
    maxiter = 500
    ncheck  = ceil(Int,0.25max(nx,ny))
    cfl     = 1.0/sqrt(2.1)
    re      = 2π
    # derived numerics
    dx,dy   = lx/nx,ly/ny
    xc,yc   = LinRange(dx/2,lx-dx/2,nx),LinRange(dy/2,ly-dy/2,ny)
    θ_dτ    = max(lx,ly)/re/cfl/min(dx,dy)
    β_dτ    = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    _1_θ_dτ = 1.0/(1.0 + θ_dτ)
    _β_dτ   = 1.0/(β_dτ)
    _dx,_dy = 1.0/dx,1.0/dy
    k_ηf_dx,k_ηf_dy = k_ηf/dx,k_ηf/dy
    # array initialisation
    Pf      = Data.Array( @. exp(-(xc-lx/2)^2 -(yc'-ly/2)^2) )
    qDx     = @zeros(nx+1,ny  )
    qDy     = @zeros(nx  ,ny+1)
    r_Pf    = @zeros(nx  ,ny  )    # visu
    if do_check
        ENV["GKSwstype"]="nul"
        if isdir("viz_out")==false mkdir("viz_out") end
        loadpath = "viz_out/"; anim = Animation(loadpath,String[])
        println("Animation directory: $(anim.dir)")
        iframe = 0
    end
    # iteration loop
    iter = 1; err_Pf = 2ϵtol
    t_tic = 0.0; niter = 0
    while err_Pf >= ϵtol && iter <= maxiter
        if (iter==11) t_tic = Base.time(); niter = 0 end
        @parallel compute_flux!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ)
        @parallel update_Pf!(Pf,qDx,qDy,_dx,_dy,_β_dτ)
        if do_check && (iter%ncheck == 0)
            r_Pf  .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            png((heatmap(xc,yc,Array(Pf)';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)),@sprintf("viz_out/%04d.png",iframe+=1))
        end
        iter += 1; niter += 1
    end
    t_toc = Base.time() - t_tic
    A_eff = (3*2)/1e9*nx*ny*sizeof(Float64)  # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                      # Execution time per iteration [s]
    T_eff = A_eff/t_it                       # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec, T_eff = %1.3f GB/s (niter = %d)\n", t_toc, round(T_eff, sigdigits=3), niter)
    return T_eff
end


function main()
    nx = ny = 16 * 2 .^ (1:10)
    Teff = []
    for i in eachindex(nx)
        push!(Teff,Pf_diffusion_2D(nx[i],ny[i]; do_check=false))
    end
    #plot the results
    plt = plot(nx,Teff,
        title = "Effective memory throughput Tesla P100 np = 1",
        xlabel="nx",
        ylabel="T_eff [GB/s]",
        marker = 2,
        markershape=:circle,
        markersize=10,
        xaxis=:log
    )
    savefig(plt,"../../../docs/StrongScaling.png")
end

main()