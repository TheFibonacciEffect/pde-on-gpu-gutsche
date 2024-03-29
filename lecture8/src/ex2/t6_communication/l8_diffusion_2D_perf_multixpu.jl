# macros to avoid array allocation
macro qx(ix,iy)  esc(:( -D_dx*(C[$ix+1,$iy+1] - C[$ix,$iy+1]) )) end
macro qy(ix,iy)  esc(:( -D_dy*(C[$ix+1,$iy+1] - C[$ix+1,$iy]) )) end

@parallel_indices (ix,iy) function compute!(C2, C, D_dx, D_dy, dt, _dx, _dy, size_C1_2, size_C2_2)
    if (ix<=size_C1_2 && iy<=size_C2_2)
        C2[ix+1,iy+1] = C[ix+1,iy+1] - dt*( (@qx(ix+1,iy) - @qx(ix,iy))*_dx + (@qy(ix,iy+1) - @qy(ix,iy))*_dy )
    end
    return
end

@views function diffusion_2D(Lx, Ly,D,nx, ny,nout,me, dims; do_visu=false,do_save=false,ttot = 1,timestep_override=0,hc_x=8, hc_y=2,hidecom=true)
    dx, dy  = Lx/nx_g(), Ly/ny_g()
    dt      = min(dx, dy)^2/D/4.1
    nt      = timestep_override==0 ? cld(ttot, dt) : timestep_override
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
    GC.gc(); GC.enable(false)
    for it = 1:nt
        if (it==11) t_tic = Base.time(); niter = 0 end  #NOTE if the time is very small this is not reached
            if hidecom
                @hide_communication (hc_x, hc_y) begin #with @hide_communication since it was the task description
                    @parallel compute!(C2, C,_dx, _dy, D_dx, D_dy, dt, size_C1_2, size_C2_2)
                    C, C2 = C2, C # pointer swap
                    update_halo!(C)
                end
            else
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
    GC.enable(true)
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
    return t_toc, me
end