using Plots,Plots.Measures,Printf
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

macro d_xa(A)  esc(:( $A[ix+1,iy]-$A[ix,iy]) ) end
macro d_ya(A)  esc(:( $A[ix,iy+1]-$A[ix,iy] )) end

@inbounds function compute_flux!(qDx, k_ηf, Pf, _dx, _1_θ_dτ, qDy, _dy)
    nx,ny=size(Pf)
    Threads.@threads for iy=1:ny
        for ix=1:nx-1
            qDx[ix+1,iy] -= (qDx[ix+1,iy] + k_ηf*((@d_xa(Pf)*_dx)))*_1_θ_dτ
        end
    end
    Threads.@threads for iy=1:ny-1
        for ix=1:nx
            qDy[ix,iy+1] -= (qDy[ix,iy+1] + k_ηf*((@d_ya(Pf))*_dy))*_1_θ_dτ
        end
    end
end
    
@inbounds function update_Pf!(Pf, qDx, _dx, qDy, _dy, _β_dτ)
    nx,ny=size(Pf)
    Threads.@threads for iy=1:ny
        for ix=1:nx
            Pf[ix,iy]  -= ((qDx[ix+1,iy]-qDx[ix,iy])*_dx + (qDy[ix,iy+1]-qDy[ix,iy])*_dy)*_β_dτ
        end
    end
end

function compute!(qDx, k_ηf, Pf, _dx, _1_θ_dτ, qDy, _dy, _β_dτ)
    compute_flux!(qDx, k_ηf, Pf, _dx, _1_θ_dτ, qDy, _dy)
    update_Pf!(Pf, qDx, _dx, qDy, _dy, _β_dτ)
    return nothing
end

function Pf_diffusion_2D_loop_fun_parralel(;do_check=false)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    # nx,ny   = 127,127
    nx,ny   = 1023,1023
    ϵtol    = 1e-8
    maxiter = 1e2#10max(nx,ny)
    ncheck  = 10#ceil(Int,0.25max(nx,ny))
    cfl     = 1.0/sqrt(2.1)
    re      = 2π
    # derived numerics
    dx,dy   = lx/nx,ly/ny
    xc,yc   = LinRange(dx/2,lx-dx/2,nx),LinRange(dy/2,ly-dy/2,ny)
    θ_dτ    = max(lx,ly)/re/cfl/min(dx,dy)
    β_dτ    = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    # array initialisation
    Pf      = @. exp(-(xc-lx/2)^2 -(yc'-ly/2)^2)
    qDx,qDy = zeros(Float64, nx+1,ny),zeros(Float64, nx,ny+1)
    r_Pf    = zeros(nx,ny)
    k_ηf_dx, k_ηf_dy = k_ηf/dx, k_ηf/dy
    _1_θ_dτ = 1.0./(1.0 + θ_dτ)
    _dx, _dy = 1.0/dx, 1.0/dy
    _β_dτ = 1.0/β_dτ
    
    # performance evaluation
    t_tic = 0.0
    # iteration loop
    iter = 1; err_Pf = 2ϵtol
    t_tic = 0.0; niter = 0
    while err_Pf >= ϵtol && iter <= maxiter
        if iter == 11 t_tic = Base.time(); niter = 0; end
        compute!(qDx, k_ηf, Pf, _dx, _1_θ_dτ, qDy, _dy, _β_dτ)
        if do_check && iter%ncheck == 0
            r_Pf  .= diff(qDx,dims=1).*_dx .+ diff(qDy,dims=2).*_dy
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            display(heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        end
        iter += 1; niter += 1
    end
    niter = iter-11
    t_toc = Base.time() - t_tic
    # read and write qDx, qDy and Pf --> 3 * read-write --> 3*2
    A_eff = 3*2*nx*ny*sizeof(eltype(Pf))/1e9   # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                        # Execution time per iteration [s]
    T_eff = A_eff/t_it                         # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec \n", t_toc)
    @printf("T_eff = %1.3f GB/sec \n", T_eff)
    @printf("niter = %1.3f \n", niter)
    return
end

Pf_diffusion_2D_loop_fun_parralel(;do_check=false)
