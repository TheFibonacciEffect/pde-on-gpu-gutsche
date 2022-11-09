using Printf,LazyArrays,Plots

const USE_GPU = true
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end

@views av1(A) = 0.5.*(A[1:end-1].+A[2:end])
@views avx(A) = 0.5.*(A[1:end-1,:].+A[2:end,:])
@views avy(A) = 0.5.*(A[:,1:end-1].+A[:,2:end])


@parallel function compute_flux!(qDx,qDy,Pf,T,k_ηf,_dx,_dy,αρgx,αρgy,θ_dτ_D)
    @inn_x(qDx) = @inn_x(qDx) - (@inn_x(qDx) + k_ηf*(@d_xa(Pf)*_dx - αρgx*@av_xa(T)))/(1.0 + θ_dτ_D)
    @inn_y(qDy) = @inn_y(qDy) - (@inn_y(qDy) + k_ηf*(@d_ya(Pf)*_dy - αρgy*@av_ya(T)))/(1.0 + θ_dτ_D)
    return nothing
end

@parallel function update_Pf!(Pf,qDx,qDy,_dx,_dy,β_dτ_D)
    @all(Pf) = @all(Pf) - (@d_xa(qDx)*_dx + @d_ya(qDy)*_dy)./β_dτ_D
    return nothing
end

@parallel function update_temperature_flux(qTx, λ_ρCp, T, _dx, θ_dτ_T, qTy, _dy)
    @all(qTx) = @all(qTx) - (@all(qTx) + λ_ρCp.*(@d_xi(T)*_dx))./(1.0 + θ_dτ_T)
    @all(qTy) = @all(qTy) - (@all(qTy) + λ_ρCp.*(@d_yi(T)*_dy))./(1.0 + θ_dτ_T)
    return
end

@parallel_indices (i,j) function update_dTdt!(dTdt,T, T_old, dt, max, qDx, _dx, min, qDy, _dy, ϕ)
            dTdt[i-1,j-1] =    (T[i,j] - T_old[i,j])/dt +
                           (max(qDx[i  ,j  ],0.0)*(T[i  ,j  ] - T[i-1,j  ])*_dx +
                            min(qDx[i+1,j  ],0.0)*(T[i+1,j  ] - T[i  ,j  ])*_dx +
                            max(qDy[i  ,j  ],0.0)*(T[i  ,j  ] - T[i  ,j-1])*_dy +
                            min(qDy[i  ,j+1],0.0)*(T[i  ,j+1] - T[i  ,j  ])*_dy)/ϕ
    return
end

@parallel function temperature_update!(T,dTdt,qTx,qTy,_dx,_dy,β_dτ_T,dt)
    @inn(T) = @inn(T)- ((@all(dTdt) + @d_xa(qTx) * _dx) + @d_ya(qTy) * _dy) / (1.0 / dt + β_dτ_T)
    return
end

@parallel_indices (iy) function bc_x!(A)
    A[1  ,iy] = A[2    ,iy]
    A[end,iy] = A[end-1,iy]
    return
end

@views function porous_convection_2D()
    # physics
    lx,ly       = 40.0,20.0
    k_ηf        = 1.0
    αρgx,αρgy   = 0.0,1.0
    αρg         = sqrt(αρgx^2+αρgy^2)
    ΔT          = 200.0
    ϕ           = 0.1
    Ra          = 1000
    λ_ρCp       = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    # numerics
    ny          = 1023
    nx          = 511
    nt          = 4000
    re_D        = 4π
    cfl         = 1.0/sqrt(2.1)
    maxiter     = 10max(nx,ny)
    ϵtol        = 1e-6
    nvis        = 20
    ncheck      = ceil(2*max(nx,ny))
    # preprocessing
    dx,dy       = lx/nx,ly/ny
    _dx,_dy       = 1.0/dx,1.0/dy
    xn,yn       = LinRange(-lx/2,lx/2,nx+1),LinRange(-ly,0,ny+1)
    xc,yc       = av1(xn),av1(yn)
    θ_dτ_D      = max(lx,ly)/re_D/cfl/min(dx,dy)
    β_dτ_D      = (re_D*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    # init
    Pf          = @zeros(nx,ny)
    r_Pf        = @zeros(nx,ny)
    qDx,qDy     = @zeros(nx+1,ny),@zeros(nx,ny+1)
    qDx_c,qDy_c = @zeros(nx,ny),@zeros(nx,ny)
    qDmag       = zeros(nx,ny)     
    T           = Data.Array(@. ΔT*exp(-xc^2 - (yc'+ly/2)^2))
    T[:,1]     .= ΔT/2; T[:,end] .= -ΔT/2
    T_old       = copy(T)
    dTdt        = @zeros(nx-2,ny-2)
    r_T         = @zeros(nx-2,ny-2)
    qTx         = @zeros(nx-1,ny-2)
    qTy         = @zeros(nx-2,ny-1)
    # vis
    st          = ceil(Int,nx/25)
    Xc, Yc      = [x for x=xc, y=yc], [y for x=xc,y=yc]
    Xp, Yp      = Xc[1:st:end,1:st:end], Yc[1:st:end,1:st:end]
    iframe = 0
    # action
    for it = 1:nt
        T_old .= T
        # time step
        dt = if it == 1
            0.1*min(dx,dy)/(αρg*ΔT*k_ηf)
        else
            min(5.0*min(dx,dy)/(αρg*ΔT*k_ηf),ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1)
        end
        re_T    = π + sqrt(π^2 + ly^2/λ_ρCp/dt)
        θ_dτ_T  = max(lx,ly)/re_T/cfl/min(dx,dy)
        β_dτ_T  = (re_T*λ_ρCp)/(cfl*min(dx,dy)*max(lx,ly))
        # iteration loop
        iter = 1; err_D = 2ϵtol; err_T = 2ϵtol
        while max(err_D,err_T) >= ϵtol && iter <= maxiter
            # hydro
            @parallel compute_flux!(qDx,qDy,Pf,T,k_ηf,_dx,_dy,αρgx,αρgy,θ_dτ_D)
            @parallel update_Pf!(Pf,qDx,qDy,_dx,_dy,β_dτ_D)
            # thermo
            @parallel update_temperature_flux(qTx, λ_ρCp, T, _dx, θ_dτ_T, qTy, _dy)
            @parallel (2:(size(T,1)-1), 2:(size(T,2)-1)) update_dTdt!(dTdt,T, T_old, dt, max, qDx, _dx, min, qDy, _dy, ϕ)
            @parallel temperature_update!(T,dTdt,qTx,qTy,_dx,_dy,β_dτ_T,dt)
            #T[[1,end],:]        .= T[[2,end-1],:]
            # Periodic boundary condition
            @parallel (1:size(T,2)) bc_x!(T)
            if iter % ncheck == 0
                r_Pf  .= Diff(qDx,dims=1)./dx .+ Diff(qDy,dims=2)./dy
                r_T   .= dTdt .+ Diff(qTx,dims=1)./dx .+ Diff(qTy,dims=2)./dy
                err_D  = maximum(abs.(r_Pf))
                err_T  = maximum(abs.(r_T))
                @printf("  iter/nx=%.1f, err_D=%1.3e, err_T=%1.3e\n",iter/nx,err_D,err_T)
            end
            iter += 1
        end
        @printf("it = %d, iter/nx=%.1f, err_D=%1.3e, err_T=%1.3e\n",it,iter/nx,err_D,err_T)
        # visualisation
        if it % nvis == 0
            qDx_c .= avx(Array(qDx))
            qDy_c .= avy(Array(qDy))
            qDmag .= sqrt.(qDx_c.^2 .+ qDy_c.^2)
            qDx_c ./= qDmag
            qDy_c ./= qDmag
            qDx_p = qDx_c[1:st:end,1:st:end]
            qDy_p = qDy_c[1:st:end,1:st:end]
            heatmap(xc,yc,T';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            display(quiver!(Xp[:], Yp[:], quiver=(qDx_p[:], qDy_p[:]), lw=0.5, c=:black))
            # save(@sprintf("anim/%04d.png",iframe),fig); iframe += 1
            if isdir("PorousConvection_2D_xpu_daint_out")==false mkdir("PorousConvection_2D_xpu_daint_out") end
            savefig("PorousConvection_2D_xpu_daint_out/PorousConvection_2D_xpu-$nx-$ny-t-$it.png")
        end
    end
    return
end

porous_convection_2D()
