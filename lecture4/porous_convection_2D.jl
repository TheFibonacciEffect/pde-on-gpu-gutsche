using Plots,Plots.Measures,Printf
default(size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views avx(A) = A[1:end-1,:] .+ A[2:end,:]./2
@views avy(A) = A[:,1:end-1] .+ A[:,2:end]./2

@views function porous_convection_2D()
    # physics
    lx        = 40.0
    ly        = 20.0
    k_ηf      = 1.0
    re        = 2π
    αρgx,αρgy = 0.0,1.0
    αρg       = sqrt(αρgx^2+αρgy^2)
    ΔT        = 200.0
    ϕ         = 0.1
    Ra        = 1000.0
    λ_ρCp     = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    # numerics
    nx        = 127
    ny        = ceil(Int,nx*ly/lx)
    nt        = 500
    ϵtol      = 1e-8
    maxiter   = 100max(nx,ny)
    ncheck    = ceil(Int,0.25max(nx,ny))
    nvis      = 10
    cfl       = 1.0/sqrt(2.1)
    # derived numerics
    dx        = lx/nx
    dy        = ly/ny
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1
    xc        = LinRange(-lx/2+dx/2,lx/2-dx/2,nx)
    yc        = LinRange(-ly+dy/2,-dy/2,ny)
    θ_dτ      = max(lx,ly)/re/cfl/min(dx,dy)
    β_dτ      = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    # visualisation init
    st        = ceil(Int,nx/25)
    Xc, Yc    = [x for x=xc, y=yc], [y for x=xc,y=yc]
    Xp, Yp    = Xc[1:st:end,1:st:end], Yc[1:st:end,1:st:end]
    # array initialisation
    T         = @. ΔT/1.5*exp(-xc^2 - (yc'+ly/2)^2)
    T[:,1]   .= ΔT/2; T[:,end] .= -ΔT/2
    Pf        = zeros(nx,ny)
    r_Pf      = zeros(nx,ny)
    qDx       = zeros(nx+1, ny)
    qDxc      = zeros(nx, ny)
    qDx_p     = qDxc[1:st:end,1:st:end]
    qDy       = zeros(nx, ny+1)
    qDyc      = zeros(nx, ny)
    qDy_p     = qDyc[1:st:end,1:st:end]
    qDmag     = zeros(nx, ny)
    # iteration loop
    anim = @animate for it in 1:nt
        iter = 1; err_Pf = 2ϵtol;
        while err_Pf >= ϵtol && iter <= maxiter
            qDx[2:end-1,:] .-= (qDx[2:end-1,:] .+ k_ηf.*(diff(Pf,dims=1)./dx .- αρgx.*avx(T)))./(1.0 + θ_dτ)
            qDy[:,2:end-1] .-= (qDy[:,2:end-1] .+ k_ηf.*(diff(Pf,dims=2)./dy .- αρgy.*avy(T)))./(1.0 + θ_dτ)
            r_Pf            .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
            Pf             .-= r_Pf./β_dτ
            if iter%ncheck == 0
                err_Pf = maximum(abs.(r_Pf))
                #p1 = heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]), ylims=(yc[1],yc[end]), aspect_ratio=1,
                #        xlabel="lx",ylabel="ly",title="iter/nx=$(round(iter/nx,sigdigits=3))",c=:turbo)
                # @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            end
            iter += 1
        end
        # recalculate time step
        dt_adv = ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1
        dt     = min(dt_diff,dt_adv)
        # update temperatur
        T[2:end-1,2:end-1] .+= dt.*λ_ρCp.*(diff(diff(T[:,2:end-1],dims=1)./dx,dims=1)./dx .+ diff(diff(T[2:end-1,:],dims=2)./dy,dims=2)./dy) # diffusion
        #T[2:end,:] .-= dt./ϕ.* (max.(0.0,qDx[2:end-1,:]) .* diff(T,dims=1)./dx)
        #T[1:end-1,:] .-= dt./ϕ.* (min.(0.0,qDx[2:end-1,:]) .* diff(T,dims=1)./dx)
        #T[:,2:end] .-= dt./ϕ.* (max.(0.0,qDy[:,2:end-1]) .* diff(T,dims=2)./dy)
        #T[:,1:end-1] .-= dt./ϕ.* (min.(0.0,qDy[:,2:end-1]) .* diff(T,dims=2)./dy)
        T[2:end-1,2:end-1] .-= dt./ϕ .* (
                                max.(0.0,qDx[2:end-2,2:end-1]) .* diff(T[1:end-1,2:end-1],dims=1)./dx .+
                                min.(0.0,qDx[2:end-2,2:end-1]) .* diff(T[2:end,2:end-1],dims=1)./dx .+
                                max.(0.0,qDy[2:end-1,2:end-2]) .* diff(T[2:end-1,1:end-1],dims=2)./dy .+
                                min.(0.0,qDy[2:end-1,2:end-2]) .* diff(T[2:end-1,2:end],dims=2)./dy
        )
        #T[:,1]   .= ΔT/2; T[:,end] .= -ΔT/2
        #T[2:end-1,2:end-1] .-= dt./ϕ.*(avx(qDx[2:end-1,2:end-1].*diff(T[:,2:end-1],dims=1)./dx) .+ avy(qDy[2:end-1,2:end-1].*diff(T[2:end-1,:],dims=2)./dy)) # advection
        # update adiabatic boundary condition
        T[[1,end],:] .= T[[2,end-1],:]
        @printf("it = %d, iter/nx=%.1f, err_Pf=%1.3e\n",it,iter/nx,err_Pf)
        if it % nvis == 0
            qDxc  .= avx(qDx)
            qDyc  .= avy(qDy)
            qDmag .= sqrt.(qDxc.^2 .+ qDyc.^2)
            qDxc  ./= qDmag
            qDyc  ./= qDmag
            qDx_p = qDxc[1:st:end,1:st:end]
            qDy_p = qDyc[1:st:end,1:st:end]
            heatmap(xc,yc,T';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            display(quiver!(Xp[:], Yp[:], quiver=(qDx_p[:], qDy_p[:]), lw=0.5, c=:black))
        end
    end
    gif(anim,"docs/l4e1t5.gif";fps=50)
end

porous_convection_2D()
