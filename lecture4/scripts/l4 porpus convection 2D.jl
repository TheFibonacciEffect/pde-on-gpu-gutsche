using Plots,Plots.Measures,Printf
using ProgressMeter
default(size=(600*2,600*3),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views avx(A) = (A[1:end-1,:] .+ A[2:end,:])./2
@views avy(A) = (A[:,1:end-1] .+ A[:,2:end])./2


@views function porous_convection_2D(bounary,nvis,timesteps)
    # physics
    lx      = 40.0
    ly      = 20.0
    k_ηf = 1
    re = 2π
    αρgx,αρgy = 0.0,1.0
    αρg       = sqrt(αρgx^2+αρgy^2)
    ΔT        = 200.0
    ϕ         = 0.1
    Ra        = 1000
    λ_ρCp     = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    # numerics
    nx      = 100
    ny      = 50
    ϵtol    = 1e-8
    maxiter   = 100max(nx,ny)
    ncheck    = ceil(Int,0.25max(nx,ny))
    cfl = 1/√2.1
    # derived numerics
    dx        = lx/nx
    dy        = ly/ny
    xc        = LinRange(-lx/2 + dx/2,lx/2-dx/2,nx)
    yc        = LinRange(-ly+dy/2,dy/2,ny)
    θ_dt      = max(lx,ly)/re/cfl/min(dx,dy)
    # θ_dt      = max(lx,ly)/re/cfl/max(dx,dy)
    β_dt      = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1
    # array initialisation
    # inital conditions
    Pf       = zeros(nx,ny)
    T         = @. ΔT*exp(-xc^2 - (yc'+ly/2)^2)
    T[:,1] .= ΔT/2; T[:,end] .= -ΔT/2
    # preallocations
    #darcy flux
    st         = ceil(Int,nx/25)
    qDx        = zeros(nx+1, ny)
    qDx_center = zeros(nx, ny)
    qDx_plot   = qDx_center[1:st:end,1:st:end]
    qDy        = zeros(nx, ny+1)
    qDy_center = zeros(nx, ny)
    qDy_plot   = qDy_center[1:st:end,1:st:end]
    qD_mag     = zeros(nx, ny)

    nvis = ceil(Int,timesteps/ nvis)
    anim = @animate for it=1:timesteps
        # iteration loop
        iter = 1; err_Pf = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        while err_Pf >= ϵtol && iter <= maxiter
            # println("($it,$iter)")
            # boundary conditions on the flux
            if bounary
                qDx[1,:]   .= 0
                qDx[end,:] .= 0
                qDx[:,1]   .= 0
                qDx[:,end] .= 0
                qDy[1,:]   .= 0
                qDy[end,:] .= 0
                qDy[:,1]   .= 0
                qDy[:,end] .= 0
            end
            # diffusion equation
            # qDx[2:end-1,:] .-= (qDx[2:end-1,:] .+ k_ηf .* diff(Pf,dims=1)./dx .- αρgx.*avx(T))./(θ_dt .+ 1)
            # qDy[:,2:end-1] .-= (qDy[:,2:end-1] .+ k_ηf .* diff(Pf,dims=2)./dy .- αρgy.*avy(T))./(θ_dt .+ 1)
            # Pf .-= (diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy)./β_dt

            qDx[2:end-1,:] .-= (qDx[2:end-1,:] .+ k_ηf.*(diff(Pf,dims=1)./dx .- αρgx.*avx(T)))./(1.0 + θ_dt)
            qDy[:,2:end-1] .-= (qDy[:,2:end-1] .+ k_ηf.*(diff(Pf,dims=2)./dy .- αρgy.*avy(T)))./(1.0 + θ_dt)
            r_Pf            = diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
            Pf             .-= r_Pf./β_dt

            if iter%ncheck == 0
                r_Pf = diff(qDx,dims=1)./dx + diff(qDy,dims=2)./dy #fluid is incompressible (continuity equation)
                err_Pf = maximum(abs.(r_Pf))
                push!(iter_evo,iter/nx); push!(err_evo,err_Pf)
            end
            iter += 1
        end
        dt_adv = ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1
        dt     = min(dt_diff,dt_adv)
        
        ∇xT = diff(T,dims=1)./dx
        ∇yT = diff(T,dims=2)./dy
        
        # diffusion
        ∂²xT = diff(diff(T[:,2:end-1],dims=1)./dx,dims=1)./dx
        ∂²yT = diff(diff(T[2:end-1,:],dims=2)./dy,dims=2)./dy
        T[2:end-1,2:end-1] .+= dt.*λ_ρCp.*(∂²xT .+ ∂²yT)

        # T[2:end-1,2:end-1] .-= dt./ϕ.*(qDx.*∇xT + qDy.*∇yT) + dt.*λ_ρCp .* ∇²T
        # eq. 7
        # advection
        # max operator for upwind
        T[2:end-1,2:end-1] .-= dt./ϕ .* (
                                max.(0.0,qDx[2:end-2,2:end-1]) .* diff(T[1:end-1,2:end-1],dims=1)./dx .+
                                min.(0.0,qDx[2:end-2,2:end-1]) .* diff(T[2:end,2:end-1],dims=1)./dx .+
                                max.(0.0,qDy[2:end-1,2:end-2]) .* diff(T[2:end-1,1:end-1],dims=2)./dy .+
                                min.(0.0,qDy[2:end-1,2:end-2]) .* diff(T[2:end-1,2:end],dims=2)./dy
        )        
        # diffusion
        # ∂²xT = diff(diff(T,dims=1),dims=1)[:,2:end-1]./dx./dx
        # ∂²yT = diff(diff(T,dims=2),dims=2)[2:end-1,:]./dy./dy
        # ∇²T = ∂²xT + ∂²yT


        # T[2:end-1,2:end-1] .+= dt.*λ_ρCp .* ∇²T
        
        # update bounary condition
        T[[1,end],:] .= T[[2,end-1],:]
        if it % nvis == 0
            @printf("it = %d, iter/nx=%.1f, err_Pf=%1.3e\n",it,iter/nx,err_Pf)
            qDx_center  .= avx(qDx)
            qDy_center  .= avy(qDy)
            qD_mag .= sqrt.(qDx_center.^2 .+ qDy_center.^2)
            qDx_center  ./= qD_mag
            qDy_center  ./= qD_mag
            qDx_plot .= qDx_center[1:st:end,1:st:end]
            qDy_plot .= qDy_center[1:st:end,1:st:end]

            Xp = xc .* ones(size(yc))'
            Yp = ones(size(xc)) .* yc'
            p1 = heatmap(xc,yc,Pf',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            title!("pressure at $it")
            p2 = heatmap(xc,yc,T',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            title!("temperature at $it")
            p3 = heatmap(xc,yc,qD_mag',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            title!("flux at $it")
            quiver!(p3,Xp[1:st:end,1:st:end], Yp[1:st:end,1:st:end], quiver=(qDx_center[1:st:end,1:st:end], qDy_center[1:st:end,1:st:end]), lw=0.5, c=:black)
            p = plot(p1,p2,p3,layout=(3,1))
            display(p)
        end
    end every nvis
    # p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",yscale=:log10,grid=true,markershape=:circle,markersize=10)
    return anim
end
a = porous_convection_2D(false,20,500)
gif(a,"figs/l4e1t4.gif",fps=15)
