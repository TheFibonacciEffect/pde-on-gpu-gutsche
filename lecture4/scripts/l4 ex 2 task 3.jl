using Plots,Plots.Measures,Printf
using ProgressMeter
default(size=(600*2,600),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views avx(A) = (A[1:end-1,:] .+ A[2:end,:])./2
@views avy(A) = (A[:,1:end-1] .+ A[:,2:end])./2


function diffusion_eq!(qDx, k_ηf, Pf, dx, αρgx, T, θ_dt, qDy, dy, αρgy, β_dt)
    qDx[2:end-1,:] .-= (qDx[2:end-1,:] .+ k_ηf .* diff(Pf,dims=1)./dx .- αρgx.*avx(T))./(θ_dt .+ 1)
    qDy[:,2:end-1] .-= (qDy[:,2:end-1] .+ k_ηf .* diff(Pf,dims=2)./dy .- αρgy.*avy(T))./(θ_dt .+ 1)
    Pf .-= (diff(qDx,dims=1)./dx + diff(qDy,dims=2)./dy)./β_dt
end


function advection!(T, dt, ϕ, max, qDx, ∇xT, min, qDy, ∇yT)
    # max operator for upwind
    T[2:end-1,2:end-1] .-= dt./ϕ .* (
                            max.(0.0,qDx[2:end-2,2:end-1]) .* ∇xT[1:end-1,2:end-1] .+
                            min.(0.0,qDx[2:end-2,2:end-1]) .* ∇xT[2:end,2:end-1] .+
                            max.(0.0,qDy[2:end-1,2:end-2]) .* ∇yT[2:end-1,1:end-1] .+
                            min.(0.0,qDy[2:end-1,2:end-2]) .* ∇yT[2:end-1,2:end]
    )        
    
end


function temperature_diffusion!(T, dx, dy, dt, λ_ρCp)
    ∂²xT = diff(diff(T[:,2:end-1],dims=1)./dx,dims=1)./dx
    ∂²yT = diff(diff(T[2:end-1,:],dims=2)./dy,dims=2)./dy
    ∇²T = ∂²xT + ∂²yT
    T[2:end-1,2:end-1] .+= dt.*λ_ρCp.*∇²T # diffusion
end


function plot_result!(qDxc, qDx, qDyc, qDy, qDmag, sqrt, st, xc, yc, Pf, it, T)
    qDxc  .= avx(qDx)
    qDyc  .= avy(qDy)
    qDmag .= sqrt.(qDxc.^2 .+ qDyc.^2)
    qDxc  ./= qDmag
    qDyc  ./= qDmag
    qDx_p = qDxc[1:st:end,1:st:end]
    qDy_p = qDyc[1:st:end,1:st:end]
    
    Xp = xc .* ones(size(yc))'
    Yp = ones(size(xc)) .* yc'
    # p1 = heatmap(xc,yc,Pf',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
    # title!("pressure at $it")
    p = heatmap(xc,yc,T',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
    title!("temperature at $it")
    # p3 = heatmap(xc,yc,qDmag',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
    # title!("flux at $it")
    quiver!(p,Xp[1:st:end,1:st:end], Yp[1:st:end,1:st:end], quiver=(qDxc[1:st:end,1:st:end], qDyc[1:st:end,1:st:end]), lw=0.5, c=:black)
    display(p)    
end



@views function porous_convection_2D(Ra,nvis,timesteps;bounary=false)
    # physics
    lx      = 40.0
    ly      = 20.0
    k_ηf = 1
    re_D = 4π
    αρgx,αρgy = 0.0,1.0
    αρg       = sqrt(αρgx^2+αρgy^2)
    ΔT        = 200.0
    ϕ         = 0.1
    
    λ_ρCp     = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    # numerics
    nx      = 100
    ny      = 50
    ϵtol    = 1e-8
    maxiter = 100nx
    ncheck  = ceil(Int,0.25nx)
    cfl = 1/√2.1
    # derived numerics
    dx      = lx/nx
    dy      = ly/ny
    xc      = LinRange(-lx/2 + dx/2,lx/2-dx/2,nx)
    yc      = LinRange(-ly+dy/2,dy/2,ny)
    θ_dt_D =max(lx,ly)/re_D/cfl/max(dx,dy)
    β_dt_D = (re_D*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1
    # array initialisation
    # inital conditions
    Pf       = zeros(nx,ny)
    T         = @. ΔT*exp(-xc^2 - (yc'+ly/2)^2)
    T[:,1] .= ΔT/2; T[:,end] .= -ΔT/2
    # preallocations
    T_old = zeros(size(T))

    qDx      = zeros(nx+1,ny) #darcy flux
    qDy      = zeros(nx,ny+1)

    st    = ceil(Int,nx/25)
    qDx       = zeros(nx+1, ny)
    qDxc      = zeros(nx, ny)
    qDx_p     = qDxc[1:st:end,1:st:end]
    qDy       = zeros(nx, ny+1)
    qDyc      = zeros(nx, ny)
    qDy_p     = qDyc[1:st:end,1:st:end]
    qDmag     = zeros(nx, ny)

    dTdt        = zeros(nx-2,ny-2)
    r_T         = zeros(nx-2,ny-2)
    qTx         = zeros(nx-1,ny-2)
    qTy         = zeros(nx-2,ny-1)

    itvis = ceil(Int,timesteps/ nvis)
    anim = @animate for it=1:timesteps
        T_old .= T
        # time step
        dt = if it == 1
            0.1*min(dx,dy)/(αρg*ΔT*k_ηf)
        else
            min(0.1*min(dx,dy)/(αρg*ΔT*k_ηf),ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1)
        end

        re_T    = π + sqrt(π^2 + ly^2/λ_ρCp/dt)
        θ_dτ_T  = max(lx,ly)/re_T/cfl/min(dx,dy)
        β_dτ_T  = (re_T*λ_ρCp)/(cfl*min(dx,dy)*max(lx,ly))

        # iteration loop
        iter = 1; err_D = 2ϵtol; err_T = 2ϵtol
        while max(err_D,err_T) >= ϵtol && iter <= maxiter
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
            diffusion_eq!(qDx, k_ηf, Pf, dx, αρgx, T, θ_dt_D, qDy, dy, αρgy, β_dt_D)

            # temperature update
            # Temperatur
            qTx            .-= (qTx .+ λ_ρCp .* diff(T[:,2:end-1],dims=1)./dx)./(1.0 + θ_dτ_T)
            qTy            .-= (qTy .+ λ_ρCp .* diff(T[2:end-1,:],dims=2)./dy)./(1.0 + θ_dτ_T)
            dTdt           .= (T[2:end-1,2:end-1] .- T_old[2:end-1,2:end-1])./dt .+ (
                                    max.(0.0,qDx[2:end-2,2:end-1]) .* diff(T[1:end-1,2:end-1],dims=1)./dx .+
                                    min.(0.0,qDx[2:end-2,2:end-1]) .* diff(T[2:end,2:end-1],dims=1)./dx .+
                                    max.(0.0,qDy[2:end-1,2:end-2]) .* diff(T[2:end-1,1:end-1],dims=2)./dy .+
                                    min.(0.0,qDy[2:end-1,2:end-2]) .* diff(T[2:end-1,2:end],dims=2)./dy
            )./ϕ .+ diff(qTx,dims=1)./dx .+ diff(qTy,dims=2)./dy
            r_T .= (dTdt .+ diff(qTx,dims=1)./dx .+ diff(qTy,dims=2)./dy)
            T[2:end-1,2:end-1] .-= r_T./(1.0/dt + β_dτ_T )
            T[[1,end],:] .= T[[2,end-1],:]
            if iter%ncheck == 0
                r_Pf = diff(qDx,dims=1)./dx + diff(qDy,dims=2)./dy #fluid is incompressible (continuity equation)
                err_D = maximum(abs.(r_Pf))
                err_T  = maximum(abs.(r_T))
            end
            iter += 1
        end
        if it % itvis == 0
            @printf("it = %d, iter/nx=%.1f, err_D=%1.3e\n",it,iter/nx,err_D)
            plot_result!(qDxc, qDx, qDyc, qDy, qDmag, sqrt, st, xc, yc, Pf, it, T)
        end
    end every itvis
    return anim
end
for Ra = [10, 40, 100, 1000]
    ani = porous_convection_2D(Ra,10,250)
    gif(ani,"figs/task3/implicit at Ra=$Ra.gif",fps=15)
end
