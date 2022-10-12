using Plots,Plots.Measures,Printf
#imports a nice progress bar
using ProgressMeter
default(size=(800,1200), framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function implicit_advection_diffusion_2D()
    # visualisation option
    advection = true
    # physics
    lx,ly   = 10.0,10.0
    dc      = 1.0
    vx      = 1.0
    vy      = -1.0
    dc      = 1.0
    da      = 1000.0
    re      = π + sqrt(π^2 + da)
    ρ       = (lx/(dc*re))^2

    # numerics
    nx,ny   = 200,201
    ϵtol    = 1e-8
    maxiter = 10nx
    ncheck  = ceil(Int,0.02nx)
    nt      = 50
    # derived numerics
    dx      = lx/nx
    dy      = ly/ny
    dt      = min(dx/abs(vx),dy/abs(vy))/2
    xc      = LinRange(dx/2,lx-dx/2,nx)
    yc      = LinRange(dy/2,ly-dy/2,ny)
    dτ      = min(dx,dy)/sqrt(1/ρ)/sqrt(2)
    # array initialisation
    C       = @. exp(-(xc-lx/4)^2 -(yc'-3ly/4)^2)
    C_old   = copy(C)
    qx      = zeros(nx-1,ny)
    qy      = zeros(nx,ny-1)
    ∇q      = zeros(nx-2,ny-2)
    # iteration loop
    p = Progress(nt, 1)
    anim = @animate for it = 1:nt
        C_old .= C
        iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        while err >= ϵtol && iter <= maxiter
            qx .-= dτ./(ρ   .+ dτ/dc).*(qx./dc .+ diff(C,dims=1) ./dx)
            qy .-= dτ./(ρ   .+ dτ/dc).*(qy./dc .+ diff(C,dims=2) ./dy)
            ∇q  .= diff(qx[:,2:end-1],dims=1)./dx .+ diff(qy[2:end-1,:],dims=2)./dy
            C[2:end-1,2:end-1] .-= dτ./(1.0 .+ dτ/dt).*((C[2:end-1,2:end-1] .- C_old[2:end-1,2:end-1])./dt .+ ∇q)
            if iter%ncheck == 0
                △yC = diff(dc.*diff(C[2:end-1,:],dims=2)./dy,dims=2)./dy
                △xC = diff(dc.*diff(C[:,2:end-1],dims=1)./dx,dims=1)./dx
                err = maximum(abs.(△xC .+ △yC .- (C[2:end-1,2:end-1] .- C_old[2:end-1,2:end-1])./dt))
                push!(iter_evo,iter/nx); push!(err_evo,err)
            end
            iter += 1
        end
        # visualisation
        p1 = heatmap(xc,yc,C';xlims=(0,lx),ylims=(0,ly),clims=(0,1),aspect_ratio=1,
                 xlabel="lx",ylabel="ly",title="iter/nx=$(round(iter/nx,sigdigits=3))")
        p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",
                yscale=:log10,grid=true,markershape=:circle,markersize=10)
        display(plot(p1,p2;layout=(2,1)))
        next!(p)
        if advection
            # x direction
            C[1:end-1,:] .-= dt.* max(vx,0) .*diff(C,dims=1)./dx
            vx > 0 && (C[1,:] = C_old[1,:])
            C[2:end,:]   .-= dt.* min(vx,0) .*diff(C,dims=1)./dx
            vx < 0 && (C[end,:] = C_old[end,:])
            # y direction
            C[:,1:end-1] .-= dt.* max(vy,0) .*diff(C,dims=2)./dy
            vy > 0 && (C[:,1] = C_old[:,1])
            C[:,2:end]   .-= dt.* min(vy,0) .*diff(C,dims=2)./dy
            vy < 0 && (C[:,end] = C_old[:,end])
        end
    end
    gif(anim,"lecture3/figs/l3e3t2.gif",fps=2)
end

implicit_advection_diffusion_2D()