using Plots,Plots.Measures,Printf
using ProgressMeter
default(size=(600*2,600),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function porous_convection_2D(bounary)
    # physics
    lx      = 40.0
    ly      = 20.0
    k_ηf = 1
    re = 2π
    # numerics
    nx      = 200
    ny      = 100
    ϵtol    = 1e-8
    maxiter = 100nx
    ncheck  = ceil(Int,0.25nx)
    cfl = 1/√2.1
    # derived numerics
    dx      = lx/nx
    dy      = ly/ny
    xc      = LinRange(dx/2,lx-dx/2,nx)
    yc      = LinRange(dy/2,ly-dy/2,ny)
    θ_dt =max(lx,ly)/re/cfl/max(dx,dy)
    β_dt = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    # array initialisation
    Pf       = @. exp(-(xc-lx/4)^2-(yc'-ly/4)^2); Pfi = copy(Pf)
    qDx      = zeros(nx-1,ny)
    qDy      = zeros(nx,ny-1)
    # iteration loop
    iter = 1; err_Pf = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    for it=1:10
        while err_Pf >= ϵtol && iter <= maxiter
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
            # differential equation
            qDx .-= (qDx .+ k_ηf .* diff(Pf,dims=1)./dx)./(θ_dt .+ 1)
            qDy .-= (qDy .+ k_ηf .* diff(Pf,dims=2)./dy)./(θ_dt .+ 1)
            Pf[2:end-1,2:end-1] .-= (diff(qDx[:,2:end-1],dims=1)./dx + diff(qDy[2:end-1,:],dims=2)./dy)./β_dt

            if iter%ncheck == 0
                r_Pf = diff(qDx[:,2:end-1],dims=1)./dx + diff(qDy[2:end-1,:],dims=2)./dy #fluid is incompressible
                err_Pf = maximum(abs.(r_Pf))
                push!(iter_evo,iter/nx); push!(err_evo,err_Pf)
            end
            iter += 1
        end
        @printf("it = %d, iter/nx=%.1f, err_Pf=%1.3e\n",it,iter/nx,err_Pf)
        p = heatmap(xc,yc,Pf',c=:turbo)
        display(p)
    end
    # p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",yscale=:log10,grid=true,markershape=:circle,markersize=10)
end

porous_convection_2D(true)
