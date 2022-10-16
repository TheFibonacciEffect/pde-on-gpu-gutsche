using Plots,Plots.Measures,Printf
using ProgressMeter
default(size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function porous_convection_2D()
    # physics
    lx      = 20.0
    ly      = 20.0
    k = 1.
    η = 1.
    # numerics
    nx      = 100
    ny      = 101
    ϵtol    = 1e-8
    maxiter = 100nx
    ncheck  = ceil(Int,0.25nx)
    re = 4π #TODO reynolds number?
    cfl = 1/√2.1
    # derived numerics
    dx      = lx/nx
    dy      = ly/ny
    xc      = LinRange(dx/2,lx-dx/2,nx)
    k_ηnf = 1 #TODO
    θ_dt =max(lx,ly)/re/cfl/max(dx,dy)
    β_dt = (re*k_ηnf)/(cfl*dx*lx)
    # array initialisation
    Pf       = @. 1.0 + exp(-(xc-lx/4)^2) - xc/lx; Pfi = copy(Pf)
    qDx      = zeros(nx-1)
    # iteration loop
    iter = 1; err_Pf = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err_Pf >= ϵtol && iter <= maxiter
        # boundary conditions on the flux
        qDx[1]   = 0
        qDx[end] = 0
        # differential equation
        qDx .-= (qDx .+ k_ηnf .* diff(Pf)./dx)./(θ_dt .+ 1)
        Pf[2:end-1] .-= (diff(qDx)./dx)./β_dt

        if iter%ncheck == 0
            r_Pf = diff(qDx)./dx #fluid is incompressible
            err_Pf = maximum(abs.(r_Pf))
            push!(iter_evo,iter/nx); push!(err_evo,err_Pf)
        end
        iter += 1
    end
    p1 = plot(xc,[Pf,Pfi])
    p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",
              yscale=:log10,grid=true,markershape=:circle,markersize=10)
    display(plot(p1,p2;layout=(2,1)))
end

porous_convection_2D()
