using Plots,Plots.Measures,Printf
default(size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function porous_convection_2D()
    # physics
    lx      = 20.0
    k = 1.
    η = 1.
    ρ       = (lx/(dc*2π))^2
    # numerics
    nx      = 100
    ϵtol    = 1e-8
    maxiter = 100nx
    ncheck  = ceil(Int,0.25nx)
    re = 4π
    cfl = 1/√
    # derived numerics
    dx      = lx/nx
    xc      = LinRange(dx/2,lx-dx/2,nx)
    dτ      = dx/sqrt(1/ρ)
    k_ηnf = 
    θ_dt =max(lx,ly)/re/cfl/max()
    β_dt = (re*k_ηnf)
    # array initialisation
    P       = @. 1.0 + exp(-(xc-lx/4)^2) - xc/lx; C_i = copy(P)
    qx      = zeros(Float64, nx-1)
    # iteration loop
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err >= ϵtol && iter <= maxiter
        # qx         .-= dτ./(ρ*k/η .+ dτ).*(qx .+ k/η.*diff(P)./dx)
        # P[2:end-1] .-= dτ.*diff(qx)./dx

        qDx[2:end-1] .-= (qDx[2:end-1] .+ k_ηnf .* diff(Pf)./dx)./(θ_dt .+ 1)
        Pf .-= (diff(qDx)./dx)./β_dt

        if iter%ncheck == 0
            err = maximum(abs.(diff(dc.*diff(P)./dx)./dx))
            push!(iter_evo,iter/nx); push!(err_evo,err)
            p1 = plot(xc,[C_i,P];xlims=(0,lx), ylims=(-0.1,2.0),
                      xlabel="lx",ylabel="Concentration",title="iter/nx=$(round(iter/nx,sigdigits=3))")
            p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",
                      yscale=:log10,grid=true,markershape=:circle,markersize=10)
            display(plot(p1,p2;layout=(2,1)))
        end
        iter += 1
    end
end

porous_convection_2D()
