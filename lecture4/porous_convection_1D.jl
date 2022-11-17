using Plots,Plots.Measures,Printf
default(size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function porous_convection_1D()
    # physics
    lx      = 20.0
    ly      = 20.0
    k_ηf    = 1.0
    re      = 2π
    # numerics
    nx      = 100
    ϵtol    = 1e-8
    maxiter = 100nx
    ncheck  = ceil(Int,0.25nx)
    cfl     = 1.0/sqrt(2.1)
    # derived numerics
    dx      = lx/nx
    xc      = LinRange(dx/2,lx-dx/2,nx)
    θ_dτ    = lx/re/cfl/dx
    β_dτ    = (re*k_ηf)/(cfl*dx*lx)
    # array initialisation
    Pf       = @. exp(-(xc-lx/4)^2); Pf_0 = copy(Pf)
    r_Pf     = zeros(nx)
    qDx      = zeros(Float64, nx+1)
    # iteration loop
    iter = 1; err_Pf = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err_Pf >= ϵtol && iter <= maxiter
        qDx[2:end-1] .-= (qDx[2:end-1] .+ k_ηf.*diff(Pf)./dx)./(1.0 + θ_dτ)
        r_Pf         .= diff(qDx)./dx
        Pf           .-= r_Pf./β_dτ
        if iter%ncheck == 0
            err_Pf = maximum(abs.(diff(k_ηf.*diff(Pf)./dx)./dx))
            push!(iter_evo,iter/nx); push!(err_evo,err_Pf)
            p1 = plot(xc,[Pf_0,Pf];xlims=(0,lx), ylims=(-0.1,2.0),
                      xlabel="lx",ylabel="Pressure",title="iter/nx=$(round(iter/nx,sigdigits=3))")
            p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err_Pf",
                      yscale=:log10,grid=true,markershape=:circle,markersize=10)
            display(plot(p1,p2;layout=(2,1)))
        end
        iter += 1
    end
end

porous_convection_1D()
