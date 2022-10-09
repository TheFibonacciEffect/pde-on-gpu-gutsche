using Plots,Plots.Measures,Printf
#imports a nice progress bar
using ProgressMeter
default(size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function steady_diffusion_1D()
    # physics and nummerics
    vx = 1.0
    lx      = 20.0
    nx      = 100
    dx      = lx/nx
    dt = dx/abs(vx)
    dc      = 1.0
    da      = lx^2/dc/dt
    re      = π + sqrt(π^2 + da)
    ρ       = (lx/(dc*re))^2
    ϵtol    = 1e-8
    maxiter = 50nx
    ncheck  = ceil(Int,0.25nx)
    nt = 10 #number of physical timesteps
    xc      = LinRange(dx/2,lx-dx/2,nx)
    dτ      = dx/sqrt(1/ρ)
    # array initialisation
    C       = @. 1.0 + exp(-(xc-lx/4)^2) - xc/lx; C_i = copy(C)
    C_old    = copy(C)
    qx      = zeros(Float64, nx-1)
    # other
    p = Progress(nt, 1)
    for it = 1:nt
        C_old .= C
        # iteration loop
        iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        while err >= ϵtol && iter <= maxiter
            qx         .-= dτ./(ρ   .+ dτ/dc).*(qx./dc                  .+ diff(C) ./dx)
            C[2:end-1] .-= dτ./(1.0 .+ dτ/dt ).*((C[2:end-1] .- C_old[2:end-1])./dt .+ diff(qx)./dx)
            if iter%ncheck == 0
                err = maximum(abs.(diff(dc.*diff(C)./dx)./dx .- (C[2:end-1] .- C_old[2:end-1])./dt))
                push!(iter_evo,iter/nx); push!(err_evo,err)
            end
            iter += 1
        end
        
        next!(p)
    end
    plot(xc,[C_i,C];xlims=(0,lx), ylims=(-0.1,2.0),
            xlabel="lx",ylabel="Concentration",title="implicit_advection_diffusion_1D")
end

steady_diffusion_1D()
