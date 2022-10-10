using Plots,Plots.Measures,Printf
#imports a nice progress bar
using ProgressMeter
default(framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function steady_diffusion_1D()
    # physics
    lx   = 10.0
    dc      = 1.0
    vx      = 10.0
    da      = 1000.0
    re      = π + sqrt(π^2 + da)
    ρ       = (lx/(dc*re))^2 #how does ly factor into this?
    # numerics
    nx   = 200
    ϵtol    = 1e-8
    maxiter = 10nx
    ncheck  = ceil(Int,0.02nx)
    nt      = 50
    # derived numerics
    dx      = lx/nx
    dt       = lx^2/dc/da
    xc      = LinRange(dx/2,lx-dx/2,nx)
    dτ      = dx/sqrt(1/ρ)
    # array initialisation
    C       = @. 1.0 + exp(-(xc-lx/4)^2) - xc/lx; C_i = copy(C)
    C_old   = copy(C) 
    qx      = zeros(nx-1)
    # other
    p = Progress(nt, 1)
    anim = @animate for it = 1:nt
        C_old .= C
        # why does this produce such a strange plot?
        # p1 = heatmap(xc,yc,C;xlabel="lx",ylabel="ly",title="Implicit transient diffusion using dual timestepping")
        p1 = plot(xc,C,title="please help! D:")
        # iteration loop
        iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        while err >= ϵtol && iter <= maxiter
            qx           .-= dτ./(ρ   .+ dτ/dc).*(qx./dc                                 .+ diff(C,dims=1) ./dx)
            diffxC = dτ./(1.0 .+ dτ/dt).*((C[2:end-1] .- C_old[2:end-1])./dt .+ diff(qx,dims=1)./dx)
            C[2:end-1] .-= diffxC
            println(sqrt.(sum((dτ./(1.0 .+ dτ/dt).*((C[2:end-1,:] .- C_old[2:end-1,:])./dt .+ diff(qx,dims=1)./dx)).^2)))
            if iter%ncheck == 0
                △xC = (diff(dc.*diff(C,dims=1)./dx,dims=1)./dx)
                err = maximum(abs.(△xC .- (C[2:end-1] .- C_old[2:end-1])./dt))
                push!(iter_evo,iter/nx); push!(err_evo,err)
            end
            iter += 1
        end
        p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",
                yscale=:log10,grid=true,markershape=:circle,markersize=10, title="error in peseudotimestep")
        display(plot(p1,p2;layout=(2,1)))
        next!(p)
    end
    anim
end

an =  steady_diffusion_1D()

gif(an,"figs/l3e1.gif",fps=2)
