using Plots,Plots.Measures,Printf
default(size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function porous_convection_1D()
    # physics
    lx      = 40.0
    ly      = 20.0
    k_ηf    = 1.0
    re      = 2π
    ρg      = 1.0
    # numerics
    nx      = 100
    ny      = ceil(Int,nx*ly/lx)
    ϵtol    = 1e-8
    maxiter = 100max(nx,ny)
    ncheck  = ceil(Int,0.25max(nx,ny))
    cfl     = 1.0/sqrt(2.1)
    # derived numerics
    dx      = lx/nx
    dy      = ly/ny
    xc      = LinRange(dx/2,lx-dx/2,nx)
    yc      = LinRange(dy/2,ly-dy/2,ny)
    θ_dτ    = max(lx,ly)/re/cfl/min(dx,dy)
    β_dτ    = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    # array initialisation
    Pf       = @. exp(-(xc-lx/4)^2 -(yc'-ly/4)^2)
    r_Pf     = zeros(nx,ny)
    qDx      = zeros(nx+1, ny)
    qDy      = zeros(nx, ny+1)
    # iteration loop
    iter = 1; err_Pf = 2ϵtol;
    while err_Pf >= ϵtol && iter <= maxiter
        qDx[2:end-1,:] .-= (qDx[2:end-1,:] .+ k_ηf.*diff(Pf,dims=1)./dx)./(1.0 + θ_dτ)
        qDy[:,2:end-1] .-= (qDy[:,2:end-1] .+ k_ηf.*diff(Pf,dims=2)./dy .+ ρg)./(1.0 + θ_dτ)
        r_Pf         .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
        Pf           .-= r_Pf./β_dτ
        if iter%ncheck == 0
            err_Pf = maximum(abs.(r_Pf))
            p1 = heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]), ylims=(yc[1],yc[end]), aspect_ratio=1,
                      xlabel="lx",ylabel="ly",title="iter/nx=$(round(iter/nx,sigdigits=3))",c=:turbo)
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
        end
        iter += 1
    end
end

porous_convection_1D()
