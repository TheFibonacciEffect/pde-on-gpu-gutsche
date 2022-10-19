using Plots,Plots.Measures,Printf
using ProgressMeter
default(size=(600*2,600*2),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@doc "this is kind of ad hoc"
function slice(A, dimension,direction, quantity=1)
    if dimension==2
        quantity == 1 && direction=='l' && return A[:,1:end-1]
        quantity == 1 && direction=='r' && return A[:,2:end]
    elseif dimension==1
        quantity == 1 && direction=='l' && return A[1:end-1,:]
        quantity == 1 && direction=='r' && return A[2:end,:]
    else
        quantity == 2 && return A[2:end-1,:]
    end
    return Nothing        
end

@views avx(A) = (A[1:end-1,:] .+ A[2:end,:])./2
@views avy(A) = (A[:,1:end-1] .+ A[:,2:end])./2


@views function porous_convection_2D(bounary,nvis)
    # physics
    lx      = 40.0
    ly      = 20.0
    k_ηf = 1
    re = 2π
    αρgx,αρgy = 0.0,1.0
    αρg       = sqrt(αρgx^2+αρgy^2)
    ΔT        = 200.0
    ϕ         = 0.1
    Ra        = 100
    λ_ρCp     = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
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
    xc      = LinRange(-lx/2 + dx/2,lx/2-dx/2,nx)
    yc      = LinRange(-ly+dy/2,dy/2,ny)
    θ_dt =max(lx,ly)/re/cfl/max(dx,dy)
    β_dt = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1
    # array initialisation
    # inital conditions
    Pf       = @. exp(-xc^2 - (yc'+ly/2)^2); Pfi = copy(Pf)
    T         = @. ΔT*exp(-xc^2 - (yc'+ly/2)^2)
    T[:,1] .= ΔT/2; T[:,end] .= -ΔT/2
    # preallocations
    qDx      = zeros(nx-1,ny) #darcy flux
    qDy      = zeros(nx,ny-1)
    qTx      = zeros(nx-1,ny)
    qTy      = zeros(nx,ny-1)

    tit = 500
    itvis = tit ÷ nvis
    for it=1:tit
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
            qDx .-= (qDx .+ k_ηf .* diff(Pf,dims=1)./dx .- αρgx.*avx(T))./(θ_dt .+ 1)
            qDy .-= (qDy .+ k_ηf .* diff(Pf,dims=2)./dy .- αρgy.*avy(T))./(θ_dt .+ 1)
            Pf[2:end-1,2:end-1] .-= (diff(qDx[:,2:end-1],dims=1)./dx + diff(qDy[2:end-1,:],dims=2)./dy)./β_dt

            if iter%ncheck == 0
                r_Pf = diff(qDx[:,2:end-1],dims=1)./dx + diff(qDy[2:end-1,:],dims=2)./dy #fluid is incompressible (continuity equation)
                err_Pf = maximum(abs.(r_Pf))
                push!(iter_evo,iter/nx); push!(err_evo,err_Pf)
            end
            iter += 1
        end
        dt_adv = ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1
        dt     = min(dt_diff,dt_adv)
        
        ∇xT = diff(T,dims=1)./dx
        ∇yT = diff(T,dims=2)./dy
        
        # T[2:end-1,2:end-1] .-= dt./ϕ.*(qDx.*∇xT + qDy.*∇yT) + dt.*λ_ρCp .* ∇²T
        # eq. 7
        # advection
        # max operator for upwind
        T[2:end-1,2:end-1] .-= dt./ϕ.*(
              (max.(qDx,0).*∇xT)[2:end,2:end-1]
            + (min.(qDx,0).*∇xT)[1:end-1,2:end-1]
            + (max.(qDy,0).*∇yT)[2:end-1,2:end]
            + (min.(qDy,0).*∇yT)[2:end-1,1:end-1]
            )
        
        # diffusion
        ∂²xT = diff(diff(T,dims=1),dims=1)[:,2:end-1]./dx./dx
        ∂²yT = diff(diff(T,dims=2),dims=2)[2:end-1,:]./dy./dy
        ∇²T = ∂²xT + ∂²yT
        
        T[2:end-1,2:end-1] .+= dt.*λ_ρCp .* ∇²T
        
        T[[1,end],:] .= T[[2,end-1],:]
        if it % itvis == 0
            @printf("it = %d, iter/nx=%.1f, err_Pf=%1.3e\n",it,iter/nx,err_Pf)
            st    = ceil(Int,nx/25)
            qDxc  = avy(qDx)
            qDyc  = avx(qDy)
            qDmag = sqrt.(qDxc.^2 .+ qDyc.^2)
            qDxc  ./= qDmag
            qDyc  ./= qDmag
            # qDx_p = qDxc[1:st:end,1:st:end]
            # qDy_p = qDyc[1:st:end,1:st:end]
            Xp = xc .* ones(size(yc))'
            Yp = ones(size(xc)) .* yc'
            @infiltrate
            p1 = heatmap(xc,yc,Pf',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            p2 = heatmap(xc,yc,T',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            p3 = heatmap(xc,yc,qDmag',xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)
            quiver!(p3,Xp[1:st:end,1:st:end], Yp[1:st:end,1:st:end], quiver=(qDxc[1:st:end,1:st:end], qDyc[1:st:end,1:st:end]), lw=0.5, c=:black)
            display(plot(p1,p2,p3,layout=(3,1)))
        end
    end
    # p2 = plot(iter_evo,err_evo;xlabel="iter/nx",ylabel="err",yscale=:log10,grid=true,markershape=:circle,markersize=10)
end

porous_convection_2D(false,10)
