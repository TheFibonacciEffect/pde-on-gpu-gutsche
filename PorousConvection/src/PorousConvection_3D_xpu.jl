
module PorousConvection_3D_xpu
using Printf,Plots
using Dates
    
const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences3D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end

@views av1(A) = 0.5.*(A[1:end-1].+A[2:end])

"""
    save_array(Aname,A)

Saves array A as binary file with name Aname.
"""
function save_array(Aname,A)
    fname = string(Aname,".bin")
    out = open(fname,"w"); write(out,A); close(out)
end

"""
    compute_Dflux!(qDx,qDy,qDz,Pf,T,k_ηf,_dx,_dy,_dz,αρg,_1_θ_dτ_D)

Calculated current Darcy flux from effective pressure and subtracts it from previous result.
"""
@parallel function compute_flux!(qDx,qDy,qDz,Pf,T,k_ηf,_dx,_dy,_dz,αρgx,αρgy,αρgz,_θ_dτ_D)
    @inn_x(qDx) = @inn_x(qDx) - (@inn_x(qDx) + k_ηf*(@d_xa(Pf)*_dx - αρgx*@av_xa(T)))*_θ_dτ_D
    @inn_y(qDy) = @inn_y(qDy) - (@inn_y(qDy) + k_ηf*(@d_ya(Pf)*_dy - αρgy*@av_ya(T)))*_θ_dτ_D
    @inn_z(qDz) = @inn_z(qDz) - (@inn_z(qDz) + k_ηf*(@d_za(Pf)*_dz - αρgz*@av_za(T)))*_θ_dτ_D
    return nothing
end

"""
    update_Pf!(Pf,qDx,qDy,qDz,_dx,_dy,_dz,_β_dτ_D)

Effective pressure update calculated from the divergence of the Darcy flux.
"""
@parallel function update_Pf!(Pf,qDx,qDy,qDz,_dx,_dy,_dz,_β_dτ_D)
    @all(Pf) = @all(Pf) - (@d_xa(qDx)*_dx + @d_ya(qDy)*_dy + @d_za(qDz)*_dz)*_β_dτ_D
    return nothing
end

"""
    update_temperature_flux(qTx,qTy,qTz,λ_ρCp,T,_θ_dτ_T,_dx,_dy,_dz)

Calculates temperature flux based on the gradient of temperature and material propoerties.
"""
@parallel function update_temperature_flux(qTx,qTy,qTz,λ_ρCp,T,_θ_dτ_T,_dx,_dy,_dz)
    @all(qTx) = @all(qTx) - (@all(qTx) + λ_ρCp.*(@d_xi(T)*_dx))*_θ_dτ_T
    @all(qTy) = @all(qTy) - (@all(qTy) + λ_ρCp.*(@d_yi(T)*_dy))*_θ_dτ_T
    @all(qTz) = @all(qTz) - (@all(qTz) + λ_ρCp.*(@d_zi(T)*_dz))*_θ_dτ_T
    return
end

"""
    update_dTdt!(dTdt,T,T_old,dt,max,qDx,_dx,min,qDy,_dy,qDz,_dz,ϕ)

Calculates the advection of the temperature field using an upwind scheme.
"""
@parallel_indices (i,j,k) function update_dTdt!(dTdt,T,T_old,dt,max,qDx,_dx,min,qDy,_dy,qDz,_dz,ϕ)
    nx,ny,nz=size(T)
    if (2<=i<=nx-1 && 2<=j<=ny-1 && 2<=k<=nz-1)
        dTdt[i-1,j-1,k-1] = (T[i,j,k] - T_old[i,j,k])/dt +
            (
            max(qDx[i  ,j  ,k  ],0.0)*(T[i  ,j  ,k  ] - T[i-1,j  ,k  ])*_dx +
            min(qDx[i+1,j  ,k  ],0.0)*(T[i+1,j  ,k  ] - T[i  ,j  ,k  ])*_dx +
            max(qDy[i  ,j  ,k  ],0.0)*(T[i  ,j  ,k  ] - T[i  ,j-1,k  ])*_dy +
            min(qDy[i  ,j+1,k  ],0.0)*(T[i  ,j+1,k  ] - T[i  ,j  ,k  ])*_dy +
            max(qDz[i  ,j  ,k  ],0.0)*(T[i  ,j  ,k  ] - T[i  ,j  ,k-1])*_dz +
            min(qDz[i  ,j  ,k+1],0.0)*(T[i  ,j  ,k+1] - T[i  ,j  ,k  ])*_dz
            )/ϕ
    end
    return
end

"""
    update_T!(T,qTx,qTy,qTz,dTdt,_dx,_dy,_dz,_1_dt_β_dτ_T)

Temperature update calculated from the divergence of the temperature flux and the advected temperature field
"""
@parallel function temperature_update!(T,dTdt,qTx,qTy,qTz,_dx,_dy,_dz,β_dτ_T,dt)
    @inn(T) = @inn(T)- (@all(dTdt) + @d_xa(qTx) * _dx + @d_ya(qTy) * _dy + @d_za(qTz) * _dz) / (1.0 / dt + β_dτ_T)
    return
end

"""
    function bc_x!(A)

Von Neumann boundary condition in x direction
"""
@parallel_indices (iy,iz) function bc_x!(A)
    A[1  ,iy,iz] = A[2    ,iy,iz]
    A[end,iy,iz] = A[end-1,iy,iz]
    return
end

"""
    function bc_y!(A)
        
Von Neumann boundary condition in y direction
"""
@parallel_indices (ix,iz) function bc_y!(A)
    A[ix ,1 ,iz] = A[ix,    2,iz]
    A[ix,end,iz] = A[ix,end-1,iz]
    return
end

"""
    porous_convection_3D(;nz=63,do_visu=false)
        
Porous convection solver using the pseudo-transient method and time evolution
"""
@views function porous_convection_3D(;nz = 127,nt= 2000,do_vis=false,save_arr=true)
    # physics
    lx,ly,lz    = 40.0,20.0,20.0
    k_ηf        = 1.0
    αρgx        = 0.0
    αρgy        = 0.0
    αρgz        = 1.0
    αρg         = sqrt(αρgx^2+αρgy^2+αρgz^2)
    ΔT          = 200.0
    ϕ           = 0.1
    Ra          = 1000
    λ_ρCp       = 1/Ra*(αρg*k_ηf*ΔT*lz/ϕ) # Ra = αρg*k_ηf*ΔT*lz/λ_ρCp/ϕ
    # numerics
    ny          = nz
    nx          = Int(ceil(255/127*nz))
    cfl         = 1.0/sqrt(3.1)
    re_D        = 4π
    maxiter     = 10max(nx,ny)
    ϵtol        = 1e-6
    nvis        = 50
    ncheck      = ceil(2*max(nx,ny,nz))
    # preprocessing
    dx,dy,dz    = lx/nx,ly/ny,lz/nz
    _dx,_dy,_dz = 1.0/dx,1.0/dy,1.0/dz
    xn,yn,zn    = LinRange(-lx/2,lx/2,nx+1),LinRange(-ly/2,ly/2,ny+1),LinRange(-lz,0,nz+1)
    xc,yc,zc    = av1(xn),av1(yn),av1(zn)
    θ_dτ_D      = max(lx,ly,lz)/re_D/cfl/min(dx,dy,dz)
    _θ_dτ_D     = 1.0/(1.0 + θ_dτ_D)
    _β_dτ_D     = 1.0/((re_D*k_ηf)/(cfl*min(dx,dy,dz)*max(lx,ly,lz)))
    # init
    Pf          = @zeros(nx  ,ny  ,nz  )
    r_Pf        = zeros(nx  ,ny  ,nz  )
    qDx         = @zeros(nx+1,ny  ,nz  )
    qDy         = @zeros(nx  ,ny+1,nz  )
    qDz         = @zeros(nx  ,ny  ,nz+1)   
    # initial conditions  
    T           = Data.Array([ΔT*exp(-xc[ix]^2 -yc[iy]^2 -(zc[iz]+lz/2)^2) for ix=1:nx,iy=1:ny,iz=1:nz])
    T[:,:,1]   .= ΔT/2 ; T[:,:,end] .= -ΔT/2
    T_old       = copy(T)
    dTdt        = @zeros(nx-2,ny-2,nz-2)
    r_T         = zeros(nx-2,ny-2,nz-2)
    qTx         = @zeros(nx-1,ny-2,nz-2)
    qTy         = @zeros(nx-2,ny-1,nz-2)
    qTz         = @zeros(nx-2,ny-2,nz-1)
    # visualisation dir
    if do_vis
        ENV["GKSwstype"]="nul"
        if isdir("viz3D_out")==false mkdir("viz3D_out") end
        loadpath = "viz3D_out/"; anim = Animation(loadpath,String[])
        println("Animation directory: $(anim.dir)")
        iframe = 0
    end
    println("nx,ny,nz,nt = $(nx),$(ny),$(nz),$(nt)")
    # action
    for it = 1:nt
        T_old .= T
        # time step
        dt = if it == 1
            0.1*min(dx,dy,dz)/(αρg*ΔT*k_ηf)
        else
            min(5.0*min(dx,dy,dz)/(αρg*ΔT*k_ηf),ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)), dz/maximum(abs.(qDz)))/3.1)
        end
        re_T    = π + sqrt(π^2 + ly^2/λ_ρCp/dt)
        θ_dτ_T  = max(lx,ly,lz)/re_T/cfl/min(dx,dy,dz)
        _θ_dτ_T = 1.0/(1.0 + θ_dτ_T)
        β_dτ_T  = (re_T*λ_ρCp)/(cfl*min(dx,dy,dz)*max(lx,ly,lz))
        # iteration loop
        iter = 1; err_D = 2ϵtol; err_T = 2ϵtol
        while max(err_D,err_T) >= ϵtol && iter <= maxiter
            # hydro
            @parallel compute_flux!(qDx,qDy,qDz,Pf,T,k_ηf,_dx,_dy,_dz,αρgx,αρgy,αρgz,_θ_dτ_D)
            @parallel update_Pf!(Pf,qDx,qDy,qDz,_dx,_dy,_dz,_β_dτ_D)
            # thermo
            @parallel update_temperature_flux(qTx,qTy,qTz,λ_ρCp,T,_θ_dτ_T,_dx,_dy,_dz)
            @parallel (2:(size(T,1)-1), 2:(size(T,2)-1),2:(size(T,3)-1)) update_dTdt!(dTdt,T, T_old, dt, max, qDx, _dx, min, qDy, _dy,qDz, _dz, ϕ)
            @parallel temperature_update!(T,dTdt,qTx,qTy,qTz,_dx,_dy,_dz,β_dτ_T,dt)
            # Periodic boundary condition
            @parallel (1:size(T,2),1:size(T,3)) bc_x!(T)
            @parallel (1:size(T,1),1:size(T,3)) bc_y!(T)
            if iter % ncheck == 0
                r_Pf  .= diff(Array(qDx),dims=1)./dx .+ diff(Array(qDy),dims=2)./dy .+ diff(Array(qDz),dims=3)./dz
                r_T   .= Array(dTdt) .+ diff(Array(qTx),dims=1)./dx .+ diff(Array(qTy),dims=2)./dy .+ diff(Array(qTz),dims=3)./dz
                err_D  = maximum(abs.(r_Pf))
                err_T  = maximum(abs.(r_T))
                @printf("  iter/nx=%.1f, err_D=%1.3e, err_T=%1.3e\n",iter/nx,err_D,err_T) #if itter/nx is 10 it hit 10 itterations
            end
            iter += 1
        end
        @printf("it = %d, iter/nx=%.1f, err_D=%1.3e, err_T=%1.3e\n",it,iter/nx,err_D,err_T)
        # visualisation
        if do_vis && (it % nvis == 0)
            p1=heatmap(xc,zc,Array(T)[:,ceil(Int,ny/2),:]';xlims=(xc[1],xc[end]),ylims=(zc[1],zc[end]),aspect_ratio=1,c=:turbo)
            png(p1,@sprintf("viz3D_out/%04d.png",iframe+=1))
        end
    end
    if save_arr
        save_array("out_T-$(nx)-$(ny)-$(nz)-$(nt)",convert.(Float32,Array(T)))
    end
    return T
end

end#module
