const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences3D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end

using Plots,Plots.Measures,Printf
default(
    size=(600,500),
    framestyle=:box,
    label=false,
    grid=false,
    margin=10mm,
    lw=6,
    labelfontsize=11,
    tickfontsize=11,
    titlefontsize=11
)

@parallel function compute_flux!(qDx,qDy,qDz,Pf,k_ηf_dx,k_ηf_dy,k_ηf_dz,_1_θ_dτ)
    @inn_x(qDx) = @inn_x(qDx) - (@inn_x(qDx) + k_ηf_dx*@d_xa(Pf))*_1_θ_dτ
    @inn_y(qDy) = @inn_y(qDy) - (@inn_y(qDy) + k_ηf_dy*@d_ya(Pf))*_1_θ_dτ
    @inn_z(qDz) = @inn_z(qDz) - (@inn_z(qDz) + k_ηf_dz*@d_za(Pf))*_1_θ_dτ
    return nothing
end

@parallel function update_Pf!(Pf,qDx,qDy,qDz,_β_dτ_dx,_β_dτ_dy,_β_dτ_dz)
    @all(Pf) = @all(Pf) - (@d_xa(qDx)*_β_dτ_dx + @d_ya(qDy)*_β_dτ_dy + @d_za(qDz)*_β_dτ_dz)
    return nothing
end

function cumpute!(qDx,qDy,qDz,Pf,k_ηf_dx,k_ηf_dy,k_ηf_dz,_1_θ_dτ,_β_dτ_dx,_β_dτ_dy,_β_dτ_dz)
    @parallel compute_flux!(qDx,qDy,qDz,Pf,k_ηf_dx,k_ηf_dy,k_ηf_dz,_1_θ_dτ)
    @parallel update_Pf!(Pf,qDx,qDy,qDz,_β_dτ_dx,_β_dτ_dy,_β_dτ_dz)
    return nothing
end

function Pf_diffusion_3D(;nx=32,ny=32,nz=32,do_check=false)
    # physics
    lx,ly,lz = 20.0,20.0,20.0
    k_ηf     = 1.0
    # numerics
    ϵtol     = 1e-8
    cfl      = 1.0/sqrt(3.1)
    re       = 2π
    maxiter  = 50 #max(nx,ny)
    ncheck   = ceil(Int,0.25max(nx,ny,nz)) 
    # derived numerics
    dx       = lx/nx
    dy       = ly/ny
    dz       = ly/nz
    xc       = LinRange(dx/2,lx-dx/2,nx)
    yc       = LinRange(dy/2,ly-dy/2,ny)
    zc       = LinRange(dz/2,lz-dz/2,nz)
    θ_dτ     = max(lx,ly,lz)/re/cfl/min(dx,dy,dz)
    β_dτ     = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    _1_θ_dτ  = 1.0/(1.0 + θ_dτ)
    _β_dτ_dx = 1.0/β_dτ/dx
    _β_dτ_dy = 1.0/β_dτ/dy
    _β_dτ_dz = 1.0/β_dτ/dz
    k_ηf_dx  = k_ηf/dx
    k_ηf_dy  = k_ηf/dy
    k_ηf_dz  = k_ηf/dz
    # array initialisation
    Pf       = Data.Array([exp(-(xc[ix]-lx/2)^2 -(yc[iy]-ly/2)^2 -(zc[iz]-lz/2)^2) for ix=1:nx,iy=1:ny,iz=1:nz])
    qDx      = @zeros(nx+1,ny  , nz  )
    qDy      = @zeros(nx  ,ny+1, nz  )
    qDz      = @zeros(nx  ,ny  , nz+1)
    r_Pf     = @zeros(nx  ,ny  , nz  )
    # visualisation dir
    if do_check
        ENV["GKSwstype"]="nul"
        if isdir("viz_out")==false mkdir("viz_out") end
        loadpath = "viz_out/"; anim = Animation(loadpath,String[])
        println("Animation directory: $(anim.dir)")
        iframe = 0
    end
    # iteration loop
    iter = 1; err_Pf = 2ϵtol
    t_tic = 0.0; niter = 0
    while err_Pf >= ϵtol && iter <= maxiter
        if (iter==11) t_tic = Base.time(); niter = 0 end
        cumpute!(qDx,qDy,qDz,Pf,k_ηf_dx,k_ηf_dy,k_ηf_dz,_1_θ_dτ,_β_dτ_dx,_β_dτ_dy,_β_dτ_dz)
        if do_check && (iter%ncheck == 0)
            r_Pf  .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy .+ diff(qDz,dims=3)./dz
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            png(
                (
                    heatmap(
                        xc,yc,Array(Pf)[:,ceil(Int,ny/2),:]';
                        xlims=(xc[1],xc[end]),
                        ylims=(zc[1],zc[end]),
                        aspect_ratio=1,c=:turbo
                    )
                ),
                @sprintf("viz_out/%04d.png",iframe+=1)
            )
        end
        iter += 1; niter += 1
    end
    t_toc = Base.time() - t_tic
    A_eff = (3*2)/1e9*nx*ny*nz*sizeof(Float64)  # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                      # Execution time per iteration [s]
    T_eff = A_eff/t_it                       # Effective memory throughput [GB/s]
    @printf(
                "Time = %1.3f sec, T_eff = %1.3f GB/s (niter = %d)\n",
                t_toc,
                round(T_eff, sigdigits=3),
                niter
            )
    return
end

Pf_diffusion_3D(do_check=true)