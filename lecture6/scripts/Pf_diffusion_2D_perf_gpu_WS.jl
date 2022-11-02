using Plots,Plots.Measures,Printf, CUDA, Test

# initilaize plotting
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)
ENV["GKSwstype"]="nul"
if isdir("viz_out")==false mkdir("viz_out") end
loadpath = "./viz_out/"; anim = Animation(loadpath,String[])
println("Animation directory: $(anim.dir)")


macro d_xa(A)  esc(:( $A[ix+1,iy]-$A[ix,iy] )) end
macro d_ya(A)  esc(:( $A[ix,iy+1]-$A[ix,iy] )) end

# triad Benchmark kernel
@inbounds function memcopy_triad_KP!(A, B, C, s)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    A[ix,iy] = B[ix,iy] + s*C[ix,iy]
    return nothing
end

# compute fluxes
function compute_flux!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ)
    nx,ny=size(Pf)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if (ix<nx && iy<ny +1)
        qDx[ix+1,iy] -= (qDx[ix+1,iy] + k_ηf_dx*@d_xa(Pf))*_1_θ_dτ
    end
    if (ix< nx + 1 && iy<ny)
        qDy[ix,iy+1] -= (qDy[ix,iy+1] + k_ηf_dy*@d_ya(Pf))*_1_θ_dτ
    end
    return nothing
end

# pressure update
function update_Pf!(Pf,qDx,qDy,_dx_β_dτ,_dy_β_dτ)
    nx,ny=size(Pf)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    if (ix<=nx && iy<=ny)
        Pf[ix,iy]  -= @d_xa(qDx)*_dx_β_dτ + @d_ya(qDy)*_dy_β_dτ
    end
    return nothing
end

function compute!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ,_dx_β_dτ,_dy_β_dτ, threads, blocks)
    CUDA.@sync @cuda blocks=blocks threads=threads compute_flux!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ)
    CUDA.@sync @cuda blocks=blocks threads=threads  update_Pf!(Pf,qDx,qDy,_dx_β_dτ,_dy_β_dτ)
    return nothing
end

function Pf_diffusion_2D(nx,ny;do_check=false)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    ϵtol    = 1e-8
    maxiter = 50
    ncheck  = ceil(Int,0.25max(nx,ny))
    ntest   = 50
    cfl     = 1.0/sqrt(2.1)
    re      = 2π
    threads = (32,8)
    blocks  = (nx÷threads[1], ny÷threads[2])
    # derived numerics
    dx,dy   = lx/nx,ly/ny
    xc,yc   = LinRange(dx/2,lx-dx/2,nx),LinRange(dy/2,ly-dy/2,ny)
    θ_dτ    = max(lx,ly)/re/cfl/min(dx,dy)
    β_dτ    = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    _1_θ_dτ = 1.0/(1.0 + θ_dτ)
    _β_dτ   = 1.0/(β_dτ)
    _dx_β_dτ = 1.0/dx/_β_dτ
    _dy_β_dτ = 1.0/dy/_β_dτ
    k_ηf_dx,k_ηf_dy = k_ηf/dx,k_ηf/dy
    # array initialisation
    # gpu arrays
    Pf      = CuArray(@. exp(-(xc-lx/2)^2 -(yc'-ly/2)^2))
    qDx,qDy = CUDA.zeros(Float64, nx+1,ny),CUDA.zeros(Float64, nx,ny+1)
    r_Pf    = CUDA.zeros(nx,ny)
    # iteration loop
    iter = 1; err_Pf = 2ϵtol
    t_tic = 0.0; niter = 0
    while err_Pf >= ϵtol && iter <= maxiter
        if (iter==11) t_tic = Base.time(); niter = 0 end
        compute!(qDx,qDy,Pf,k_ηf_dx,k_ηf_dy,_1_θ_dτ,_dx_β_dτ,_dy_β_dτ,threads,blocks)
        if do_check && (iter%ncheck == 0)
            Pf_cpu = Array(Pf)
            r_Pf  .= diff(qDx,dims=1).*_dx .+ diff(qDy,dims=2).*_dy # leave r_Pf on GPU
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            display(heatmap(xc,yc,Pf_cpu';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo))
        end
        iter += 1; niter += 1
    end

    t_toc = Base.time() - t_tic
    A_eff = (3*2)/1e9*nx*ny*sizeof(Float64)  # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                      # Execution time per iteration [s]
    T_eff = A_eff/t_it                       # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec, T_eff = %1.3f GB/s (niter = %d)\n", t_toc, round(T_eff, sigdigits=3), niter)
    return T_eff
end

function Triad!()
    # numerics
    nx,ny   = 127,127
    threads = (32,8)
    blocks  = (nx÷threads[1], ny÷threads[2])
    s       = rand()
    # array initialisation
    A       = CUDA.rand(Float64, nx, ny)
    B       = CUDA.rand(Float64, nx, ny)
    C       = CUDA.rand(Float64, nx, ny)
    
     # find T_peak
     t_it   = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy_triad_KP!($A, $B, $C, $s); synchronize() end
     T_peak = 3*1/1e9*nx*ny*sizeof(Float64)/t_it
     return T_peak
end

# get data for weak scaling
function main()
    ni      = []
    T_eff   = []
    nx = ny = 32 .* 2 .^ (0:8) .- 1
    T_peak  = Triad!()

    for i=1:size(nx)[1]
        push!(ni, nx[i]*ny[i])
        push!(T_eff,Pf_diffusion_2D(nx[i],ny[i];do_check=false))
    end

    #plot the results
    plt = plot(ni,T_eff,[ni[1],ni[end]],T_peak,
        title = "Effective memory throughput Tesla P100",
        xlabel="nx*ny",
        ylabel="T_eff [GB/s]",
        xaxis=:log,
        yaxis=:log,
        marker = 2,
        markershape=:circle,
        markersize=10
    )
    return
end


main()