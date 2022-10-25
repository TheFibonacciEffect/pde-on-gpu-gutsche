# TODO physics computations in compute!() function, derivatives done with macros, and multi-threading
using Plots,Plots.Measures,Printf
using BenchmarkTools
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

macro d_xa(A)  esc(:( $A[ix+1,iy]-$A[ix,iy]) ) end
macro d_ya(A)  esc(:( $A[ix,iy+1]-$A[ix,iy] )) end

@inbounds function compute_flux!(qDx, k_ηf, Pf, _dx, _1_θ_dτ, qDy, _dy)
    nx,ny=size(Pf)
    for iy=1:ny
        for ix=1:nx-1
            qDx[ix+1,iy] -= (qDx[ix+1,iy] + k_ηf*((@d_xa(Pf)*_dx)))*_1_θ_dτ
        end
    end
    for iy=1:ny-1
        for ix=1:nx
            qDy[ix,iy+1] -= (qDy[ix,iy+1] + k_ηf*((@d_ya(Pf))*_dy))*_1_θ_dτ
        end
    end
end
    
@inbounds function update_Pf!(Pf, qDx, _dx, qDy, _dy, _β_dτ)
    nx,ny=size(Pf)
    for iy=1:ny
        for ix=1:nx
            Pf[ix,iy]  -= ((qDx[ix+1,iy]-qDx[ix,iy])*_dx + (qDy[ix,iy+1]-qDy[ix,iy])*_dy)*_β_dτ
        end
    end
end

function compute_ap!(C2,C,A)
    C2 .= C + A 
end

function compute_kp!(C2,C,A)
    nx,ny = size(C2)
    for ix ∈ 1:nx
        for iy ∈ 1:ny
            C2[ix,iy] = C[ix,iy] + A[ix,iy]
        end
    end
end

function memcopy(;do_check=false,bench=:loop)
    # Numerics
    nx, ny  = 512, 512
    nt      = 2e4
    # array initialisation
    C       = rand(Float64, nx, ny)
    C2      = copy(C)
    A       = copy(C)
    # performance evaluation
    t_tic = 0.0
    # iteration loop
    iter = 1
    t_tic = 0.0; niter = 0
    t_toc = 0
    if bench==:loop
    for iter=1:nt
        if iter == 11 t_tic = Base.time(); niter = 0; end
        compute_ap!(C2,C,A)
        if do_check && iter%ncheck == 0
            r_Pf  .= diff(qDx,dims=1).*_dx .+ diff(qDy,dims=2).*_dy
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            display(heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        end
        iter += 1; niter += 1
        niter = iter-11
    end
    t_toc = Base.time() - t_tic
    elseif bench == :btool
        t_toc = @belapsed compute_ap!($C2,$C,$A)
        niter = 1
    end
    # read and write qDx, qDy and Pf --> 3 * read-write --> 3*2
    A_eff = 3*2*nx*ny*sizeof(eltype(C2))/1e9   # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                        # Execution time per iteration [s]
    T_eff = A_eff/t_it                         # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec \n", t_toc)
    @printf("Time iter = %1.3f sec \n", t_it)
    @printf("T_eff = %1.3f GB/sec \n", T_eff)
    @printf("niter = %1.3f \n", niter)
    @infiltrate
end

memcopy(;do_check=false,bench = :btool)
memcopy(;do_check=false,bench = :loop)
