using Plots,Plots.Measures,Printf
using BenchmarkTools
using ProgressMeter

function compute_ap!(C2,C,A)
    C2 .= C + A 
end

function compute_kp!(C2,C,A)
    nx,ny = size(C2)
    Threads.@threads for ix ∈ 1:nx
        for iy ∈ 1:ny
            C2[ix,iy] = C[ix,iy] + A[ix,iy]
        end
    end
end

function memcopy(f,n;do_check=false,bench=:loop)
    # Numerics
    nt      = Int(1e2)
    # array initialisation
    C       = rand(Float64, n, n)
    C2      = copy(C)
    A       = copy(C)
    # performance evaluation
    t_tic = 0.0
    # iteration loop
    iter = 1;niter = 0
    t_toc = 0
    if bench==:loop
        for iter=1:nt
            # allow 11 iterations for warmup
            if iter == 11 t_tic = Base.time(); niter = 0; end
            f(C2,C,A)
            if do_check && iter%ncheck == 0
                r_Pf  .= diff(qDx,dims=1).*_dx .+ diff(qDy,dims=2).*_dy
                err_Pf = maximum(abs.(r_Pf))
                @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/n,err_Pf)
                display(heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
            end
            iter += 1; niter += 1
        end
        t_toc = Base.time() - t_tic
    elseif bench == :btool
        @infiltrate
        (f == compute_kp!) && (t_toc = @belapsed compute_kp!($C2,$C,$A))
        (f == compute_ap!) && (t_toc = @belapsed compute_ap!($C2,$C,$A))
        niter = 1
        @infiltrate
    end
    # read 2, write 1
    A_eff = 3*n*n*sizeof(eltype(C2))/1e9   # Effective main memory access per iteration [GB]
    t_it  = t_toc/niter                        # Execution time per iteration [s]
    T_eff = A_eff/t_it                         # Effective memory throughput [GB/s]
    @printf("Time = %1.3f sec \n", t_toc)
    @printf("Time iter = %1.3f sec \n", t_it)
    @printf("T_eff = %1.3f GB/sec \n", T_eff)
    @printf("niter = %i \n", niter)
    @infiltrate
    return T_eff
end

memcopy(compute_kp!,512;do_check=false,bench = :btool)
memcopy(compute_ap!,512;do_check=false,bench = :btool)
memcopy(compute_kp!,512;do_check=false,bench = :loop)
memcopy(compute_ap!,512;do_check=false,bench = :loop)

#task 2
T_kp_btool = []
T_ap_btool = []
T_kp_loop  = []
T_ap_loop  = []
N = 16 * 2 .^ (1:8)
p = Progress(sum(N.^2),1.0)
for n in N
    push!(T_kp_btool,memcopy(compute_kp!,n;do_check=false,bench = :btool))
    push!(T_ap_btool,memcopy(compute_ap!,n;do_check=false,bench = :btool))
    push!(T_kp_loop ,memcopy(compute_kp!,n;do_check=false,bench = :loop ))
    push!(T_ap_loop ,memcopy(compute_ap!,n;do_check=false,bench = :loop ))
    next!(p;step=n^2)
end

@show findmax(T_kp_btool)
@show findmax(T_ap_btool)
@show findmax(T_kp_loop )
@show findmax(T_ap_loop)

default( marker = :diamond,seriestype=:line)
plot( N,T_kp_btool,xaxis=:log, label="T_kp_btool")
plot!(N,T_kp_loop ,xaxis=:log, label="T_kp_loop" )
plot!(N,T_ap_btool,xaxis=:log, label="T_ap_btool")
plot!(N,T_ap_loop ,xaxis=:log, label="T_ap_loop" )
ylabel!("\$T_{eff}\$ in [GB/s]")
xlabel!("n")
title!(raw"memory thourghput of $n \times n$ array")
savefig("/figs/memcopy.png")
