# task 2
using BenchmarkTools

include("Pf_diffusion_2D/Pf_diffusion_2D_Teff.jl")
include("Pf_diffusion_2D/l5_Pf_Diffusion_2D_perf.jl")
include("Pf_diffusion_2D/Pf_diffusion_2D_loop_fun_parralel.jl")

Nx = 16 * 2 .^ (1:8)

Arr_T_eff      = []
Arr_T_eff_perf = []
Arr_T_eff_loop = []

# TODO use BenchmarkTools
for nx in Nx
    ny = nx
    T_eff      = Pf_diffusion_2D_Teff(nx,ny;do_check=false)
    T_eff_perf = l5_Pf_Diffusion_2D_perf(nx,ny;do_check=false)
    T_eff_loop = Pf_diffusion_2D_loop_fun_parralel(nx,ny)
    push!( Arr_T_eff     ,T_eff      )
    push!( Arr_T_eff_perf,T_eff_perf )
    push!( Arr_T_eff_loop,T_eff_loop )
end

default( marker = :diamond,seriestype=:line,xaxis=:log,yaxis=:log)
plot( Nx,Arr_T_eff      , label= "T_eff      ")
plot!(Nx,Arr_T_eff_perf , label= "T_eff_perf ")
plot!(Nx,Arr_T_eff_loop , label= "T_eff_loop ")
hline!([15.33],markeralpha = 0.,linestyle=:dash,label="best value for array programming")
hline!([5.89],markeralpha = 0.,linestyle=:dash,label="best value for loop kernel programming")
hline!([45.8],markeralpha = 0.,linestyle=:dash,label="value supplied by vendor")
ylabel!("\$T_{eff}\$ in [GB/s]")
xlabel!("n")
title!(raw"time to evalueate Pf_diffusion_2D using different methods")
savefig("./figs/Pf_diffusion_2D_Teff.png")

