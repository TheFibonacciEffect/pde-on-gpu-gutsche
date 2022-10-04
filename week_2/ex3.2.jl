using Plots
#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

"Struct to initialize all variables needed for the computation"
Base.@kwdef mutable struct Setup
    # physics
    lx   = 20.0
    vx   = 1.0
    n    = 2
    # numerics
    nx   = 1000
    nvis = 15
    # derived numerics
    dx   = lx/nx
    dt   = dx/abs(vx)/2
    nt   = 2nx
    xc   = LinRange(dx/2,lx-dx/2,nx)
    # inital condition
    ic = x-> exp(-(x-lx/4)^2)
    # C0 = ic.(xc)
    C0 = ic.(xc)
    C = copy(C0)
    q = zeros(nx-1)
end


function finite_step!(s::Setup)
    # inviscid Burgers
    s.C[2:end] .+= -  s.dt.*diff(s.C.^s.n)./s.dx
end

"modefies the Setup type in the input so that it contains the solution in the Setup.C field"
function solve_pde!(s::Setup)
    p = Progress(s.nt, 0.2)
    anim = @animate for i in 1:s.nt
        finite_step!(s)
        # visualisation
        plot(s.xc,s.C,label="concentration at t=$(round(s.dt*i,digits=1))",xlabel="distance",ylabel="concentration")
        plot!(s.xc,s.C0,label="inital concentration")
        next!(p)
    end every s.nvis
    return anim
end

function main()
    s = Setup()
    anim = solve_pde!(s)
    anim,s
end

anim,s = main();
gif(anim,"week_2/figs/burgers.gif",fps=15)

m,n = findmax(s.C)
plot(s.xc,s.C,label="final concentration  (max:$(round(m,digits=2)) at $(s.xc[n])",xlabel="distance",ylabel="concentration")
plot!(s.xc,s.C0,label="inital concentration",xlabel="distance",ylabel="concentration")
