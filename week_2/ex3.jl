using Revise

using Plots
#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

Base.@kwdef mutable struct Setup
    # physics
    lx   = 20.0
    dc   = 1.0
    n    = 4
    # numerics
    nx   = 200
    nvis = 50
    # derived numerics
    dx   = lx/nx
    dt   = dx^2/dc/10
    nt   = nx^2 ÷ 5
    # inital condition
    ic = x-> exp(-(x-lx/4)^2)
    # arrays
    xc   = LinRange(s.dx/2,s.lx-s.dx/2,s.nx)
    # C0 = ic.(xc)
    C0 = ic.(xc)
    C = copy(C)
    q = zeros(nx-1)
    # collections for easy use
    dgl_params=[dc,n]
end


function finite_step!(s::Setup)
    # diffusion
    s.q          .= -s.dc.*diff(s.C.^s.n)./s.dx
    s.C[2:end-1] .+= -  s.dt.*diff(q)./s.dx
end

function reaction_diff!(s::Setup)
p = Progress(s.nt, 0.2)
anim = @animate for i in 1:s.nt
    finite_step!(s)
    # visualisation
    plot(s.xc,s.C,label="concentration at t=$(round(dt*i,digits=1))",xlabel="distance",ylabel="concentration")
    plot!(s.xc,s.C0,label="inital concentration")
    next!(p)
end every s.nvis
return anim
end


# In numpy there is
# `np.set_printoptions(formatter={'float_kind':float_formatter})`
function main()
    s = Setup()
    anim = reaction_diff!(s,20,3)
    anim
end

anim = main();
gif(anim,"week_2/figs/reaction_diff.gif",fps=15)

m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(x[n])",xlabel="distance",ylabel="concentration")
plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")
