using Plots
#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

abstract type PDE_system end

Base.@kwdef mutable struct Acoustic <: PDE_system
    nx = 128
    ny = 129
end    



function finite_step!(s::PDE_system)
    # inviscid Burgers
    s.C[2:end] .+= -  max(s.v,0.) .* s.dt.*diff(s.C.^s.n)./s.dx
    s.C[1:end-1] .+= -  min(s.v,0) .* s.dt.*diff(s.C.^s.n)./s.dx
end

"modefies the Setup type in the input so that it contains the solution in the Setup.C field"
function solve_pde!(s::PDE_system)
p = Progress(s.nt, 0.2)
    anim = @animate for i in 1:s.nt
        finite_step!(s)
        # visualisation
        plot(s.xc,s.C,label="concentration at t=$(round(s.dt*i,digits=1))",xlabel="distance",ylabel="concentration")
        plot!(s.xc,s.C0,label="inital concentration")
        next!(p)
        if i > s.nx #change velocyity after half the time
            s.v = -1.
        end
    end every s.nvis
    return anim
end