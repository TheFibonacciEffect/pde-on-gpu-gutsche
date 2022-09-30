using Plots

#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

function finite_step!(C,q,dt,dx,ξ,dc,Ceq)
    # I am a little confused about what boundary conditions to use
    # reaction
    C -= (C.-Ceq) / ξ *dt
    # diffusion
    q          .= -dc.*diff(C)./dx
    C[2:end-1] .-=   dt.*diff(q)./dx
end

function reaction_diff(ic::Function,ttot,lx,dc,xi,Ceq)
    # Nummerics
    nx = 200
    dx = lx/nx
    @show dt   = dx^2/dc/2
    nt = Int(round(ttot/dt))
    nvis = 3
    x = range(0,lx,nx)
    C0 = ic.(x); C = copy(C0)
    q = zeros(nx-1)

    p = Progress(nt, 0.2)
    anim = @animate for i in 1:nt
        finite_step!(C,q,dt,dx,xi,dx,Ceq)
        # visualisation
        plot(x,C,label="concentration at t=$(round(dt*i,digits=1))",xlabel="distance",ylabel="concentration")
        plot!(x,C0,label="inital concentration")
        next!(p)
    end every nvis
    return C,C0,x,anim
end

function main()
    # Physics
    lx   = 20.0  # domain length
    dc   = 0.1   # diffusion coefficient
    xi   = 10.0  # reaction rate
    ttot = 20.0  # total simulation time
    C_eq = 0.4
    gauss(x) = exp(-(x-lx/4)^2)
    Cfinal,C0,x,anim = reaction_diff(gauss,ttot,lx,dc,xi,C_eq)
    anim,Cfinal,C0,x
end

anim,Cfinal,C0,x = main()
gif(anim,"week_2/tmp/reaction_diff.gif",fps=15)


m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(x[n])",xlabel="distance",ylabel="concentration")
plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")