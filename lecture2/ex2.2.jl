using Plots

#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

function finite_step!(C,q,dt,dx,ξ,dc,Ceq)
    # reaction
    C .+= -(C.-Ceq) / ξ *dt
    # diffusion
    q          .= -dc.*diff(C)./dx
    C[2:end-1] .+= -  dt.*diff(q)./dx
end

function reaction_diff(ic::Function,ttot,lx,Da,xi,Ceq,animation_length_seconds)
    # derived physics
    dc   = lx^2/Da/xi   # diffusion coefficient
    # Nummerics & derived nummerics
    nx = 200
    dx = lx/nx
    dt   = dx^2/dc/2

    nt = Int(round(ttot/dt))
    nvis = Int(round(nt/15/animation_length_seconds))
    x = range(0,lx,nx)
    C0 = ic.(x); C = copy(C0)
    q = zeros(nx-1)

    p = Progress(nt, showspeed=true)
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
    xi   = 10.0  # reaction rate
    ttot = 20.0  # total simulation time
    C_eq = 0.4
    gauss(x) = exp(-(x-lx/4)^2)
    Da = 0.1
    Cfinal,C0,x,anim = reaction_diff(gauss,ttot,lx,Da,xi,C_eq,3)
    Da = 1000
    Cfinal2,C02,x2,anim2 = reaction_diff(gauss,ttot,lx,Da,xi,C_eq,3)
    anim,Cfinal,C0,x , Cfinal2,C02,x2,anim2
end

# unfortunatly it is still unstable, even so I double checked, if I used the correct dt.

A = main();
anim,Cfinal,C0,x = A[1:4]
a1 = gif(anim,fps=15)
m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(round(x[n],digits=2)))",xlabel="distance",ylabel="concentration")
p1 = plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")
title!(p1,"Da = 0.1")

Cfinal,C0,x,anim = A[5:8]
a2 = gif(anim,fps=15)
m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(round(x[n],digits=2))",xlabel="distance",ylabel="concentration")
p2 = plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")
title!(p2,"Da = 1000")

p = plot(p1,p2 ,layout=(2,1))
savefig(p,"week_2/figs/ex2.2.png")