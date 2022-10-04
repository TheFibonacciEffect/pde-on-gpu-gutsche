using Plots

#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

function finite_step!(C,q,dt,dx,vx,dc)
    # advection
    C[2:end  ] .-= dt.*max(vx,0.0).*diff(C)./dx # always take the derivative in the direction of movement
    C[1:end-1] .-= dt.*min(vx,0.0).*diff(C)./dx
    # diffusion
    q           .=  -dc.*diff(C)./dx
    C[2:end-1] .-=   dt.*diff(q)./dx
end

function avec_diff(ic::Function,lx,Pe,vx,ttot)
    # derived Physics
    dc = lx*vx/Pe
    # Nummerics & derived Nummerics
    nx = 200
    dx = lx/nx
    dta  = dx/abs(vx)
    dtd  = dx^2/dc/2
    dt   = min(dtd, dta)
    @show nt = Int(round(ttot/dt))
    @assert nt < 10_000
    nvis = 3
    x = range(0,lx,nx)
    C0 = ic.(x); C = copy(C0)
    q = zeros(nx-1)
    p = Progress(nt, 0.2)
    anim = @animate for i in 1:nt
        finite_step!(C,q,dt,dx,vx,dc)
        # visualisation
        plot(x,C,label="concentration at t=$(round(dt*i,digits=1))",xlabel="distance",ylabel="concentration")
        plot!(x,C0,label="inital concentration")
        next!(p)
        if i*dt > ttot/2
            vx = -1.
        end
    end every nvis
    return C,C0,x,anim
end

function main()
    # Physics
    lx   = 20.0  # domain length
    vx   = 1.0   # advection velocity
    ttot = 20.0  # total simulation time
    gauss(x) = exp(-(x-lx/4)^2)
    Pe = 10
    Cfinal,C0,x,anim = avec_diff(gauss,lx,Pe,vx,ttot)
    Pe = 1000
    Cfinal2,C02,x2,anim2 = avec_diff(gauss,lx,Pe,vx,ttot)
    anim,Cfinal,C0,x , Cfinal2,C02,x2,anim2
end

A = main();
anim,Cfinal,C0,x = A[1:4]
a1 = gif(anim,fps=15)
m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(round(x[n],digits=2)))",xlabel="distance",ylabel="concentration")
p1 = plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")
title!(p1,"Pe = 10")

Cfinal,C0,x,anim = A[5:8]
a2 = gif(anim,fps=15)
m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(round(x[n],digits=2))",xlabel="distance",ylabel="concentration")
p2 = plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")
title!(p2,"Pe = 1000")

p = plot(p1,p2 ,layout=(2,1))
savefig(p,"week_2/figs/ex1.2.png")