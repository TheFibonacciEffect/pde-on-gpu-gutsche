using Plots

#imports a nice progress bar
using ProgressMeter

# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

function finite_step!(C,q,dt,dx,vx,dc)
    # I am a little confused about what boundary conditions to use
    # advection
    C[2:end  ] .-= dt.*max(vx,0.0).*diff(C)./dx # always take the derivative in the direction of movement
    C[1:end-1] .-= dt.*min(vx,0.0).*diff(C)./dx
    # diffusion
    q          .= -dc.*diff(C)./dx
    C[2:end-1] .-=   dt.*diff(q)./dx
end

function avec_diff(ic::Function,lx,dc,vx,ttot)
    # Nummerics
    nx = 200
    dx = lx/nx
    dta  = dx/abs(vx)
    dtd  = dx^2/dc/2
    dt   = min(dtd, dta)
    nt = Int(round(ttot/dt))
    nvis = 3
    x = range(0,lx,nx)
    C0 = ic.(x); C = copy(C0)
    q = zeros(nx-1)

    p = Progress(nt, 0.2)
    anim = @animate for i in 1:nt
        finite_step!(C,q,dt,dx,vx,dx)
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
    vx   = 1.0   # advection velocity
    ttot = 20.0  # total simulation time
    gauss(x) = exp(-(x-lx/4)^2)
    Cfinal,C0,x,anim = avec_diff(gauss,lx,dc,vx,ttot)
    anim
end

anim = main()
gif(anim,"week_2/tmp/advection_diff.gif",fps=15)
m,n = findmax(Cfinal)
plot(x,Cfinal,label="final concentration  (max:$(round(m,digits=2)) at $(x[n])",xlabel="distance",ylabel="concentration")
plot!(x,C0,label="inital concentration",xlabel="distance",ylabel="concentration")
