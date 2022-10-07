using Plots
#imports a nice progress bar
using ProgressMeter
using PropertyUtils


# supress warnings about No strict ticks found
import Logging
Logging.disable_logging(Logging.Warn)

abstract type PDE_system end


Base.@kwdef mutable struct Acoustic <: PDE_system
    # physics
    c = 1
    ρ = (c/2)^2
    σ = (c/2)^2
    # inital condition
    ic = (x,y) -> exp(-(x^2+ y^2))
    Lx = 10
    Ly = 10
    # nummerics
    nx = 128
    ny = 129
    dx = Lx/nx
    dy = Ly/ny
    dt = min(dx,dy) /c
    # derived nummerics
    x = range(-Lx/2,Lx/2,nx)
    y = range(-Ly/2,Ly/2,ny)
    Vx = zeros(nx-1,ny)
    Vy = zeros(nx,ny-1)
    P0 = ic.(x,y')
    P = copy(P0)
end

function finite_step!(s::PDE_system)
    # wave equation
    @infiltrate
    @with s begin
        s.Vx .+= dt .* diff(P,dims=1)./dx ./ρ
        s.Vy .+= dt .* diff(P,dims=2)./dx ./ρ
        s.P[2:end-1,:] .+= dt .* diff(Vx,dims=1) ./ dx ./σ
        s.P[:,2:end-1] .+= dt .* diff(Vy,dims=2) ./ dx ./σ
    end 
end

"modefies the PDE_system in the input so that it contains the solution in corrsponding field"
function solve_pde!(s::PDE_system,anim_time,ttot)
    fps = 15
    nt = Int(ceil(ttot/s.dt))
    nvis = Int(ceil(nt/(anim_time*fps)))
    p = Progress(nt, 0.2)
    anim = @animate for i in 1:nt
        finite_step!(s)
        # visualisation
        heatmap(s.x,s.y,s.P)
        @infiltrate
        next!(p)
    end every nvis
    return anim
end

# for testing
s = Acoustic()
finite_step!(s)
an = solve_pde!(s,4,20)
gif(an,fps=15)
