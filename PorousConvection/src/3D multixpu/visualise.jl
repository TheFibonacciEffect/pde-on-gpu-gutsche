using GLMakie

function load_array(Aname,A)
    fname = string(Aname)
    fid=open(fname,"r"); read!(fid,A); close(fid)
end

function visualise(aname)
    lx,ly,lz = 40.0,20.0,20.0
    _, nx,ny,nz,nt,_ = split(aname,"-")
    nx = parse(Int,nx); ny = parse(Int,ny); nz = parse(Int,nz); nt= parse(Int,nt)
    T  = zeros(Float32,nx,ny,nz)
    load_array(aname,T)
    @assert maximum(T) != 0.0
    xc,yc,zc = LinRange(0,lx,nx),LinRange(0,ly,ny),LinRange(0,lz,nz)
    fig      = Figure(resolution=(1600,1000),fontsize=24)
    ax       = Axis3(fig[1,1];aspect=(1,1,0.5),title="Temperature after $nt",xlabel="lx",ylabel="ly",zlabel="lz")
    surf_T   = contour!(ax,xc,yc,zc,T;alpha=0.05,colormap=:turbo)
    return fig
end

function main()
    indir = "../../docs/l9ex1t2_out/"
    outdir = "../../docs/visualisation_3D"
    if !isdir(outdir) mkdir(outdir) end
    for (it,file) in enumerate(readdir(indir))
        if occursin(".bin",file)
            println("working on ",file, "...")
            fig = visualise(indir*file)
            save(outdir*"/$(lpad(it,4,"0")).png",fig)
        end
    end
end

main()
println("Done")