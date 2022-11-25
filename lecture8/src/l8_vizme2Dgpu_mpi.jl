# Visualisation script for the 2D MPI solver
using Plots, MAT

nprocs = (2, 2) # nprocs (x, y) dim

@views function vizme2D_mpi(dir,nprocs,timestep)
    C  = []
    ip = 1
    for ipx = 1:nprocs[1]
        for ipy = 1:nprocs[2]
            file = matopen(string(dir,"$(ip-1)_$timestep.mat")); C_loc = read(file, "C"); close(file)
            nx_i, ny_i = size(C_loc,1)-2, size(C_loc,2)-2
            ix1, iy1   = 1+(ipx-1)*nx_i, 1+(ipy-1)*ny_i
            if (ip==1)  C = zeros(nprocs[1]*nx_i, nprocs[2]*ny_i)  end
            C[ix1:ix1+nx_i-1,iy1:iy1+ny_i-1] .= C_loc[2:end-1,2:end-1]
            ip += 1
        end
    end
    fontsize = 12
    opts = (aspect_ratio=1, yaxis=font(fontsize, "Courier"), xaxis=font(fontsize, "Courier"),
        ticks=nothing, framestyle=:box, titlefontsize=fontsize, titlefont="Courier", 
        xlabel="Lx", ylabel="Ly", xlims=(1, size(C,1)), ylims=(1, size(C,2)) )
    display(heatmap(C'; c=:turbo, title="diffusion 2D MPI", opts...))
    return
end

if (length(ARGS) < 1)
    # println("Usage: julia l8_vizme2D_mpi.jl <dir> <timestep>")
    push!(ARGS,"../docs/l8ex1t3_out/mpi2Dgpu_out_")
end
    
an = @animate for i in 5:5:100
    vizme2D_mpi(ARGS[1],nprocs,lpad(i,3,'0'))
end

gif(an, "../docs/$(split(ARGS[1],"_")[1]).gif", fps = 5)