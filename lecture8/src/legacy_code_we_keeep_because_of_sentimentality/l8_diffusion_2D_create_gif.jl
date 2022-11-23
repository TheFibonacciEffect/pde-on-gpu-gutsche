using Plots
using Printf

function load_array!(Aname,A)
    fname = string(Aname)
    # TODO use filesize("file.dat") to determine array dimensions
    fid=open(fname,"r"); read!(fid,A); close(fid)
    println("Loaded $(Aname) from file $(fname)")
end

process_location = Dict([
    0 => [1,1],
    1 => [1,2],
    2 => [2,1],
    3 => [2,2],
])

function sort_by_timestep(l)
    A = similar(l,Int)
    for (i,x) in enumerate(l)
        A[i] = split(x,"_")[end][1:end-4] |> x -> parse(Int,x)
    end
    # sort l by A
    A,l = transpose(sort(collect(zip(A, l)); by=first))
    @infiltrate
    return l
end

# list files in directory
function list_files(dir)
    files = readdir(dir)
    files_by_process = Array{String}(undef,2,2,length(files) ÷ 4)
    for p ∈ 0:3
        files_process = filter(x -> occursin("C_$p",x), files)
        files_process = sort(files_process)
        files_by_process[process_location[p]...,:] = files_process
    end
    return files_by_process
end

function combine_arrays(n,files_by_process,dir)
    A11 = Array{Float64}(undef,32,32)
    load_array!(string(dir,"/",files_by_process[1,1,n]),A11)
    A12 = Array{Float64}(undef,32,32)
    load_array!(string(dir,"/",files_by_process[1,2,n]),A12)
    A21 = Array{Float64}(undef,32,32)
    load_array!(string(dir,"/",files_by_process[2,1,n]),A21)
    A22 = Array{Float64}(undef,32,32)
    load_array!(string(dir,"/",files_by_process[2,2,n]),A22)
    A = [A11 A12; A21 A22]
    return A
end

# create gif
function create_gif(dir, gifname)
    files_by_process = list_files(dir)
    anim = @animate for i in 1:(size(files_by_process)[end])
        A = combine_arrays(i,files_by_process,dir)
        p = heatmap(A', aspect_ratio=:equal, c=:viridis, title=@sprintf("t=%d",i))
        plot!(p, size=(400,400))
    end
    gif(anim, string(dir,"/",gifname), fps = 5)
end

try
    create_gif("../docs/l8ex1t2_out/", "l8ex1t2.gif")
catch IOError
    create_gif("./lecture8/docs/l8ex1t2_out", "l8ex1t2.gif")
end 

