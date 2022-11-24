ENV["GKSwstype"]="nul";
using Plots

# read csv
fname = "weak_scaling.dat"
fid = open(fname,"r")
lines = readlines(fid)
close(fid)

# parse csv
nprocs = Array{Int}(undef,length(lines))
t = Array{Float64}(undef,length(lines))
for (i,line) in enumerate(lines)
    nprocs[i],t[i] = split(line,",") |> x -> parse.(Float64,x)
end
t_norm = t./t[1]
# plot
p = plot(nprocs,t_norm./nprocs,marker=:circle,legend=false)
plot!(p, xlabel="nprocs", ylabel="normlized t/nprocs", title="Weak scaling")
plot!(p, size=(400,400))
savefig(p, "weak_scaling.png")