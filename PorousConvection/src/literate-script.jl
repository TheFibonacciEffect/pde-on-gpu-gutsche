using Literate

# directory where the markdown files are put
md_dir = joinpath(@__DIR__, "md")
if md_dir==false mkdir(md_dir) end
Literate.markdown("bin_io_script.jl", md_dir, execute=true, documenter=false, credit=false)