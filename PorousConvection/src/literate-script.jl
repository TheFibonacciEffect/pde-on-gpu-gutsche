using Literate

# directory where the markdown files are put
md_dir = "../docs/md"
if md_dir==false mkdir(md_dir) end
Literate.markdown("bin_io_script.jl", md_dir, execute=true, documenter=false, credit=false,mdstrings=true,comments=false)