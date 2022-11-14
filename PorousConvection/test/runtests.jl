# push!(LOAD_PATH, "../src")
push!(LOAD_PATH, "../src")

using PorousConvection
using Test
function runtests()
    exename = joinpath(Sys.BINDIR, Base.julia_exename())
    testdir = pwd()
    print(pwd())
    printstyled("Testing PorousConvection.jl\n"; bold=true, color=:white)
    try
        run(`$exename -O3 --startup-file=no --check-bounds=no $(joinpath(testdir, "test2D.jl"))`)
        # run(`$exename -O3 --startup-file=no --check-bounds=no $(joinpath(testdir, "test3D.jl"))`)
    catch e
        printstyled("Error in tests: $(e)\n"; bold=true, color=:red)
        stacktrace()
        printstyled("Test failed\n"; bold=true, color=:red)
        return 1 # failure
    end
    return 0 # success
end

exit(runtests())