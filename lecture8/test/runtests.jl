using Test, MAT

file = matopen("../docs/l8ex2t2/gputrue_out.mat")
gputrue = read(file, "C"); close(file)
file = matopen("../docs/l8ex2t2/gpufalse_out.mat")
gpufalse = read(file, "C"); close(file)
file = matopen("../docs/l8ex2t2/mpigpu_out.mat")
mpigpu = read(file, "C"); close(file)

#test boundary conditions
@testset "Unit test" begin
    @test all(gputrue ≈ gpufalse) 
    @test all (gputrue ≈ mpigpu)
end
