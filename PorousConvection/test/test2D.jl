using Test

push!(LOAD_PATH, "../src")

function load_array(Aname,A)
    fname = string(Aname,".bin")
    fid=open(fname,"r"); read!(fid,A); close(fid)
end

include("../src/PorousConvection_2D_xpu_daint.jl")

ny=50;nx=50;nt=10
T = PorousConvection_2D_xpu.porous_convection_2D(ny,nx,nt;do_vis=false)

#test boundary conditions
@testset "Test BC" begin
    @test all(T[1,:] .== T[2,:])
    @test all(T[end,:] .== T[end-1,:])
end

# compare arrays
T_ref = zeros(Float32,nx,ny)
load_array("../bins/out_T_2D-$(nx)-$(ny)-$(nt)",T_ref)
@testset "Reference Test" begin
    @test T_ref ≈ T
end

