using Test

push!(LOAD_PATH, "../src")

function load_array(Aname,A)
    fname = string(Aname,".bin")
    fid=open(fname,"r"); read!(fid,A); close(fid)
end

include("../src/PorousConvection_3D_xpu.jl")
# _, nx,ny,nz,_ = split(aname,"-")
nz=25
ny          = nz
nx          = Int(ceil(255/127*nz))
nt = 10
T = PorousConvection_3D_xpu.porous_convection_3D(;nz=nz,nt=nt,do_vis=false,save_arr=false)

# compare arrays
T_ref = zeros(Float32,nx,ny,nz)
load_array("../bins/out_T-$(nx)-$(ny)-$(nz)-$(nt)",T_ref)
@test T_ref ≈ T

