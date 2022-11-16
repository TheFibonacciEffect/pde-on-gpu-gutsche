include("./PorousConvection_3D_xpu.jl")

nz=25;nt=10
cd("../bins/")
PorousConvection_3D_xpu.porous_convection_3D(;nz=nz,nt=nt,do_vis=false,save_arr=true)
cd("../src/")