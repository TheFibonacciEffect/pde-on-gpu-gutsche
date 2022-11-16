include("./PorousConvection_2D_xpu_daint.jl")

function save_array(Aname,A)
    fname = string(Aname,".bin")
    out = open(fname,"w"); write(out,A); close(out)
end

ny=50;nx=50;nt=10
T = porous_convection_2D(ny,nx,nt;do_vis=false)
save_array("bins/out_T_2D-$(nx)-$(ny)-$(nt)",convert.(Float32,Array(T)))
