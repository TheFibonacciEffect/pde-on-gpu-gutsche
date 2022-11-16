if [ ! -d "./PorousConvection_2D_xpu_daint_out" ]; then
  mkdir -p "./PorousConvection_2D_xpu_daint_out";
fi
scp -r class202@daint:/users/class202/github-repo/pde-on-gpu-gutsche/lecture7/scripts/PorousConvection_2D/PorousConvection_2D_xpu_daint_out "."