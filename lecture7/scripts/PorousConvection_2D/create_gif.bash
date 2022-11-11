interval=20
figs_path="../../figs"
imgs_path="$figs_path/PorousConvection_2D_xpu_daint_out"
convert -delay 10 -loop 0 $imgs_path/*.png $figs_path/out.gif