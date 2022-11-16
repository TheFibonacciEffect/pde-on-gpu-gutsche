interval=20
figs_path="../../figs"
imgs_path="./PorousConvection_2D_xpu_daint_out"
# cd $imgs_path
convert -delay 5 -loop 0 $imgs_path/*.png "$figs_path/out.gif"
# convert -delay 10 -loop 0 *.png out.gif