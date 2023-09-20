interval=20
outpath="../../docs"
inpath="../../docs/visualisation_3D_out"
# cd $inpath
convert -delay 5 -loop 0 $inpath/*.png "$outpath/porous_conv_multixpu.gif"
# convert -delay 10 -loop 0 *.png out.gif