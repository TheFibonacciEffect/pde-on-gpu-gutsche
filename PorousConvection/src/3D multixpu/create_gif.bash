interval=20
outpath="../../docs"
inpath="../../docs/visualisation_3D"
# cd $inpath
convert -delay 5 -loop 0 $inpath/*.png "$outpath/out.gif"
# convert -delay 10 -loop 0 *.png out.gif