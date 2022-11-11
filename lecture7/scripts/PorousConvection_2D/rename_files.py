#!/usr/bin/env python
import sys,os,re

figs_path="../../figs"
imgs_path=f"{figs_path}/PorousConvection_2D_xpu_daint_out"

for file in os.listdir(imgs_path):
    if __file__ in file: continue
    words = re.split('-|_|\.',file)
    words[-2] = str(int(words[-2])//20).zfill(4)
    new_name = "_".join(words[:-1]) + "." + words[-1]
    print(f"{file} \t-->{new_name}")
    os.rename(f"{imgs_path}/{file}",f"{imgs_path}/{new_name}")