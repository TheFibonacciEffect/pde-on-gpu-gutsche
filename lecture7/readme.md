# Lecture 7

## Diffusion 2D code
Computes the 2D diffusion using an xpu implementation. Here is the output:
![diff2D](figs/diffusion_2D_xpu.gif)

## Diffusion 3D code
Same thing in 3D

## Porous convection 2D code
Computes the porous convection in 2D using an xpu implementation. Here is the resulting temperature field ans flux.

![pconvection2D](figs/out.gif)

## Porous convection 3D code
Computes the porous convection in 3D using an xpu implementation. 
Here is the cross section:
![cross sectio 3D](figs/cross_section_3D_code/0006.png)

We unfortunatly only ran the 3D code with the testing paramerters (and not the final parameters). But now we corrected the mistake, here is the final output after 2000 timsteps (instead of only 100 we used for testing).
![out 3D](figs/T_3D.png)
