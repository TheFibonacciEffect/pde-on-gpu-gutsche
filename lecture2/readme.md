# Week 2

## Plots
### Advenction Diffusion
![gif](figs/advection_diff.gif)
#### Question 1
The plot with the inital and final concenttation is generated when running the script.
#### Question 2
We can see that the fluid is fully diffused for small $Pe$ while for large $Pe$ is is not. This matches with the intended behaviour.
![Q2](figs/ex1.2.png)

### Reaction Diffusion
$$ \frac{∂C}{∂t} = -\frac{(C-C_\mathrm{eq})}{ξ}$$

![gif](figs/reaction_diff.gif)
#### Question 1
The plot with the inital and final concenttation is generated when running the script
#### Question 2
Unfortunatly my code becomes unstable when running with a $Da > 1000$. I was unable to determine the reason for this.
![Q2](figs/ex2.2.png)

## Nonlinear Problems
### Task 1
![T1](figs/power-law.gif)

### Task 2
![T2](figs/burgers.gif)

# Ex 4
See https://github.com/TheFibonacciEffect/pde-on-gpu-gutsche