# Welcome 
to the Homework Directory of Caspar Gutsche for the course [Solving Solving partial differential equations in parallel on GPUs](https://pde-on-gpu.vaw.ethz.ch/) by Ludovic Räss,   Mauro Werder,   Samuel Omlin and Ivan Utkin

# Usage
Activate the environment using
```
(v1.0) pkg> activate .

(SomeProject) pkg> instantiate
```

NOTE: When you get an error that says
` Module XY with build ID 12345 is missing from the cache.` it is because every worker [tries to precomile the packege on it's own](https://stackoverflow.com/questions/55410326/module-does-not-support-precompilation-but-is-imported-by-a-module-that-does). Fix it by prcomiling beforehand using `]` precompile.

The homework for every week can be found in the corresponding folder

If you have any questions write a mail to `cgutsche@ethz.ch`

## TODO
- [X] You'll find a version of the PorousConvection_3D_xpu.jl code in the solutions folder on Polybox after exercises deadline if needed to get you started.
- [ ] Create a multi-xPU version of your thermal porous convection 3D xPU code you finalised in lecture 7
  - [X] Copy your working PorousConvection_3D_xpu.jl code developed for the exercises in Lecture 7 and rename it PorousConvection_3D_multixpu.jl.
  - [X] Add at the beginning of the code ```using ImplicitGlobalGrid
import MPI```
  - [X] Also add global maximum computation using MPI reduction function
  - [X] max_g(A) = (max_l = maximum(A); MPI.Allreduce(max_l, MPI.MAX, MPI.COMM_WORLD))
  - [X] In the # numerics section, initialise the global grid right after defining nx,ny,nz and use now global grid nx_g(),ny_g() and nz_g() for defining maxiter and ncheck, as well as in any other places when needed.
  - [X] Modify the temperature initialisation using ImplicitGlobalGrid's global coordinate helpers (x_g, etc...), including one internal boundary condition update (update halo):
  - [X] Prepare for visualisation, making sure only me==0 creates the output directory. Also, prepare an array for storing inner points only (no halo) T_inn as well as global array to gather subdomains T_v
  - [X] Use the max_g function in the timestep dt definition (instead of maximum) as one now needs to gather the global maximum among all MPI processes.
  - [X] Moving to the time loop, add halo update function update_halo! after the kernel that computes the fluid fluxes. You can additionally wrap it in the @hide_communication block to enable communication/computation overlap (using b_width defined above)
  - [X] Apply a similar step to the temperature update, where you can also include boundary condition computation as following (⚠️ no other construct is currently allowed)
  - [X] Use now the max_g function instead of maximum to collect the global maximum among all local arrays spanning all MPI processes.
  - [X] Make sure all printing statements are only executed by me==0 in order to avoid each MPI process to print to screen, and use nx_g() instead of local nx in the printed statements when assessing the iteration per number of grid points.
  - [X] Update the visualisation and output saving part
  - [X] Finalise the global grid before returning from the main function
  - [X] Make sure to have set following parameters:
    ```lx,ly,lz    = 40.0,20.0,20.0
    Ra          = 1000
    nz          = 63
    nx,ny       = 2*(nz+1)-1,nz
    b_width     = (8,8,4) # for comm / comp overlap
    nt          = 500
    nvis        = 50

    ```
  - [X] Then, launch the script on Piz Daint on 8 GPU nodes upon adapting the the runme_mpi_daint.sh or sbatch sbatch_mpi_daint.sh scripts (see here) using CUDA-aware MPI 🚀
    - [X] The final 2D slice (at ny_g()/2) produced should look similar as the figure depicted in Lecture 9.
  - [X] Now that you made sure the code runs as expected, launch PorousConvection_3D_multixpu.jl for 2000 steps on 8 GPUs at higher resolution (global grid of 508x252x252) setting:

    ```
    nz          = 127
    nx,ny       = 2*(nz+1)-1,nz
    nt          = 2000
    nvis        = 100    
    ```

    -[X] Use sbtach command to launch a non-interactive job which may take about 5h30-6h to execute.
    - [X] Produce a figure or animation showing the final stage of temperature distribution in 3D and add it to a new section titled ## Porous convection 3D MPI in the PorousConvection project subfolder's README.md. You can use the Makie visualisation helper script from Lecture 7 for this purpose (making sure to adapt the resolution and other input params if needed).

- [X] Keep it xPU compatible using ParallelStencil.jl

- [X] Deploy it on multiple xPUs using ImplicitGlobalGrid.jl

### Exercise 2 - Automatic Documentation
- [X] write some documentation
  - [X] using [doc-strings]([https://](https://docs.julialang.org/en/v1/manual/documentation/))
    - [X] Add doc-string to the functions of following scripts:
      - PorousConvection_3D_xpu.jl
      - PorousConvection_3D_multixpu.jl
  - [X] using Literate.jl
    - [X] Add to the PorousConvection folder a Literate.jl script called bin_io_script.jl that contains and documents following save_array and load_array functions you may have used in your 3D script
    - [X] Add to the bin_io_script.jl a main() function
    - [X] Make the Literate-based workflow to automatically build on GitHub using GitHub Actions. For this, you need to add to the .github/workflow folder the Literate.yml script from the lecuture

## Project
Make sure to have following items in your private GitHub repository:
- [X] a PorousConvection folder containing the structure proposed in Lecture 7
- [X] the 2D and 3D scripts from Lecture 7
- [X] the CI set-up to test the 2D and 3D porous convection scripts
- [X] a lecture_8 folder (different from the PorousConvection folder) containing the codes, README.md and material listed in Exercises - Lecture 8
  - [X] Why does weak scaling not work?
    - [X] Maybe the simulation time is too short
  - [X] strong scaling from `lecture8/src/ex2/t4/l8_diffusion_2D_pref_multixpu._SC.jl`
    - [ ] log scale strong scaling
    - [ ] ideal problem size for strong scaling
- [X] the 3D multi-xPU thermal porous convection script and output as per directions from Exercises - Lecture 9.

In addition enhance the README.md within the PorousConvection folder to include:
- [X] a short motivation/introduction
- [X] concise information about the equations you are solving
- [ ] concise information about the numerical method and implementation
- [ ] the results, incl. figures with labels, captions, etc...
- [ ] a short discussion/conclusion section about the performed work, results, and outlook
