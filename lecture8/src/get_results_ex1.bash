# !bin/bash -l

scp -r daint:/users/class202/github-repo/pde-on-gpu-gutsche/lecture8/docs/l8ex1t3_out/ ../docs/
julia --project=.. l8_vizme2Dgpu_mpi.jl ../docs/l8ex1t3_out/mpi2Dgpu_out_C_