#!/bin/bash
#SBATCH -J CSUr2
#SBATCH -o CSUr2.log
#SBATCH -e CSUr2.err
#SBATCH -N 1 -n 107
mpirun -np 107 Mod3DMT.iMPI.BICG -I NLCG Modular_NLCG_072.rho 27-freq-run-watex01.dat inv.ctrl fwd.ctrl 27-freq-run-watex01.cov
