#!/bin/bash

#! System to be used
#SBATCH -p knl
#! Number of nodes to be used
#SBATCH -N 1
#! Number of MPI tasks allocated to each node
#SBATCH --ntasks-per-node=8
#! Wall clock time required
#SBATCH --time=00:30:00
#! Name of job
#SBATCH -J hdf5off

srun ./Main_ScalarField3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ICPX2022.0.ex params.txt
