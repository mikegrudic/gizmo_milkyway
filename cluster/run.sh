#!/bin/bash
# Run GIZMO in the current directory on one node. Usage: run.sh <GIZMO binary> <params> [ranks] [restart flag]
# Restart flag 1 resumes from output/restartfiles. The modules are for Flatiron's Rusty cluster.
unset ${!SLURM_@}   # inside an allocation, mpirun would otherwise size itself from the Slurm task count
module --force purge
module load modules openmpi gsl hdf5/mpi-1.12.3
export OMP_NUM_THREADS=1
hostname
mpirun -np ${3:-96} --bind-to core "$1" "$2" $4 > gizmo.out 2> gizmo.err
