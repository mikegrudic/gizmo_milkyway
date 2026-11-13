#!/bin/bash
# One run per node as a Slurm job array. Usage: sbatch --array=0-<N-1> run_grid.sbatch <run dir 0> ... <run dir N-1>
# Submitting again resumes every run from its restart files.
# The directives and modules are for Flatiron's Rusty cluster: change them for your site.
#SBATCH --job-name=mw_grid
#SBATCH --partition=cca
#SBATCH --constraint=genoa&rocky9
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
RUNS=("$@")
cd "${RUNS[$SLURM_ARRAY_TASK_ID]}" || exit 1
unset ${!SLURM_@}   # mpirun would otherwise size itself from the one-task allocation
module --force purge
module load modules openmpi gsl hdf5/mpi-1.12.3
export OMP_NUM_THREADS=1
flag=
if [ -e output/restartfiles/restart.0 ]; then flag=1; mv gizmo.out gizmo.out.$(date +%s); fi
echo "$(hostname) $PWD restart_flag=${flag:-0}"
mpirun -np 96 --bind-to core ./GIZMO mw.params $flag > gizmo.out 2> gizmo.err
