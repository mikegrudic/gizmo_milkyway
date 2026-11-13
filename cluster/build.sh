#!/bin/bash
# Build GIZMO with a run's Config.sh. Usage: build.sh <GIZMO source dir> <Config.sh>
# The modules and SYSTYPE are for the Genoa nodes of Flatiron's Rusty cluster: change both for your machine.
set -e
module --force purge
module load modules openmpi gsl hdf5/mpi-1.12.3
export SYSTYPE=${SYSTYPE:-RUSTY_GENOA}
cp "$2" "$1/Config.sh"
cd "$1"
make clean > /dev/null
make -j 32 > make.log 2>&1 || { tail -30 make.log; exit 1; }
ls -la GIZMO
