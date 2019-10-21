#!/bin/tcsh
#
#PBS -S /bin/tcsh
#PBS -N L-galaxies
#PBS -l nodes=1
#PBS -j oe
#PBS -o  ./L-galaxies_output.dat

# Load required modules
source /etc/profile.d/modules.csh
module load openmpi/intel_11.1/1.4.3

# Change to directory where job was submitted
cd $PBS_O_WORKDIR

# Run the program
$MPIROOT/bin/mpirun ./L-Galaxies input/input.par 
