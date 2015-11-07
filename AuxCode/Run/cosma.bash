#!/bin/tcsh

#BSUB -L /bin/tcsh
#BSUB -q cosma
#BSUB -n 120
#BSUB -J senna_1
#BSUB -oo ./job.log
#BSUB -eo ./job.err
#BSUB -P durham
#BSUB -R "span[ptile=12]"

module purge

#module load gsl
module load gsl/c4/1.15
#module load platform_mpi/intel_2012.0.032/8.2.0
module load platform_mpi/c4/intel_11.1/8.2.0

#module load platform_mpi/8.2.1
#set MPI_ROOT='/opt/platform_mpi/'



echo "job = $LSB_JOBID"

$MPIROOT/bin/mpirun -IBV ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_PLANCK.par>prog.out
#$MPIROOT/bin/mpirun ./L-Galaxies ./input/input_mcmc_MR_plus_MRII.par>prog.out




#/gpfs/COSMA/openmpi/pgi_11.10/1.4.4//bin/mpirun ./L-Galaxies ./input/input_mcmc.par>prog.out

