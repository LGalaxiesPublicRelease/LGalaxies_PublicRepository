#!/bin/tcsh                                                                                                                                                    

#BSUB -L /bin/tcsh                                                                                                                                             
#BSUB -q cosma                                                                                                                                                 
#BSUB -n 120                                                                                                                                                   
#BSUB -J senna_0                                                                                                                                               
#BSUB -oo ./job.log                                                                                                                                            
#BSUB -eo ./job.err                                                                                                                                            
#BSUB -P durham                                                                                                                                                
#BSUB -R "span[ptile=4]"                                                                                                                                       


module load gsl
module load platform_mpi

set MPI_ROOT='/opt/platform_mpi/'


echo "job = $LSB_JOBID"

$MPIROOT/bin/mpirun -IBV ./L-Galaxies ./input/MCMC_inputs/input_mcmc_MR_plus_MRII_W1_PLANCK.par>prog.out
#$MPIROOT/bin/mpirun ./L-Galaxies ./input/input_mcmc_MR_plus_MRII.par>prog.out                                                                                 

