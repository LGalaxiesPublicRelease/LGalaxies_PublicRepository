#!/bin/tcsh                                                                                                                                                                           

#BSUB -L /bin/tcsh                                                                                                                                                                    
#BSUB -q bench1                                                                                                                                                                       
#BSUB -n 12                                                                                                                                                                           
#BSUB -J senna_1                                                                                                                                                                      
#BSUB -oo ./job.log                                                                                                                                                                   
#BSUB -eo ./job.err                                                                                                                                                                   
#BSUB -P durham                                                                                                                                                                       
##BSUB -R "span[hosts=1]"
#BSUB -R "span[ptile=4 ]"  
#BSUB -W 0:30 
 
module load gsl
module load platform_mpi
set MPI_ROOT='/opt/platform_mpi/'


echo "job = $LSB_JOBID"

$MPIROOT/bin/mpirun -IBV ./L-Galaxies ./input/MCMC_inputs/input_mcmc_MR_plus_MRII_W1_PLANCK.par>prog.out
