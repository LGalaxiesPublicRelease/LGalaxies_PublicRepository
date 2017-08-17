################################                                                                                                                                              
#$ -S /bin/bash                                                                                                                                                               
#$ -j y                                                                                                                                                                       
#$ -cwd                                                                                                                                                                       
#$ -m n                                                                                                                                                                       
### must be a multiple of 12 up to 72                                                                                                                                         
#$ -pe impi4 72                                                                                                                                                                                                                                                                                                                      
### up to 24 hours                                                                                                                                                            
#$ -l h_rt=24:00:00                                                                                                                                                           
#$ -M bhenriques@mpa-garching.mpg.de                                                                                                                                          
module load impi
ulimit -aH
#export I_MPI_DEBUG=2                                                                                                                                                         

###export OMP_NUM_THREADS=$(4)                                                                                                                                                

time mpiexec  -ppn 6 ./L-Galaxies ./input/MCMC_inputs/input_mcmc_MR_plus_MRII_W1_PLANCK.par > prog.out

#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_MR_W1_PLANCK.par  > prog01.out $TMPDIR                                                                             
#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_nifty_MR_PLANCK.par > prog03.out                                                                                          
#time mpiexec -n $NSLOTS ./L-Galaxies ./input//input_guo_MRII_W1_PLANCK.par > prog02.out 
#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_guo11_MR_W1_W1.par > prog01.out    

#time mpiexec -perhost 2 -np 6 ./L-Galaxies ./input/input_mcmc_halomodel_guo11_MR_W1_W1.par > prog01.out                                                                      
#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_PLANCK.par > prog11.out                                                                              
date
####$TMPDIR                                                                                                                                                                   
################################                         
