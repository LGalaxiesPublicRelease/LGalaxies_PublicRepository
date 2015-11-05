################################
#$ -S /bin/bash
#$ -j y 
#$ -cwd
#$ -m n
### must be a multiple of 12 up to 72
#$ -pe impi4 12
### up to 24 hours
#$ -l h_rt=24:00:00
#$ -M bhenriques@mpa-garching.mpg.de
module load impi
ulimit -aH
#export I_MPI_DEBUG=2

#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_MR_W1_W7.par > prog01.out 

#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_MR_W1_PLANCK.par  > prog01.out $TMPDIR

time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_nifty_MR_PLANCK.par > prog03.out 

#time mpiexec -n $NSLOTS ./L-Galaxies ./input//input_guo_MRII_W1_PLANCK.par > prog02.out 


#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_MR_W1_W1.par > prog01.out 

#time mpiexec -np $NSLOTS ./L-Galaxies ./input/input_Henriques2014_MR_W1_W1.par $TMPDIR > prog01.out 

#time mpiexec -np $NSLOTS ./L-Galaxies ./input/input_Henriques2013_MR_W1_W7.par > prog01.out

#mpiexec -np $NSLOTS ./L-Galaxies ./input/input_guo13_MR_W1_W7.par > prog01.out
#mpiexec -np $NSLOTS ./L-Galaxies ./input/input_MR_W1_W1.par > prog01.out

#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_guo11_MR_W1_W1.par > prog01.out


 
#time mpiexec -perhost 2 -np 6 ./L-Galaxies ./input/input_mcmc_halomodel_guo11_MR_W1_W1.par > prog01.out 
#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_PLANCK.par > prog11.out
date
####$TMPDIR
################################




