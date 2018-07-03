################################
#$ -S /bin/bash
#$ -N Senna_02
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 16 
#$ -pe impi_hydra 128
### up to 24 hours
#$ -l h_rt=24:00:00
#$ -M bhenriques@mpa-garching.mpg.de

module load intel
module load impi
module load fftw/2.1.5
module load gsl/1.14
module load hdf5-serial

mpiexec -np $NSLOTS ./L-Galaxies ./input/input_Henriques2014_MRII_W1_PLANCK.par > prog02.out 
##mpiexec -np $NSLOTS ./L-Galaxies ./input/input_MR_W1_PLANCK.par > prog01.out 
##mpiexec -np $NSLOTS ./L-Galaxies ./input/input_MR_w1_planck.par > prog.out
##time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_PLANCK.par > prog00.out

####$TMPDIR
################################
