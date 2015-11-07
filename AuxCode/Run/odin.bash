################################
#$ -S /bin/bash
#$ -N Senna_02
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 16 
##$ -pe impi_hydra 512
#$ -pe impi_hydra 64
### up to 24 hours
#$ -l h_rt=24:00:00
#$ -M bhenriques@mpa-garching.mpg.de

module load intel
module load impi
module load fftw/2.1.5
module load gsl/1.14
module load hdf5-serial

#mpiexec -np $NSLOTS ./L-Galaxies ./input/input_Henriques2014_MR_W1_PLANCK.par > prog01.out 
mpiexec -np $NSLOTS ./L-Galaxies ./input/input_MRII_W1_PLANCK.par > prog02.out 
#mpiexec -np $NSLOTS ./L-Galaxies ./input/input_MR_W1_W1.par > prog01.out 
#mpiexec -perhost 2 -np 64 ./L-Galaxies ./input/input_mcmc_halomodel_guo11_MR_W1_W1.par > prog02.out 
#time mpiexec -perhost 4 -np 128 ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_PLANCK.par > prog11.out
#time mpiexec -perhost 4 -np 16 ./L-Galaxies ./input/input_mcmc_halomodel_guo11_MR_W1_W1.par > prog22.out

####$TMPDIR
################################



