################################
#$ -S /bin/bash
#$ -N Senna_Halo_1
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 12 up to 48
#$ -pe impi_hydra 32
### up to 24 hours
#$ -l h_rt=24:00:00
#$ -M bhenriques@mpa-garching.mpg.de

module load intel
module load impi
module load fftw/2.1.5
module load gsl/1.14
module load hdf5-serial

time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_mcmc_guo11_MR_W1_W1.par > prog.out


####$TMPDIR
################################
