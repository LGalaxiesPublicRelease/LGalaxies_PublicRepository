#!/usr/local/bin/tcsh
#$ -cwd
#$ -pe mpich 4
#$ -m be
#$ -M myemailaddress@mpa-garching.mpg.de
#$ -N galaxies

mpirun -np $NSLOTS ./L-Galaxies input/input.par

