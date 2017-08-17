#!/usr/local/bin/tcsh
#$ -cwd
#$ -q amd
#$ -m be
#$ -l h_vmem=2000M
#$ -M myemailaddress@mpa-garching.mpg.de
#$ -N galaxies

./L-Galaxies input/input.par

