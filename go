### run using ./go Nprocessors Simulation: 
### eg: ./go 16 MRII
### make sure you make go an executable file first: chmod u+x ./go

#$ -S /bin/bash

if [ $# -ne 2 ]; then
    echo wrong input - check the batch file!
    echo use: "go <processors> <simulation>"
    exit
else

	processors=$1
	simulation=$2
	var1="./input/input_Henriques2015_"
	var2="_W1_PLANCK.par"
	echo $simulation
	
	echo running SA model: on $processors processors
	mpirun -np $processors ./L-Galaxies $var1$simulation$var2

fi

echo finished 
