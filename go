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
	#var1="./input/Guo_inputs/input_"
	#var2="_W1_W1.par"
	var1="./input/Hen15_inputs/input_"
	#var1="./input/input_"
	var2="_W1_PLANCK.par"
	echo $simulation
	
	echo running SA model: on $processors processors
	nice -n -10 mpirun -np $processors ./L-Galaxies $var1$simulation$var2

fi

echo finished 
