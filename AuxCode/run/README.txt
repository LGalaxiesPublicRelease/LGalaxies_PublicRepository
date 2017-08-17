GO
==

The shell script go is an example script to run the code provided by Darren
Croton. Here some information he provided:

Firstly, edit the script and change "base_output_dir" to point to the place
where your galaxies are to be stored.  To start with you will perhaps want to
just use "./output/$1".  We will use this value here to illustrate, although
you may want to change it to "/scratch/$1" or whatever if disk space is tight.
Either way, this should coordinate with "OutputDir" in the input parameter
file.

> go test

This creates the directory ./output/test/ (*** which should match "OutputDir"
in the input.par file ***) and then runs the model using the input file
./input/input.par.  Everything that is written to the screen will also be
written in the file test.out and put in the just created ./output/test/
directory (this way you can keep track of your different model settings, since
the parameter choices are written to the screen when the code is run).  Of
course the semi-analytic galaxies will also appear here for each snapshot you
have requested.

An additional feature has been added.  When on the RZG computers the "go"
script will run the code on multiple processors.  You can specify as many as
you want, upto 32, with:

> go test 16

(say).  Each tree file will then be run on a seperate processor.  If you have
more procesors specified than tree files then the left-over processors will
exit straight away.

*** To use the parallel option you must set the -DPARALLEL flag in the
makefile. ***
