# Script to call other scripts to read in, mass trim and pickle the data.
# This is the python2 version.

# The parameters that can be changed in this script are those that the user
# will change most frequently.
# The locations of the data and metadata for each model are set in script_read;
# also the desired properties (which cannot be set until after the template is
# read in).

# Optionally set model.  Defaults to Hen15
model='new'
#model='Hen15'
#model='Guo10'
#model='Hen14_MRII'

# Optionally set file range.  Defaults to 0-511
firstfile = 5
lastfile = 5

# Select desired redshifts
#redshifts = 0,0.15,0.25,0.5,1,1.5,2,3,4,5,6
redshifts = 0,

#-----------------------------------------------------------

for redshift in redshifts:
    execfile('script_read.py')
    execfile('script_masstrim.py')
    execfile('script_pickle2.py')
