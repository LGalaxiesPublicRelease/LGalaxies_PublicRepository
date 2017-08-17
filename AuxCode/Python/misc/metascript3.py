# Script to call other scripts to read in, mass trim and pickle the data.
# This is the python3 version.

# The parameters that can be changed in this script are those that the user
# will change most frequently.
# The locations of the data and metadata for each model are set in script_read;
# also the desired properties (which cannot be set until after the template is
# read in).

# Optionally set model.  Defaults to Hen15
#model='Guo11'
#model='Hen15_spin'
#model='Guo11_spin'
#model='Hen14_MRII'
#model='Rob'
model='new'

# Optionally set file range.  Defaults to 0-511
firstfile = 0
lastfile = 9

# Select desired redshifts
#redshifts = 0,0.15,0.25,0.5,1,1.5,2,3,4,5,6
redshifts = 0.,

#-----------------------------------------------------------

# Needed in python3 to mimic python2 execfile.
# The globals() parameter copies variables between namespaces.
def execfile(file):
    with open(file,'r') as fid:
        exec(fid.read(),globals())

for redshift in redshifts:
    execfile('script_read.py')
    execfile('script_masstrim.py')
    execfile('script_pickle3.py')
