import h5py
from matplotlib.pyplot import *
import re

# Haven't worked out why but show() does not replot the data.
# So need to close the plotting window between calls.

# Parameters for Hen14
hubble=0.673
boxside=480.28  # Units Mpc/h 
datafile='data_pkl/SA_output_5.h5'

# Define limits of plot
xmin,xmax=[1e9,1e12]
ymin,ymax=[1e5,1e12]

#--------------------------------------------------------------------

fin=h5py.File(datafile,'r')

close()

# Loop over snapshots
for dir in fin:
    if not re.match(r'\d\d',dir): break
    StellarMass=fin[dir]['StellarMass']*1e10/hubble
    BlackHoleMass=fin[dir]['BlackHoleMass']*1e10/hubble
    loglog(StellarMass,BlackHoleMass/StellarMass,'.',label=str(dir))

show()
legend()

fin.close()


