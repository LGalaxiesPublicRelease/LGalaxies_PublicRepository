import h5py
from matplotlib.pyplot import *
import re
import seaborn as sb

# Haven't worked out why but show() does not replot the data.
# So need to close the plotting window between calls.

# Parameters for Hen14
hubble=0.673
boxside=480.28  # Units Mpc/h 
prefix_Hen15='data_pkl/Hen15/SA_output_'
prefix_DI='data_pkl/SA_output_'
postfix='.h5'

# snapshots to plot
#snaps=['30','58']
snaps=['58']

# files to include
files=np.arange(10)

# Define limits of plot
xmin,xmax=[1e9,1e12]
ymin,ymax=[1e5,1e12]

#--------------------------------------------------------------------


close()

fin_Hen15=[]
# Open files for reading
for fid in files:
    fin_Hen15.append(h5py.File(prefix_Hen15+str(fid)+postfix,'r'))
# Loop over snapshots
for dir in snaps:
    # First count how many galaxies there are over all files.
    ngal=0
    for fin in fin_Hen15:
        ngal+=len(fin[dir]['StellarMass'])
    # Next declare space for desired arrays.
    StellarMass=np.empty(ngal)
    BlackHoleMass=np.empty(ngal)
    # Now do second loop to populate arrays.
    igal=0
    for fin in fin_Hen15:
        ngal=len(fin[dir]['StellarMass'])        
        StellarMass[igal:igal+ngal]=fin[dir]['StellarMass']*1e10/hubble
        BlackHoleMass[igal:igal+ngal]=fin[dir]['BlackHoleMass']*1e10/hubble
        igal+=ngal
    assert igal==len(StellarMass)
    loglog(StellarMass,BlackHoleMass/StellarMass,'.',label='Hen15 - '+str(dir))
# Close files
for fin in fin_Hen15: fin.close()

fin_DI=[]
# Open files for reading
for fid in files:
    fin_DI.append(h5py.File(prefix_DI+str(fid)+postfix,'r'))
# Loop over snapshots
for dir in snaps:
    # First count how many galaxies there are over all files.
    ngal=0
    for fin in fin_DI:
        ngal+=len(fin[dir]['StellarMass'])
    # Next declare space for desired arrays.
    StellarMass=np.empty(ngal)
    BlackHoleMass=np.empty(ngal)
    # Now do second loop to populate arrays.
    igal=0
    for fin in fin_DI:
        ngal=len(fin[dir]['StellarMass'])        
        StellarMass[igal:igal+ngal]=fin[dir]['StellarMass']*1e10/hubble
        BlackHoleMass[igal:igal+ngal]=fin[dir]['BlackHoleMass']*1e10/hubble
        igal+=ngal
    assert igal==len(StellarMass)
    loglog(StellarMass,BlackHoleMass/StellarMass,'.',label='DI - '+str(dir))
# Close files
for fin in fin_DI: fin.close()

show()
legend()
savefig('figs/bh_mstar.png')
