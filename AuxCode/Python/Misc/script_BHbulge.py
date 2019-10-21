from matplotlib.pyplot import *
import numpy as np


# Haven't worked out why but show() does not replot the data.
# So need to close the plotting window between calls.

# Parameters for Hen14
hubble_Hen14=0.673
boxside_Hen14=480.28  # Units Mpc/h 
datafile_Hen14='data/Hen14_sfh2/gal_Hen14.pkl'

# Parameters for Guo10
hubble_Guo10=0.73
boxside_Guo10=500.  # Units Mpc/h 
datafile_Guo10='data/Guo10_sfh2/gal_Guo10.pkl'

# Plot scale
xrange=np.array([10**8.5,10**12.5])
yrange=np.array([1e5,10**10.5])

#--------------------------------------------------------------------

# Loading in the galaxy data is slow, so we only want to do it once.
# This checks to see whether or not we have done so.
# Load in Hen14
try:
    BHMass_Hen14
    BulgeMass_Hen14
except:
    pickleFile=datafile_Hen14
    execfile('script_unpickle.py')
    BHMass_Hen14=gals['BlackHoleMass']*1e10
    BulgeMass_Hen14=gals['BulgeMass']*1e10*hubble_Hen14
# Load in Guo10
try:
    BHMass_Guo10
    BulgeMass_Guo10
except:
    pickleFile=datafile_Guo10
    execfile('script_unpickle.py')
    BHMass_Guo10=gals['BlackHoleMass']*1e10
    BulgeMass_Guo10=gals['BulgeMass']*1e10*hubble_Guo10

# Produce log-log plot
clf() # Clear existing plot
#loglog(BulgeMass_Hen14,BHMass_Hen14,'r,',BulgeMass_Guo10,BHMass_Guo10,'b,')
loglog(BulgeMass_Guo10,BHMass_Guo10,'b,',BulgeMass_Hen14,BHMass_Hen14,'r,')
xlabel(r'$M_\mathrm{Bulge}/h^{-2}M_\odot$') # Can use latex in labels
ylabel(r'$M_\mathrm{BH}/h^{-1}M_\odot$')
# We only believe galaxies more massive than 1e9 Msun
xlim(xrange)
ylim(yrange)
# Plot a linear line just for comparison
ratio=1e-3
loglog(xrange,ratio*xrange)
#legend(['Hen14','Guo10'],2)
legend(['Guo10','Hen14'],2)
#show()
savefig('bhbulge.png')
