from matplotlib.pyplot import *

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

# Define limits of plot
# Because I want to use these as plotting data, they need to be arrays, 
# not lists.  The following seems to work
xrange=np.array([5,11])
binperdex=4
nbin=binperdex*(xrange[1]-xrange[0])

#--------------------------------------------------------------------

# Loading in the galaxy data is slow, so we only want to do it once.
# This checks to see whether or not we have done so.
# To force a reread do "del gals"
# Load in Hen14
try:
    BHMass_Hen14
except:
    pickleFile=datafile_Hen14
    execfile('script_unpickle.py')
    BHMass_Hen14=gals['BlackHoleMass']*1e10/hubble_Hen14
# Load in Guo10
try:
    BHMass_Guo10
except:
    pickleFile=datafile_Guo10
    execfile('script_unpickle.py')
    BHMass_Guo10=gals['BlackHoleMass']*1e10/hubble_Guo10

# Put into bins and normalise to number per unit volume in Mpc^3 per dex
#Hen14
nobj,bins,junk=hist(np.log10(BHMass_Hen14), bins=nbin, range=xrange, log=True)
y_Hen14=nobj/(boxside_Hen14/hubble_Hen14)**3*binperdex
#Guo10
nobj,bins,junk=hist(np.log10(BHMass_Guo10), bins=nbin, range=xrange, log=True)
y_Guo10=nobj/(boxside_Guo10/hubble_Guo10)**3*binperdex

# Plot at centre of bins
x=0.5*(bins[:-1]+bins[1:])

# Plot
close()
semilogy(x,y_Hen14,'*r',x,y_Guo10,'*b')
xlabel(r'$\log_{10}(M_\mathrm{BH}/M_\odot)$')
ylabel(r'$\log_{10}(N_\mathrm{BH}/(\mathrm{dex\ Mpc}^3)$')
grid(True)
legend(['Hen14','Guo10'])
show()
savefig('bhmf.png')

