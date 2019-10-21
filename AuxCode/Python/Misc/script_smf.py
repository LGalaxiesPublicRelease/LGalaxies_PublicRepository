from matplotlib.pyplot import *

# Needed in python3 to mimic python2 execfile.
# The globals() parameter copies variables between namespaces.
def execfile(file):
    with open(file,'r') as fid:
        exec(fid.read(),globals())

# Haven't worked out why but show() does not replot the data.
# So need to close the plotting window between calls.

# Parameters for Hen15
hubble=0.673
boxside=480.28  # Units Mpc/h 
datafile='data_pkl/new/snap_58.pkl'
nfile=10

# Define limits of plot
# Because I want to use these as plotting data, they need to be arrays, 
# not lists.  The following seems to work
binperdex=10
xrange=np.array([7,12])
nbin=(xrange[1]-xrange[0])*binperdex

#--------------------------------------------------------------------

# Loading in the galaxy data is slow, so we only want to do it once.
# This checks to see whether or not we have done so.
# To force a reread do "del gals"
# Load in data
try:
    mstar
except:
    pickleFile=datafile
    execfile('script_unpickle.py')
    # Mass in units of Msun/h^2
    mstar=(gals['DiskMass']+gals['BulgeMass'])*1e10*hubble

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj,bins,junk=hist(np.log10(mstar), bins=nbin, range=xrange, log=True)
y=nobj*512/(nfile*boxside**3)*binperdex

# Plot at centre of bins
x=0.5*(bins[:-1]+bins[1:])

# Plot
close()
semilogy(x,y,'+r')
axis([12.4,7,10.**(-5.9),10.**0.5])
xlabel(r'$\log_{10}(M_*/h^{-2}M_\odot)$')
ylabel(r'$\log_{10}(N/(\mathrm{dex}\ (h^{-1}\mathrm{Mpc})^3)$')
grid(True)
show()
savefig('smf.png')

