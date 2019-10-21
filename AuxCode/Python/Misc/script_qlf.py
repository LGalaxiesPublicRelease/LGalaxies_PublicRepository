from matplotlib.pyplot import *

# Plots my MR and Scott's MRII on the same graph

hubble=0.673
# Parameters for Hen14
boxside_Hen14=480.28  # Units Mpc/h 
datafile_Hen14='data/Hen14_sfh2/gal_Hen14.pkl'
# Parameters for Guo10
boxside_Guo10=510  # Units Mpc/h 
datafile_Guo10='data/Guo10_sfh2/gal_Guo10.pkl'

# Define limits of plot
# Because I want to use these as plotting data, they need to be arrays, 
# not lists.  The following seems to work
binperdex=10
xrange=np.array([7,13])
nbin=(xrange[1]-xrange[0])*binperdex

#--------------------------------------------------------------------

# Conversion factor from Msun/yr to Lsun
epsilon=0.1
# Not astropy on apollo so need to write numbers in by hand:-(
factor=epsilon*2.0e30*3.00e8**2/(3.16e7*3.85e26)

#--------------------------------------------------------------------

try:
    qlf_Hen14
except:
    pickleFile=datafile_Hen14
    execfile('script_unpickle.py')
    # Mass in units of Msun/h^2
    qlf_Hen14=(gals['QuasarAccretionRate']+gals['RadioAccretionRate'])*factor
try:
    qlf_Guo10
except:
    pickleFile=datafile_Guo10
    execfile('script_unpickle.py')
    # Mass in units of Msun/h^2
    qlf_Guo10=(gals['QuasarAccretionRate']+gals['RadioAccretionRate'])*factor

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
nobj,bins,junk=hist(np.log10(qlf_Hen14), bins=nbin, range=xrange, log=True)
y_Hen14=nobj/(boxside_Hen14)**3*binperdex
nobj,bins,junk=hist(np.log10(qlf_Guo10), bins=nbin, range=xrange, log=True)
y_Guo10=nobj/(boxside_Guo10)**3*binperdex

# Plot at centre of bins
x=0.5*(bins[:-1]+bins[1:])

# Plot
close()
semilogy(x,y_Hen14,'+r',x,y_Guo10,'xb')
axis([xrange[0],xrange[1],1e-7,0.1])
xlabel(r'$\log_{10}(L/L_\odot)$')
ylabel(r'$\log_{10}(N/(\mathrm{dex}\ (h^{-1}\mathrm{Mpc})^3)$')
legend(['Hen14','Guo10'],1)
grid(True)
show()
savefig('qlf.png')

