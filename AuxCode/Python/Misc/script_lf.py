# Python 3 script to plot stellar mass function

import matplotlib
from matplotlib.pyplot import *
import numpy as np

# Magnitude that is desired: 0 - u, 1 - g, 2 - r, 3 - i, 4 - v
mag_type=2

# Parameters for Hen14
hubble=0.673

# MR
boxside_MR=480.28  # Units Mpc/h 
datafile_MR='data/Hen14_sfh2/snap_58.pkl'

# Parameters for Hen14 MRII
boxside_MRII=96.0558  # Units Mpc/h 
datafile_MRII='data/Hen14_MRII/snap_58.pkl'

# Define limits of plot
# Because I want to use these as plotting data, they need to be arrays, 
# not lists.  The following seems to work
xrange=np.array([-27,-17])
binperdex=5
nbin=binperdex*(xrange[1]-xrange[0])

#--------------------------------------------------------------------

if mag_type==0: Xlabel='u'
if mag_type==1: Xlabel='g'
if mag_type==2: Xlabel='r'
if mag_type==3: Xlabel='i'
if mag_type==4: Xlabel='z'
Ylabel=r'$N/$dex$\,Mpc^3$'
Title='Luminosity function'
pngfile='lf.png'

#--------------------------------------------------------------------

# Loading in the galaxy data is slow, so we only want to do it once.
# This checks to see whether or not we have done so.
# Load in MR
try:
    lum_MR
except:
    pickleFile=datafile_MR
    exec(open('script_unpickle.py').read())
    lum_MR=gals['Mag'][:,mag_type]
# Load in MR
try:
    lum_MRII
except:
    pickleFile=datafile_MRII
    exec(open('script_unpickle.py').read())
    lum_MRII=gals['Mag'][:,mag_type]

#--------------------------------------------------------------------

# Put into bins and normalise to number per unit volume in Mpc^3 per dex
#MR
nobj,bins,junk=hist(lum_MR, bins=nbin, range=xrange, log=True)
y_MR=nobj*binperdex/(boxside_MR/hubble)**3
#MRII
nobj,bins,junk=hist(lum_MRII, bins=nbin, range=xrange, log=True)
y_MRII=nobj*binperdex/(boxside_MRII/hubble)**3

# Plot at centre of bins
x=0.5*(bins[:-1]+bins[1:])

# Plot
clf()
fig1=semilogy(x,y_MR,'dr',x,y_MRII,'sb')
grid(True)
xlabel(Xlabel)
ylabel(Ylabel)
title(Title)
legend(['MR','MRII'])
show()
savefig(pngfile)

