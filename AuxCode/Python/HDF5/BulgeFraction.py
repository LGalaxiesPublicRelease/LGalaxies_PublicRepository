
# coding: utf-8

# # Bulge fraction

# In[5]:


# Imports
import astropy.constants as c
import astropy.units as u
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.table import Table
get_ipython().run_line_magic('matplotlib', 'inline')
matplotlib.rc('font',size=20)


# In[6]:


# Parameters
Compare=True

# Simulation data
outputDirRoot='../../../output/'
outputDir1='cooling_old/'
outputDir2='cooling_new/'
#outputDir='../../../output/merge/hack1/'
obsDir='../../../../Obsdata/'
filePrefix='SA_output_'
filePostfix='.h5'
snap='58'
redshift=0. # Read from file
firstFile=0
lastFile=9
maxFile=512

# Observational data
#Loading the 3 different bulge fractions obs files
obsdir='../../../../Obsdata/'
obsBulgeFile='conselice2006_bulge_fract.txt'
obsDiskFile='conselice2006_disk_fract.txt'
obsIrrFile='conselice2006_irr_fract.txt'
obsHubble=0.70

# Divisions between bulge classes
logRatio1=-0.154902
logRatio2=-2

# Bins for hostogram and plotting
binwidth=0.25
xrange=np.array([8.8,11.5])
bins=np.arange(xrange[0],xrange[1]+0.001,binwidth)

# Simulation parameters (read from file!)
hubble=0.673
boxside=480.28  # Units Mpc/h 


# In[7]:


# First determine the size of the arrays that we need to hold the data
nGal1=0
for iFile in range(firstFile,lastFile+1):
    # The following line closes the file at the end of the loop
    with h5py.File(outputDirRoot+outputDir1+filePrefix+'%i'%iFile+filePostfix,'r') as f:
        nGal1+=len(f[snap])
print('nGal1=',nGal1)

# Declare numpy arrays to hold the data
Type1=np.empty(nGal1)
BulgeMass1=np.empty(nGal1)
DiskMass1=np.empty(nGal1)
StellarMass1=np.empty(nGal1)

# Now read in the data
iGal=0
for iFile in range(firstFile,lastFile+1):
    # The following line closes the file at the end of the loop
    with h5py.File(outputDirRoot+outputDir1+filePrefix+'%i'%iFile+filePostfix,'r') as f:
        nGalFile=len(f[snap])
        Type1[iGal:iGal+nGalFile]=f[snap]['Type']
        BulgeMass1[iGal:iGal+nGalFile]=f[snap]['BulgeMass']
        DiskMass1[iGal:iGal+nGalFile]=f[snap]['DiskMass']
        StellarMass1[iGal:iGal+nGalFile]=f[snap]['StellarMass']
        iGal+=nGalFile
assert np.all(abs(StellarMass1-BulgeMass1-DiskMass1)<1e-5*StellarMass1)

# Put into observer units and add scatter to stellar mass estimate
offset=10+np.log10(hubble)
logBulge1=np.log10(BulgeMass1)+offset
logDisk1=np.log10(DiskMass1)+offset
logStellarMass1=np.log10(StellarMass1)+offset
logStellarMassObs1=logStellarMass1+np.random.randn(nGal1)*0.08*(1+redshift)

if Compare:
    nGal2=0
    for iFile in range(firstFile,lastFile+1):
        # The following line closes the file at the end of the loop
        with h5py.File(outputDirRoot+outputDir2+filePrefix+'%i'%iFile+filePostfix,'r') as f:
            nGal2+=len(f[snap])
    print('nGal2=',nGal2)

    # Declare numpy arrays to hold the data
    Type2=np.empty(nGal2)
    BulgeMass2=np.empty(nGal2)
    DiskMass2=np.empty(nGal2)
    StellarMass2=np.empty(nGal2)

    # Now read in the data
    iGal=0
    for iFile in range(firstFile,lastFile+1):
        # The following line closes the file at the end of the loop
        with h5py.File(outputDirRoot+outputDir2+filePrefix+'%i'%iFile+filePostfix,'r') as f:
            nGalFile=len(f[snap])
            Type2[iGal:iGal+nGalFile]=f[snap]['Type']
            BulgeMass2[iGal:iGal+nGalFile]=f[snap]['BulgeMass']
            DiskMass2[iGal:iGal+nGalFile]=f[snap]['DiskMass']
            StellarMass2[iGal:iGal+nGalFile]=f[snap]['StellarMass']
            iGal+=nGalFile
    assert np.all(abs(StellarMass2-BulgeMass2-DiskMass2)<1e-5*StellarMass2)

    # Put into observer units and add scatter to stellar mass estimate
    offset=10+np.log10(hubble)
    logBulge2=np.log10(BulgeMass2)+offset
    logDisk2=np.log10(DiskMass2)+offset
    logStellarMass2=np.log10(StellarMass2)+offset
    logStellarMassObs2=logStellarMass2+np.random.randn(nGal2)*0.08*(1+redshift)


# In[8]:


# Initialise plot
if Compare:
    plt.figure(figsize=[16,6])
    plt.subplot(1,2,1)
else:
    plt.figure(figsize=[8,6])

# Put galaxies into bins
indBin=np.digitize(logStellarMassObs1,bins)
nBin=len(bins)-1
x=np.empty(nBin)
yBulge=np.empty(nBin)
yDisk=np.empty(nBin)

yIrr=np.empty(nBin)
# Loop over bins, couting fractions in each class
for iBin in range(nBin):
    x[iBin]=0.5*(bins[iBin]+bins[iBin+1])
    indThisBin=np.where(indBin==iBin+1)[0]
    allBin=len(indThisBin)
    # Bulges
    yBulge[iBin]=len(np.where((logBulge1[indThisBin]-logStellarMass1[indThisBin])>logRatio1)[0])/allBin
    # Disks
    yIrr[iBin]=len(np.where((logBulge1[indThisBin]-logStellarMass1[indThisBin])<logRatio2)[0])/allBin
    # Intermediates
    yDisk[iBin]=1.-yBulge[iBin]-yIrr[iBin]
plt.xlabel(r'$M_*/h^{-2}\mathrm{M}_\odot$')
plt.ylabel(r'Fraction')
plt.ylim([0,1])
plt.plot(x,yBulge,'r-',label='Bulge')
plt.plot(x,yDisk,'b-',label='Disk')
plt.plot(x,yIrr,'g-',label='Irr')
plt.legend()
plt.title(outputDir1[:-1])
# Plot observations
plt.gca().set_prop_cycle(None) # Doesn't make errorbar use the correct colours
obsBulge=np.loadtxt(obsdir+obsBulgeFile)
obsMass=obsBulge[:,0]+2*np.log10(obsHubble)
plt.errorbar(obsMass,obsBulge[:,1],yerr=obsBulge[:,2],marker='o',linestyle='None',color='r')
obsDisk=np.loadtxt(obsdir+obsDiskFile)
obsMass=obsDisk[:,0]+2*np.log10(obsHubble)
plt.errorbar(obsMass,obsDisk[:,1],yerr=obsDisk[:,2],marker='o',linestyle='None',color='b')
obsIrr=np.loadtxt(obsdir+obsIrrFile)
obsMass=obsIrr[:,0]+2*np.log10(obsHubble)
plt.errorbar(obsMass,obsIrr[:,1],yerr=obsIrr[:,2],marker='o',linestyle='None',color='g')

if Compare:
    plt.subplot(1,2,2)

    # Put galaxies into bins
    indBin=np.digitize(logStellarMassObs2,bins)
    nBin=len(bins)-1
    x=np.empty(nBin)
    yBulge=np.empty(nBin)
    yDisk=np.empty(nBin)

    yIrr=np.empty(nBin)
    # Loop over bins, couting fractions in each class
    for iBin in range(nBin):
        x[iBin]=0.5*(bins[iBin]+bins[iBin+1])
        indThisBin=np.where(indBin==iBin+1)[0]
        allBin=len(indThisBin)
        # Bulges
        yBulge[iBin]=len(np.where((logBulge2[indThisBin]-logStellarMass2[indThisBin])>logRatio1)[0])/allBin
        # Disks
        yIrr[iBin]=len(np.where((logBulge2[indThisBin]-logStellarMass2[indThisBin])<logRatio2)[0])/allBin
        # Intermediates
        yDisk[iBin]=1.-yBulge[iBin]-yIrr[iBin]
    plt.xlabel(r'$M_*/h^{-2}\mathrm{M}_\odot$')
    plt.ylabel(r'Fraction')
    plt.ylim([0,1])
    plt.plot(x,yBulge,'r-',label='Bulge')
    plt.plot(x,yDisk,'b-',label='Disk')
    plt.plot(x,yIrr,'g-',label='Irr')
    plt.legend()
    plt.title(outputDir2[:-1])
    # Plot observations
    plt.gca().set_prop_cycle(None) # Doesn't make errorbar use the correct colours
    obsBulge=np.loadtxt(obsdir+obsBulgeFile)
    obsMass=obsBulge[:,0]+2*np.log10(obsHubble)
    plt.errorbar(obsMass,obsBulge[:,1],yerr=obsBulge[:,2],marker='o',linestyle='None',color='r')
    obsDisk=np.loadtxt(obsdir+obsDiskFile)
    obsMass=obsDisk[:,0]+2*np.log10(obsHubble)
    plt.errorbar(obsMass,obsDisk[:,1],yerr=obsDisk[:,2],marker='o',linestyle='None',color='b')
    obsIrr=np.loadtxt(obsdir+obsIrrFile)
    obsMass=obsIrr[:,0]+2*np.log10(obsHubble)
    plt.errorbar(obsMass,obsIrr[:,1],yerr=obsIrr[:,2],marker='o',linestyle='None',color='g')

if Compare:
    plt.savefig('figs/BulgeFraction_'+outputDir1[:-1]+'_'+outputDir2[:-1]+'.png')
else:
    plt.savefig('figs/BulgeFraction_'+outputDir1[:-1]+'.png')

