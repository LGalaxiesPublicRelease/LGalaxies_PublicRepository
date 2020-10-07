#!/usr/bin/env python
# coding: utf-8

# # Script to test reading in of SFHs and example usage

# In[3]:


# Imports
import astropy.constants as c
import astropy.units as u
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
matplotlib.rc('font',size=20)


# In[31]:


# Parameters
outputDir='../../../output/'
filePrefix='SA_output_'
filePostfix='.h5'
snap='58'
firstFile=5
lastFile=5
# Size of SFH array.  This can probably be read in HDF5 from the dtype
nSFHbins=20
# Plotting parameters
xmin=1e12
xmax=1e15
ymin0=1e9
ymax0=1e14
ymin1=0.
ymax1=0.2


# In[53]:


# import template
#sys.path.insert(0,idatadir)
import sfh_template
# I presume that the following offset is needed because the struct size is not a multiple of 8 bytes.
itemsize_correction=4

# Open data file for binary reading
with open(outputDir+'SFH_Bins','rb') as f:
    Nbin = np.fromfile(f,np.int32,1)[0]
    # Position file pointer at end of header
    f.seek(sfh_template.sfh_struct_dtype.itemsize+itemsize_correction,0)
    # Create record array to hold data and read it in
    #SFHbins=np.array(Nbin,dtype=sfh_struct_dtype)
    SFHbins=np.fromfile(f,sfh_template.sfh_struct_dtype,Nbin)
print(SFHbins.dtype)
print(SFHbins[np.where(SFHbins['Snapnum']==int(snap))])
nSFHbins_used=SFHbins[np.where(SFHbins['Snapnum']==int(snap))]['Bin'][-1]
print(nSFHbins_used)
dt=SFHbins[np.where(SFHbins['Snapnum']==int(snap))]['dt']
print('dt=',dt)
LookbackTime=SFHbins[np.where(SFHbins['Snapnum']==int(snap))]['LookbackTime']
print('LookbackTime=',LookbackTime)


# In[62]:


# First determine the size of the arrays that we need to hold the data
nGal=0
for iFile in range(firstFile,lastFile+1):
    # The following line closes the file at the end of the loop
    with h5py.File(outputDir+filePrefix+'%i'%iFile+filePostfix,'r') as f:
        nGal+=len(f[snap])
print('nGal=',nGal)

# Declare numpy arrays to hold the data
Type=np.empty(nGal,dtype=np.int32)
DiskMass=np.empty(nGal)
sfh_DiskMass=np.empty([nGal,nSFHbins])

# Now read in the data
iGal=0
for iFile in range(firstFile,lastFile+1):
    # The following line closes the file at the end of the loop
    with h5py.File(outputDir+filePrefix+'%i'%iFile+filePostfix,'r') as f:
        nGalFile=len(f[snap])
        Type[iGal:iGal+nGalFile]=f[snap]['Type']
        DiskMass[iGal:iGal+nGalFile]=f[snap]['DiskMass']
        sfh_DiskMass[iGal:iGal+nGalFile,:]=f[snap]['sfh_DiskMass']
        iGal+=nGalFile
        
# Select galaxies in mass range 3e9-1e10 Msun/h and convert to Msun/h
index=np.where((DiskMass>0.3) & (DiskMass<=1.0))[0]
Type=Type[index]
DiskMass=DiskMass[index]*1e10
sfh_DiskMass=sfh_DiskMass[index,:nSFHbins_used+1]*1e10


# In[66]:



# Sanity check on mass
# Note that SFH stores initial mass and DiskMass stores remaining mass for an evolved stellar population
plt.figure(figsize=[12,8])
plt.semilogx(DiskMass,np.sum(sfh_DiskMass,1)/DiskMass,'.')
plt.xlabel(r'DiskMass/$h^{-1}$M$_\odot$')
plt.ylabel(r'np.sum(sfh_DiskMass,1)/DiskMass')
plt.savefig('figs/sfh_massCheck.png')
plt.show()

# Average star formation history over all galaxies
sfh_DiskMass_mean=np.mean(sfh_DiskMass,0)
# Average star formation rate in each bin
sfr_mean=sfh_DiskMass_mean/dt
# Plot star formation rate as a function of Lookback time
plt.figure(figsize=[12,8])
plt.semilogx(LookbackTime,sfr_mean,'o')
plt.semilogx(LookbackTime,sfr_mean,'-')
plt.xlabel(r'Lookback time / yr')
plt.ylabel(r'Star formation rate / (Msun/yr)')
plt.savefig('figs/sfh_sfr.png')
plt.show()

