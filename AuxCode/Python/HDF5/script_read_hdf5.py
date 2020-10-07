# Example python script to read in HDF5 data and produce plot

# Imports
import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('poster')
sns.set_style('whitegrid')

# Parameters
dataDir='/research/astrodata/virgo/FLARES/data/'
filenr=0

# Open HDF5 file in a way that closes at end of loop
with h5py.File(dataDir+'GEAGLE_%2.2i_sp_info.hdf5' %filenr,'r') as f:
    # View contents
    print(list(f.items()))
    print()
    # Open one of the listed groups and view its contents
    g=f['006_z009p000']
    print(list(g.items()))
    print()
    # Likewise for the contents of each of those groups
    Particle=g['Particle']
    print(list(Particle.items()))
    print()
    Subhalo=g['Subhalo']
    print(list(Subhalo.items()))
    print()
    # Extract some pertinent quantities
    # This seems to be a list of subgroup positions and masses
    COP=Subhalo['COP']
    Mstar_30=Subhalo['Mstar_30']
    # This seems to be a list of properties of each star particle
    # position, mass, metallicity and age would seem to be useful
    S_Coordinates=Particle['S_Coordinates'][:]
    S_Mass=Particle['S_Mass']
    S_Age=Particle['S_Age']
    S_Z=Particle['S_Z']

plt.figure(figsize=[12,12])
plt.plot(S_Coordinates[0,:],S_Coordinates[1,:],',')
plt.axis('Equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Stellar positions')
plt.savefig('figs/testplot.png',bbox_inches='tight')
plt.show()


