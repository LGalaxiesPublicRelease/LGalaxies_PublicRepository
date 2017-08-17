# Script to compare 2 L-Galaxies outputs

import numpy as np
import pickle

#-----------------------------------------------------

datadir1="/Users/petert/lgalaxies/Development_Branch_Apr21/output/Guo11_old/"
datadir2="/Users/petert/lgalaxies/Development_Branch_Apr21/output/Guo11/"
datafile="snap_57.pkl"

#-----------------------------------------------------

with open(datadir1+datafile, 'rb') as fin:
    gals1=pickle.load(fin)
with open(datadir2+datafile, 'rb') as fin:
    gals2=pickle.load(fin)

# Number of galaxies
if len(gals1) != len(gals2): 
    print('len1=',len(gals1))
    print('len2=',len(gals2))
    raise ValueError("Different galaxy numbers")
else:
    print('Same number of galaxies')

# Galaxy types
index=np.where(gals1['Type'] != gals2['Type'])
if len(index[0]) != 0: 
    raise ValueError("Different Type")
else:
    print('Types agree')

# Halo indices
index=np.where(gals1['HaloIndex'] != gals2['HaloIndex'])
if len(index[0]) != 0: 
    raise ValueError("Different HaloIndex")
else:
    print('HaloIndices agree')

ratio=gals1['StellarMass']/gals2['StellarMass']
print('min,max[mass ratio]={:0.3f},{:0.3f}'.format(np.min(ratio),np.max(ratio)))
if (np.min(ratio)<0.999 or np.max(ratio)>1.001):
    raise ValueError("Masses differ")
else:
    print("Masses agree")

print("Tests successful")
