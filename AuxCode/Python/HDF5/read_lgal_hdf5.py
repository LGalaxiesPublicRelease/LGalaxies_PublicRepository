import h5py 
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#========================================================================
#
#   Sample code to read in the L-galaxies hdf5 data 
#
#
#=======================================================================

start_time=time.time()

datadir='../../../output/'
snap='17'
filenr=40

with h5py.File(datadir+'SA_output_%i.h5' %filenr,'r') as f:
    
    #Get the data from the snapshot
    data=f[snap]
    #Get an array of the data returned
    Mvir=data['Mvir']

    # want to only plot galaxies which have a central viral mass than 50
    selection=[f[snap]['Mvir']>20.0] #10^10 msun
    
    # Get the coordinates of the data that follow the above selection
    x=f[snap]['Pos'][:,0][selection]
    y=f[snap]['Pos'][:,1][selection]  
    z=f[snap]['Pos'][:,2][selection]
    SfrRings=f[snap]['SfrRings']
    Mag=f[snap]['Mag']

    #Output the possible data labels where the first 3 are: table class, table version, table name
    print(list(f[snap].attrs.values()))


