
# coding: utf-8

# # Python Plots for LGalaxies

# ## Import Libraries and Read Catalogs

# <p>Use functions read_snap or read_tree to read catalogs. These are both defined in procedures.py. In case of read_snap, SnapshotList will be returned containing the list of snapshots read (usefull to later select galaxies in a given redshift).<p>

# In[83]:

import numpy as np
get_ipython().magic('matplotlib inline')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from importlib import reload
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *

from LGalaxies_Henriques2015a_struct import LGalaxiesStruct
from LGalaxies_Henriques2015a_struct import PropertiesToRead


FirstFile = 40
LastFile = 79

Volume_MR = (BoxSize_MR**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 
Volume_MRII = (BoxSize_MRII**3.0) * (LastFile - FirstFile + 1) / MaxTreeFiles 

(G_MR, SnapshotList) = read_snap(DirName_MR,FirstFile,LastFile,
                 PropertiesToRead,LGalaxiesStruct,RedshiftsToRead,RedshiftList)
    
#print (np.log10(gal['StellarMass'][1:5]*1.e10))
#help(gal)


# ## Plots

# In[126]:

plt.rcParams.update({'font.size': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.linewidth': 2, 
                     'xtick.major.size': 6, 'xtick.major.width': 1.5, 
                     'ytick.major.size': 6, 'ytick.major.width': 1.5, 
                     'xtick.minor.size': 3, 'xtick.minor.width': 1.,                   
                     'ytick.minor.size': 3, 'ytick.minor.width': 1.})


# In[149]:

xmin=7.0
xmax=12.5
ymin=-6.5
ymax=0.5
bin=0.1


plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                     'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

fig = plt.figure(figsize=(9,9))
grid = gridspec.GridSpec(2, 2)
grid.update(wspace=0.0, hspace=0.0)

for ii in range(0,4):
    
    subplot=plt.subplot(grid[ii])

    subplot.set_ylim([ymin, ymax]), subplot.set_xlim([xmin, xmax])
    if ii==2 or ii == 3:
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
    else:
        xlab=''
    if ii==0 or ii == 2:
        ylab='$log_{10}(\phi [h^3 Mpc^{-3} log_{10}(M^{-1})])$'
    else:
        ylab=''      
    subplot.set_xlabel(xlab, fontsize=16)
    subplot.set_ylabel(ylab, fontsize=16)
    
   
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
    subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
    if ii==1 or ii == 3:
        plt.tick_params(axis='y', which='both', left='on', labelleft='off')
    
    #OBSERVATIONS
    if ii == 0 :
        char_redshift="%0.0f" % RedshiftList[ii]
    else :  #avoid using z=0.4
        char_redshift="%0.0f" % RedshiftList[ii+1]
        
    file = Datadir + '/ObsConstraints/StellarMassFunction_z'+char_redshift+'.00.txt'
    obs = Table.read(file, format='ascii')
    subplot.errorbar(obs['col1'], np.log10(obs['col3']),yerr=obs['col4'], 
             fmt='o', markersize=5, ecolor='blue', color='blue')
    #sub = plt.subplot(111)

    #MODEL
    if ii == 0 :
        sel=G_MR['SnapNum']==SnapshotList[ii]
    else :  #avoid using z=0.4    
        sel=G_MR['SnapNum']==SnapshotList[ii+1]
    G0_MR=G_MR[sel]
    StellarMass=np.log10(G0_MR['StellarMass']*1.e10*Hubble_h)

    bin_arr=np.arange(7.0,12.0+bin,bin)
    hist=np.histogram(StellarMass, bins=bin_arr, range=(7.0,12.0))
    subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
             color='red', linewidth=2)


    #LABELS
    if ii==0:
        subplot.text(7.3, 0.0, 'Observations used in MCMC')
        subplot.errorbar(7.2, 0.16, yerr=0.15, fmt='o', markersize=5, color='blue')
#endfor


plt.tight_layout()
plt.savefig('./fig/plots.pdf')


# In[134]:




# In[ ]:



