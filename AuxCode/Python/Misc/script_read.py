#=========================================================================
#
#  Script to read in L-galaxies snapshot data
#
#  To force a re-read of the data do gals=None
#
#-------------------------------------------------------------------------

# Imports

import os
import sys

# Template structure for L-Galaxies data
import read_lgal       # function to read in data

#-------------------------------------------------------------------------

# Parameters

# Decide what data set you want
# Allow for calling with model set from metascript
try:
    model
except:
    #model='Hen15'
    #model='Hen15_spin'
    model='Guo11'
    #model='Hen15_MRII'

# Decide what redshift you want
# Allow for calling with redshift set from metascript
try:
    redshift
except:
    redshift=0.15

# Define which files you want to read in
# 0-511 gives all the files (but takes a long time)
try:
    firstfile
except:
    firstfile = 0
    lastfile = 511

# Path to input data
if model=='Guo11': idatadir = 'data_raw/Guo11'
if model=='Guo11_spin': idatadir = 'data_raw/Guo11'
if model=='Hen15': idatadir = 'data_raw/Hen15'
if model=='Hen15_spin': idatadir = 'data_raw/Hen15'
if model=='new': idatadir = 'data_raw/new'
if model=='Rob': idatadir = 'data_raw/Rob'

# Path for output data
if model=='Guo11': odatadir = 'data_pkl/Guo11'
if model=='Guo11_spin': odatadir = 'data_pkl/Guo11'
if model=='Hen15': odatadir = 'data_pkl/Hen15'
if model=='Hen15_spin': odatadir = 'data_pkl/Hen15'
if model=='new': odatadir = 'data_pkl/new'
if model=='Rob': odatadir = 'data_pkl/Rob'

# import template
sys.path.insert(0,idatadir)
import snap_template

# Define what properties you want to read in
props = snap_template.properties_used
props['Type'] = True
props['HaloIndex'] = True
props['Rvir'] = True
props['Mvir'] = True
props['CentralMvir'] = True
props['DistanceToCentralGal'] = True
props['Pos'] = True
props['Vel'] = True
props['Vmax'] = True
props['InfallVmax'] = True
props['InfallVmaxPeak'] = True
if 'spin' in model:
#    props['BulgeSpin'] = True
#    props['BulgeSpinMax'] = True
    props['DiskSpin'] = True
#    props['DiskSpinMax'] = True
    props['ColdGasSpin'] = True
props['BulgeMass'] = True
props['DiskMass'] = True
props['StellarMass'] = True
props['ColdGas'] = True
props['HotGas'] = True
props['EjectedMass'] = True
props['MetalsDiskMass'] = True
#props['MetalsStellarMass'] = True
props['MetalsBulgeMass'] = True
props['MetalsColdGas'] = True
props['MetalsHotGas'] = True
props['MetalsEjectedMass'] = True
props['BulgeSize'] = True
if 'spin' in model:
    props['DiskRadius'] = True
    props['ColdGasRadius'] = True
props['Sfr'] = True
props['BlackHoleMass'] = True
props['QuasarAccretionRate'] = True
props['RadioAccretionRate'] = True
if model == 'Rob': props['Mag'] = True
if model == 'Guo11' or model == 'Hen15': props['ObsMag'] = True

#-------------------------------------------------------------------------

# Working body of the program

# Matching between redshift and snapshot provided by dictionary
# In future create this dictionary from appropriate files
snapz_Hen15={0:58,0.15:53,0.25:50,0.5:45,1:38,1.5:34,2:30,2.5:28,3:25,4:22,5:19,6:17}
snapz_Guo10={0:63,0.15:57,0.25:54,0.5:48,0.75:52,1:41,1.5:36,2:32,3:27,4:24,5:21,6:18}
snapz_Hen15_MRII={0:62,0.25:54,0.5:49,1:42,1.5:38,2:34,3:29,4:26,5:23,6:21}
if model=='Hen15': snapshot = snapz_Hen15[redshift]
if model=='Hen15_spin': snapshot = snapz_Hen15[redshift]
if model=='new': snapshot = snapz_Hen15[redshift]
if model=='Guo11': snapshot = snapz_Guo10[redshift]
if model=='Guo11_spin': snapshot = snapz_Guo10[redshift]
if model=='Hen15_MRII': snapshot = snapz_Hen15_MRII[redshift]
if model=='Rob': snapshot = snapz_Hen15[redshift]

# Snaplist file
for word in os.listdir(idatadir):
    if 'snaplist' in word:
        snaplist_file=idatadir+'/'+word

# Read in redshift of snapshot and create file prefix
with open(snaplist_file) as f:
    lines = f.readlines()
    for this_line in lines:
        words = this_line.split()
        #print words[0],words[2]
        if words[0]==str(snapshot):
            if model=='Rob':
                file_prefix = "SA1_z"+words[2]
            else:
                file_prefix = "SA_z"+words[2]

# Read in galaxy output
(nTrees,nHalos,nTreeHalos,gals) = \
    read_lgal.read_snap(idatadir,file_prefix,firstfile,lastfile,\
                            props,snap_template.struct_dtype)
