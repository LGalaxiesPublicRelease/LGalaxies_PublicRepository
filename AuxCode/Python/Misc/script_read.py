#=========================================================================
#
#  Script to read in L-galaxies snapshot data
#
#  To force a re-read of the data do Gal=None
#
#-------------------------------------------------------------------------

# Imports

import sys

# Path to data
datadir = '/mnt/lustre/scratch/virgo/SAM_output/Hen14_sfh2'
sys.path.insert(0,datadir)

# Template structure for L-Galaxies data
import snap_template   # structure temple for data
import read_lgal       # function to read in data

#-------------------------------------------------------------------------

# Parameters

# Snaplist file
snaplist_file = datadir+'/MRPlancksnaplist.txt'

# Define what snapshot you want
snapshot = 58

# Define which files you want to read in
firstfile = 0
lastfile = 511

# Define what properties you want to read in
props = snap_template.properties_used
props['Type'] = True
props['Mvir'] = True
props['DiskMass'] = True
props['BulgeMass'] = True

#-------------------------------------------------------------------------

# Working body of the program

# Read in redshift of snapshot and create file prefix
f = open(snaplist_file)
lines = f.readlines()
f.close()
for this_line in lines:
    words = this_line.split()
    #print words[0],words[2]
    if words[0]==str(snapshot):
        file_prefix = "SA_z"+words[2]

# Read in galaxy output
(nTrees,nHalos,nTreeHalos,gals) = \
    read_lgal.read_snap(datadir,file_prefix,firstfile,lastfile,\
                            props,snap_template.struct_dtype)
