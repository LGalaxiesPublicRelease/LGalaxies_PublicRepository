import sys

# Template structure for L-Galaxies data
import read_lgal       # function to read in data

#-------------------------------------------------------------------------

file_prefix='SA_z0.00'
firstfile = 40
lastfile = 40

# import template
idatadir = '../../../output/Hen15'
sys.path.insert(0,idatadir)
import snap_template

# Define what properties you want to read in
props = snap_template.properties_used
props['Type'] = True
props['StellarMass'] = True
props['MetalsStellarMass'] = True

#-------------------------------------------------------------------------

# Read in galaxy outputs
idatadir = '../../../output/Hen15_MinGalOutputMass'
(nTrees2,nHalos2,nTreeHalos2,gals2) = \
    read_lgal.read_snap(idatadir,file_prefix,firstfile,lastfile,\
    props,snap_template.struct_dtype)
idatadir = '../../../output/Hen15'
(nTree1,nHalo1,nTreeHalos1,gals1) = \
    read_lgal.read_snap(idatadir,file_prefix,firstfile,lastfile,\
    props,snap_template.struct_dtype)

#--------------------------------------------------------------------------

StellarMass1=gals1['StellarMass']*1e10
MetalsStellarMass1=gals1['MetalsStellarMass']*1e10
StellarMass2=gals2['StellarMass']*1e10
MetalsStellarMass2=gals2['MetalsStellarMass']*1e10

close()
semilogx(StellarMass1,MetalsStellarMass1/StellarMass1,'b,',label='All masses')
semilogx(StellarMass2,MetalsStellarMass2/StellarMass2,'r,',label='With min mass cut')
legend(loc=2)
show()
savefig('figs/test_minmass.png')
