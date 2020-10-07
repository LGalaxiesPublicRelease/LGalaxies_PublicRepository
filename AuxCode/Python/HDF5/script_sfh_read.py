# Script to read in SFHbin file comtaining information about SFH bins

import numpy as np

# Directory where SFHbin data is stored
idatadir = '/mnt/lustre/scratch/virgo/SAM_output/Hen14_sfh4/'
#idatadir='/mnt/lustre/scratch/petert/L-Galaxies/output/'

# import template
#sys.path.insert(0,idatadir)
import sfh_template

# Open data file for binary reading
f = open(idatadir+'SFH_Bins','rb')
Nbin = np.fromfile(f,np.int32,1)[0]
# Position file pointer at end of header
f.seek(sfh_template.sfh_struct_dtype.itemsize,0)

# Create record array to hold data and read it in
#SFHbins=np.array(Nbin,dtype=sfh_struct_dtype)
SFHbins=np.fromfile(f,sfh_template.sfh_struct_dtype,Nbin)

# Tidy up
f.close()


