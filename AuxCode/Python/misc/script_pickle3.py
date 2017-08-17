# Script to pickle (i.e. dump to disk) galaxy properties.
# Don't use this script if you want to read the data back in using
# python2 as the string encoding in numpy dtypes is incompatible.

import pickle

pickleFile=odatadir+'/snap_'+str(snapshot)+'.pkl'
fout = open(pickleFile, 'wb')
pickle.dump(gals,fout)
# Make python2 readable (but doesn't seem to work)
# pickle.dump(gals,fout,protocol=2,fix_imports=True)   
fout.close()
