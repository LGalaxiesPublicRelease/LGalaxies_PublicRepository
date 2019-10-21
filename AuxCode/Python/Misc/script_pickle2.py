# Script to pickle (i.e. dump to disk) galaxy properties
# The read back in python3, use pickle.load(fin,encoding='latin')
# Requires pickleFile to contin the name of the output file,
# and gals to hold the data to be pickled.

import cPickle

pickleFile=odatadir+'/snap_'+str(snapshot)+'.pkl'
fout = open(pickleFile, 'wb')
cPickle.dump(gals,fout,cPickle.HIGHEST_PROTOCOL)
fout.close()
