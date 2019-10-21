# Script to unpickle (i.e. read from disk) galaxy properties.
# This is the python3 version

import pickle

python2data=False
#pickleFile='data_pkl/Hen15/snap_58.pkl'

fin = open(pickleFile, 'rb')
if python2data:
    gals=pickle.load(fin,encoding='latin')
else:
    gals=pickle.load(fin)
fin.close()
