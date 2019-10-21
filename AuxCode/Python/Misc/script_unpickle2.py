# Script to unpickle (i.e. read from disk) galaxy properties.
# This is the python2 version

import cPickle

datafile='data_pkl/Hen15/snap_58.pkl'

fin = open(datafile, 'rb')
gals=cPickle.load(fin)
fin.close()
