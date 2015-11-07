# script to pickle (i.e. dump to disk) galaxy properties

import cPickle

fin = open('gal.pkl', 'rb')
gals=cPickle.load(fin)
fin.close()
