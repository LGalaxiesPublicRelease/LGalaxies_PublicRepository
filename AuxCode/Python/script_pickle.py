# script to pickle (i.e. dump to disk) galaxy properties

import cPickle

fout = open('gal.pkl', 'wb')
cPickle.dump(gals,fout,cPickle.HIGHEST_PROTOCOL)
fout.close()
