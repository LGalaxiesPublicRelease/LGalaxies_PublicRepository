# script to remove galaxies below some mass limit

import numpy as np

massCut=0.1
if model=='Hen15_MRII': massCut = 0.0008

gals=gals[np.where(gals['StellarMass'] >= massCut)]
