import numpy as np

DZ=0.1
ZMAX=10.0001

def read_snapshot_table(file):
    """Function to read snapshot table and return disctionaries of
       snapshot versus redshift and vice versa"""
    snap,a,z,t,tyr=np.loadtxt(file,skiprows=1,unpack=True,
                              dtype=[('Snap','i4'),('a','f4'),('z','f4'),('t','f4'),('tyr','f4')])
    nsnap=len(snap)
    snap_to_z={snap[i]:z[i] for i in np.arange(nsnap)}
    z_to_snap={}
    for zz in np.arange(0,ZMAX,DZ):
        # Lots of rounding error, so this is a fudge to fix it.
        ilow=np.where(z > zz)[0][-1] # ilow is index of snaps with higher redshift
        if ilow==nsnap-1: ilow=ilow-1
        if z[ilow]+z[ilow+1]>2*zz:
            z_to_snap.update([(str(DZ*(zz/DZ)),snap[ilow])])
        else:
            z_to_snap.update([(str(DZ*(zz/DZ)),snap[ilow+1])])
    return (snap_to_z,z_to_snap)

