c-----------------------------------------------------------------------
c © A J S Hamilton 2001
c-----------------------------------------------------------------------
      real*8 function felp(epoch)
      real*8 epoch
c
c        parameters
      include 'frames.par'
c        local (automatic) variables
      real*8 t
c *
c * Ecliptic latitude of NCP = Dec of ecliptic NP
c * as a function of epoch (e.g. 1950, 2000).
c *
c        RA & Dec epoch in centuries since 1900
      t=(epoch-1900.d0)/100.d0
c        ecliptic latitude of NCP = Dec of ecliptic NP
      felp=90.d0-(E1+t*(E2+t*(E3+t*E4)))
      return
      end
c
