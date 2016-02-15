c-----------------------------------------------------------------------
      subroutine gphi(angle,rp,cm,np,rpi,cmi,tol,phi,iord)
      integer np,iord(2*np)
      real*8 angle,rp(3,np),cm(np),rpi(3),cmi,tol,phi(2,np)
c
c        parameters
      include 'pi.par'
      real*8 TWOPI
      parameter (TWOPI=2.d0*PI)
c        externals
      integer gsegij,gzeroar
c        data variables
      real*8 big
c        local variables
      integer i,iseg,j,jm,jml,jmu,jp,jpl,jpu,km,kp,ni,scmi
      logical warn
      real*8 dph,p,ph,phm,php,xi(3),yi(3)
c *
c * Angle along circle about unit direction rpi satisfying
c *    1 - r.rpi = cmi
c * and bounded by
c *    1 - r.rp(j) <= cm(j)  (if cm(j).ge.0)
c *    1 - r.rp(j) > -cm(j)  (if cm(j).lt.0)
c * for j=1,np where rp(j) are unit directions.
c *
c  Input: rp(3,j),j=1,np
c         cm(j),j=1,np
c         np
c         rpi(3)
c         cmi
c         tol
c Output: angle
c Work arrays: phi and iord should be dimensioned at least 2*np
c
      data big /1.d6/
c
c        check for null circle
      if (cmi.lt.0.d0) goto 410
      if (cmi.gt.2.d0) goto 410
c        check for zero angle because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
      angle=0.d0
      scmi=1
c........construct cartesian axes with z-axis along rpi
      call gaxisi(rpi,xi,yi)
c........angles phi about z-axis rp(i) of intersection of i & j circles
      call gphij(rp,cm,np,0,rpi,scmi,cmi,xi,yi,big,ni,phi)
c        i circle lies outside polygon 
      if (ni.le.-1) goto 410
c        i circle lies outside (or at edge of) polygon
c........angles phi about z-axis rp(i) of intersection of i & j circles
      call gphij(rp,cm,np,0,rpi,scmi,cmi,xi,yi,big,ni,phi)
c        i circle lies outside (or at edge of) polygon
      if (ni.le.-1) goto 410
c........order angles around circle
c        circle has no intersections
      if (ni.eq.0) then
        angle=TWOPI
c        circle has intersections
      elseif (ni.gt.0) then
c        find ordering of intersection angles around circle
        call findbot(phi,2*np,iord,ni)
c........vertices around i circle
        jpl=0
c        come here to do another segment
  200   continue
c........is segment edge of polygon?
          iseg=gsegij(rp,cm,np,rpi,scmi,cmi,ni,tol,phi,iord,
     *      jml,jmu,jpl,jpu,jm,jp,km,kp,phm,php,ph,dph)
c        error
          if (iseg.eq.-1) goto 420
c        not an edge
          if (iseg.eq.0) goto 200
c        gone full circle
          if (iseg.eq.2) goto 270
c........segment satisfies conditions
c         print *,'segment',km,kp,' dph/(2*pi)=',real(dph/TWOPI)
          angle=angle+dph
c        do another segment
        goto 200
      endif
  270 continue
c........finish off
c        check angle is between 0 and 2*pi
      p=angle/TWOPI
c     print *,real(rpi(1)),real(rpi(2)),real(rpi(3)),
c    *  'angle/(2*pi)=',p
      if (real(p).lt.0.) then
        print *,'*** from gphi: angle/(2*pi)=',p,
     *    '  should be .ge. 0'
        warn=.true.
      elseif (real(p).gt.1.) then
        print *,'*** from gphi: angle/(2*pi)=',p,
     *    '  should be .le. 1'
        warn=.true.
      endif
      if (warn) then
        write (*,'(a3,a10,4a14)')
     *    ' ','x','y','z','r','1-c'
        write (*,'(i3,5g14.6)')
     *    (j,(rp(i,j),i=1,3),sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *    cm(j),j=1,np),
     *    (rpi(i),i=1,3),sqrt(rpi(1)**2+rpi(2)**2+rpi(3)**2),cmi
      endif
      return
c
c        zero angle
  410 angle=0.d0
      return
c
  420 print *,'*** from gvphi: total failure'
      return
c
      end
c
