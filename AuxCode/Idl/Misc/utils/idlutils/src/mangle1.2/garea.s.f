c-----------------------------------------------------------------------
      subroutine garea(area,rp,cm,np,tol,verb,phi,iord,ldegen)
      integer np,verb,iord(2*np)
      logical ldegen
      real*8 area,rp(3,np),cm(np),tol,phi(2,np)
c
c        parameters
      include 'pi.par'
      real*8 TWOPI
      parameter (TWOPI=2.d0*PI)
c        externals
      integer gsegij,gzeroar,ikrand
c        data variables
      real*8 arearnd,big
c     real*8 dphmin
c        local variables
      integer i,icmmin,ik,iseg,j,jm,jml,jmu,jp,jpl,jpu,k,km,kp,l,
     *  nbd,ni,scmi
c     logical warn
      logical whole
      real*8 cmi,cmik,cmmin,cmk,cpsi,darea,dph,ikchk,
     *  p,ph,phm,php,psi,si,sk,xi(3),yi(3)
c *
c * Area of surface of sphere of unit radius bounded by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c * See AJSH notes Multfn C115.
c * Cautions:
c * (1) This subroutine underestimates the area by 2*pi
c *     if 2*pi <= area < 4*pi
c *     and the area is bounded by more than one arc.
c * (2) This subroutine will usually work correctly when there are near
c *     multiple (.ge. 3) intersections of boundaries, but in rare
c *     instances it may fail.  If so, it should flag the failure with
c *     ldegen=.true.  This error condition should NOT be ignored.
c *
c  Input: rp(3,i),i=1,np
c         cm(i),i=1,np
c         np
c         verb
c Output: area
c         ldegen = .true. means there's a problem with multiply
c                  intersecting boundary.
c Input/Output: tol
c Work arrays: phi and iord should be dimensioned at least 2*np
c
c        set azimuthal angle of non-intersection to big
      data big /1.d6/
c        ok if area tests not too far outside [0,max]
      data arearnd /1.d-10/
c        warn about multiple intersection when dph < dphmin
c     data dphmin /1.d-8/
c
c        come here with enlarged tolerance to multiple intersections
C     print *,'--------------------'
  100 continue
c        initialise error flag to no error
      ldegen=.false.
c     warn=.false.
c        zero area
      area=0.d0
c        check for zero area because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        no constraints at all will mean area is whole sphere
      whole=.true.
c        number of arc segments bounding area
      nbd=0
c        error check on evaluation of vertex terms
      ikchk=0.d0
c--------identify boundary segments around each circle i in turn
      do 280 i=1,np
c        cm(i).ge.2 means include whole sphere, which is no constraint
        if (cm(i).ge.2.d0) goto 280
c        there is a constraint, so area is not whole sphere
        whole=.false.
c        scmi * cmi = 1-cos th(i)
        if (cm(i).ge.0.d0) then
          scmi=1
        else
          scmi=-1
        endif
        cmi=abs(cm(i))
c        si = sin th(i)
        si=sqrt(cmi*(2.d0-cmi))
c........construct cartesian axes with z-axis along rp(i)
        call gaxisi(rp(1,i),xi,yi)
c........angles phi about z-axis rp(i) of intersection of i & j circles
        call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,ni,phi)
c        i circle lies outside polygon 
        if (ni.eq.-1) goto 280
c        area of polygon is zero
        if (ni.eq.-2) goto 410 
c........i circle has no intersections
        if (ni.eq.0) then
          nbd=nbd+1
          darea=cm(i)*TWOPI
          area=area+darea
C         print *,'at',i,': full circle area +=',darea,' =',area
c........i circle has intersections
        elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
          call findbot(phi,2*np,iord,ni)
c........contribution to area from each segment of i circle
          jpl=0
c        come here to do another segment
  200     continue
c........is segment edge of polygon?
            iseg=gsegij(rp,cm,np,rp(1,i),scmi,cmi,ni,tol,phi,iord,
     *        jml,jmu,jpl,jpu,jm,jp,km,kp,phm,php,ph,dph)
c        error
            if (iseg.eq.-1) goto 420
c        not an edge
            if (iseg.eq.0) goto 200
c        gone full circle
            if (iseg.eq.2) goto 280
c        warn about near multiple intersection
c           if (dph.lt.dphmin) then
c             print *,
c    *          '*** warning from garea: near multiple intersection at',
c    *          i,': edge',km,kp,' dph=',dph
c             warn=.true.
c           endif
c........segment satisfies conditions
            nbd=nbd+1
c        there's a contribution to area from the segment...
            if (cm(i).lt.0.d0) dph=-dph
            darea=cmi*dph-dph
            area=area+darea
C           print *,'at',i,': edge',km,kp,
C    *       ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')'
C           print *,'dph =',dph,' area +=',darea,' =',area
c        ...and from the end points of the segment
            do 240 l=1,2
c        end point is intersection of i circle with k circle
              if (l.eq.1) then
                k=km
              elseif (l.eq.2) then
                k=kp
              endif
c        only do ik intersection once,
c        but check both ik and ki intersections made it here
c        from segment k to segment i right-handedly through vertex
              if (l.eq.1.and.cm(i).gt.0.d0.or.l.eq.2.and.cm(i).lt.0.d0)
     *          then
                ik=k+np*i
                if (k.lt.i) ik=ik+1
c        from segment i to segment k right-handedly through vertex
              else
                ik=i+np*k
                if (i.lt.k) ik=ik+1
              endif
c        pseudo-random number from ik
              ik=ikrand(ik)
              if (i.gt.k) then
                ikchk=ikchk-ik
                goto 240
              endif
              ikchk=ikchk+ik
              cmk=abs(cm(k))
              sk=sqrt(cmk*(2.d0-cmk))
c        cmik = 1-cos th(ik)
              cmik=((rp(1,i)-rp(1,k))**2+(rp(2,i)-rp(2,k))**2
     *          +(rp(3,i)-rp(3,k))**2)/2.d0
c        cpsi = (cik-ci*ck)/(si*sk)
              cpsi=(cmi+cmk-cmi*cmk-cmik)/(si*sk)
              if ((cm(i).ge.0.d0.and.cm(k).lt.0.d0)
     *          .or.(cm(i).lt.0.d0.and.cm(k).ge.0.d0)) cpsi=-cpsi
c        psi = exterior angle (pi-interior angle) at intersection point
              psi=acos(cpsi)
              area=area-psi
C             print *,'     intersect',i,k,' area +=',-psi,' =',area
  240       continue
c        do another segment
          goto 200
        endif
  280 continue
c--------finish off
c        whole=.true. means area is whole sphere
      if (whole) then
        area=2.d0*TWOPI
        goto 410
c        just one boundary segment
      elseif (nbd.eq.1) then
        if (area.lt.0.d0) area=area+2.d0*TWOPI
        goto 410
      endif
c        add/subtract 2*pi's to area to ensure 0.le.area.lt.2*pi
      i=area/TWOPI
      area=area-i*TWOPI
C     if (i.ne.0) print *,'area +=',i,' * TWOPI =',area
      if (area.lt.0.d0) then
        area=area+TWOPI
C       print *,'area += TWOPI =',area
      endif
c        check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0.d0) then
        ldegen=.true.
c       warn=.true.
C       print *,'*** from garea: at tol =',tol,
C    *    ', ikchk=',ikchk,' should be 0'
c       write (*,'(a3,a20,3a24)')
c    *    ' ','x','y','z','1-c'
c       write (*,'(i3,4g24.16)')
c    *    (j,(rp(i,j),i=1,3),
c    *    cm(j),j=1,np)
c        retry with enlarged tolerance
        if (tol.le.0.d0) then
          tol=1.d-15
        else
          tol=tol*2.d0
        endif
        goto 100
      elseif (tol.gt.0.d0) then
C       print *,'... from garea: success at tol =',tol
      endif
c        check area is positive
      if (area.lt.0.d0) then
c        negative area should `never' happen
        ldegen=.true.
c       warn=.true.
        if (verb.ge.1) then
          print *,'*** from garea: area=',area,' should be .ge. 0'
        endif
      endif
c        check area does not exceed area within any one circle
      cmmin=2.d0
      do i=1,np
        if (cm(i).ge.0.d0) then
          if (cm(i).lt.cmmin) then
            cmmin=cm(i)
            icmmin=i
          endif
        elseif (cm(i).lt.0.d0) then
          if (2.d0+cm(i).lt.cmmin) then
            cmmin=2.d0+cm(i)
            icmmin=i
          endif
        endif
      enddo
      darea=TWOPI*cmmin
      p=area/darea
C     print *,'area/area(',icmmin,')=',area,' /',darea,' =',p
      if (p.gt.1.d0) then
c        check if discrepancy is from numerical roundoff
        if (abs(area-TWOPI).le.arearnd) then
          area=0.d0
          goto 410
        endif
        if (area.le.darea+arearnd) then
          area=darea
          goto 410
        endif
        ldegen=.true.
c       warn=.true.
        if (verb.ge.1) then
          print *,'*** from garea: area/area(',icmmin,')=',
     *      area,' /',darea,' =',p,'  should be .le. 1'
        endif
        goto 420
      endif
c     if (warn) then
c       write (*,'(a3,a20,3a24)')
c    *    ' ','x','y','z','1-c'
c       write (*,'(i3,4g24.16)')
c    *    (j,(rp(i,j),i=1,3),
c    *    sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
c    *    cm(j),j=1,np)
c     endif
  410 continue
C     print *,'final area =',area
C     print *,'....................'
      return
c
  420 continue
      ldegen=.true.
      if (verb.ge.1) then
        print *,'*** from garea: total failure'
        write (*,'(a3,a20,3a24)')
     *    ' ','x','y','z','1-c'
        write (*,'(i3,4g24.16)')
     *    (j,(rp(i,j),i=1,3),
c    *    sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *    cm(j),j=1,np)
      endif
C     print *,'....................'
      return
c
      end
c
