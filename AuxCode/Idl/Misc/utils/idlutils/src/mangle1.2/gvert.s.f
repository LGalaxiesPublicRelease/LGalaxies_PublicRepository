c-----------------------------------------------------------------------
      subroutine gvert(ve,angle,ipv,ev,nvmax,nv,nve,nev,nev0,
     *  rp,cm,np,vcirc,tol,phi,iord,wk,iwk,ldegen)
      integer nvmax,ipv(nvmax),ev(nvmax),nv,nve,nev,nev0,
     *  np,vcirc,iord(2*np),iwk(nvmax,3)
      logical ldegen
      real*8 ve(3,nve,nvmax),angle(nvmax),
     *  rp(3,np),cm(np),tol,phi(2,np),wk(nvmax)
c
c        parameters
      include 'pi.par'
      real*8 TWOPI
      parameter (TWOPI=2.d0*PI)
c        externals
      integer gsegij,gzeroar
c        data variables
      real*8 big
c     real*8 dphmin
c        local (automatic) variables
      integer i,ii,ik,iseg,iv,ive,j,jm,jml,jmu,jp,jpl,jpu,km,kp,ni,scmi
      logical warn
      real*8 cmi,dph,ikchk,ph,phm,php,si,xi(3),xv,yi(3),yv,zv
c *
c * Vertices, plus points on edges, of the polygon defined by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c *
c * Since each pair of circles can intersect at at most 2 points,
c * an upper limit on the number of possible vertices is np*(np+1),
c * twice the number of distinct pairs of np objects.
c * I can find configurations which exhaust this bound for np <= 3.
c * For np >=4, tighter upper limits undoubtedly exist.
c *
c * The vertices are ordered right-handedly about the polygon.
c * The array ev records the end index of each connected
c * sequences of vertices.
c * Intersecting boundaries come 1st, then non-intersecting boundaries.
c * If the bounded region is simply-connected, the usual case,
c * then ev(1) = nv on output, and ev(i) = 0 for i >= 2.
c * If the bounded region is not simply-connected,
c * then the number of non-zero elements of ev will equal the number
c * of distinct connected boundaries of the region.
c *
c  Input: nvmax
c         rp(3,i),i=1,np
c         cm(i),i=1,np
c         np
c         nve = desired number of points per edge
c             = 1 for vertices only
c             = 2 for vertices + midpoints of edges.
c         vcirc = 1 to return vertices and edge points also for bounding
c                 circles which do not intersect any other circle;
c               = 0 otherwise.
c Output: ve(3,ive,iv),ive=1,nve,iv=1,nv are points uniformly spaced
c                   along the edges of the polygon,
c		    starting with vertices in ve(3,1,iv).
c         angle(iv),i=1,nv are azimuthal lengths of edges of polygon.
c         ipv(iv),iv=1,nv is the index of the iv'th edge, that is,
c                 points ve(*,ive,iv) lie on circle number ipv.
c         ev(i),i=1,nev = end index of each connected sequence of edges;
c                 the number of non-zero elements of ev is the number
c                 nev of connected boundaries of the region,
c                 and the last nev0 of these are non-intersecting.
c         nv = number of points;
c              if this exceeds nvmax, then you should call gvert again
c              with a larger nvmax.
c         nev = number of connected sequences of vertices.
c         nev0 = number of non-intersecting boundary circles;
c                non-intersecting boundaries come last in ve,
c                and the last nev0 entries of ev(i) refer
c                to non-intersecting boundaries.
c         ldegen = .true. means there's a problem with multiply
c                  intersecting boundary.
c Input/Output: tol
c Work arrays: phi and iord should be dimensioned at least 2*np.
c              iwk should be dimensioned at least 3*nvmax.
c              wk should be dimensioned at least nvmax.
c
c     data dphmin /1.d-8/
      data big /1.d6/
c
c        come here with enlarged tolerance
  100 continue
c        initialise error flag to no error
      ldegen=.false.
      warn=.false.
c        zero number of vertices
      nv=0
      nev=0
      nev0=0
c        initialise ev to zero
      do iv=1,nvmax
        ev(iv)=0
      enddo
c        check for zero area because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        error check on evaluation of vertex terms
      ikchk=0.d0
c        initialise iwk to inadmissible values
      do iv=1,nvmax
        iwk(iv,2)=-1
      enddo
      do iv=1,nvmax
        iwk(iv,3)=-1
      enddo
c--------identify boundary segments around each circle i in turn
      do 280 i=1,np
c        cm(i).ge.2 means include whole sphere, which is no constraint
        if (cm(i).ge.2.d0) goto 280
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
c        want vertices and midpoints for circles without intersections
          if (vcirc.eq.1) then
c        introduce pretend circle 0
            do j=1,2
              if (j.eq.1) then
                ii=0
                km=i
                kp=i
                phm=0.d0
                php=PI
              elseif (j.eq.2) then
                ii=i
                km=0
                kp=0
                phm=PI
                php=TWOPI
              endif
              dph=php-phm
              nv=nv+1
              if (nv.le.nvmax) then
c        edge index of this vertex
                ipv(nv)=i
c        azimuthal length of edge
                angle(nv)=dph
c        edge points
                zv=1.d0-cmi
                do ive=1,nve
                  if (cm(i).ge.0.d0) then
                    ph=phm+dph*dble(ive-1)/dble(nve)
                  else
                    ph=php-dph*dble(ive-1)/dble(nve)
                  endif
                  xv=si*cos(ph)
                  yv=si*sin(ph)
                  ve(1,ive,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
                  ve(2,ive,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
                  ve(3,ive,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
                enddo
              endif
c        record endpoints of edge
              call gvtrail(scmi,np,ii,km,kp,iwk,iwk(1,2),nvmax,nv,
     *          ik,ikchk)
            enddo
            nev0=nev0+1
          endif
c........i circle has intersections
        elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
          call findbot(phi,2*np,iord,ni)
c........vertices around i circle
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
c........segment is edge
c        warn about near multiple intersection
c           if (dph.lt.dphmin) then
c             print *,'*** warning from gvert: near multiple intersectio
c    *n at',i,': edge',km,kp,' dph=',real(dph)
c             warn=.true.
c           endif
c           print *,'at',i,': edge',km,kp,
c    *        ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')',
c    *        ' ph=',real(ph),' dph=',real(dph)
            nv=nv+1
            if (nv.le.nvmax) then
c        edge index of this vertex
              ipv(nv)=i
c        azimuthal length of edge
              angle(nv)=dph
c        edge points
              zv=1.d0-cmi
              do ive=1,nve
                if (cm(i).ge.0.d0) then
                  ph=phm+dph*dble(ive-1)/dble(nve)
                else
                  ph=php-dph*dble(ive-1)/dble(nve)
                endif
                xv=si*cos(ph)
                yv=si*sin(ph)
                ve(1,ive,nv)=zv*rp(1,i)+xv*xi(1)+yv*yi(1)
                ve(2,ive,nv)=zv*rp(2,i)+xv*xi(2)+yv*yi(2)
                ve(3,ive,nv)=zv*rp(3,i)+xv*xi(3)+yv*yi(3)
              enddo
            endif
c        record endpoints of edge
            call gvtrail(scmi,np,i,km,kp,iwk,iwk(1,2),nvmax,nv,ik,ikchk)
c        do another segment
          goto 200
        endif
  280 continue
c--------order vertices right-handedly about region
      if (nv.gt.0.and.nv.le.nvmax) then
        call gvord(np,nv,iwk,iwk(1,2),iwk(1,3),ev,nev)
        call vpermi(ipv,nv,iwk(1,2),iwk(1,3))
        call vpermdd(ve,3*nve,nv,iwk(1,2),wk)
        call vpermd(angle,nv,iwk(1,2),wk)
      endif
c--------finish off
c        check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0.d0) then
c       print *,'*** from gvert: at tol =',tol,
c    *    ', ikchk=',ikchk,' should be 0'
        ldegen=.true.
        warn=.true.
c        retry with enlarged tolerance
        if (tol.le.0.d0) then
          tol=1.d-15
        else
          tol=tol*2.d0
        endif
        goto 100
      endif
      if (warn) then
        write (*,'(a3,a20,4a24)')
     *    ' ','x','y','z','r','1-c'
        write (*,'(i3,5g24.16)')
     *    (j,(rp(i,j),i=1,3),sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *    cm(j),j=1,np)
      endif
c     if (nv.gt.nvmax) then
c       print *,'*** from gvert: number of vertices =',nv,
c    *    ' exceeds maximum',nvmax
c     endif
  410 continue
      return
c
  420 print *,'*** from gvert: total failure'
      ldegen=.true.
      return
c
      end
c
