c-----------------------------------------------------------------------
      integer function gzeroar(cm,np)
      integer np
      real*8 cm(np)
c
c        local (automatic) variables
      integer i
c *
c * Check for zero area because one circle is null.
c *
c  Input: cm
c         np
c Return value: 0 if area is zero because one circle is null,
c               1 otherwise
c
      do i=1,np
        if (cm(i).eq.0.d0) goto 200
        if (cm(i).le.-2.d0) goto 200
      enddo
      gzeroar=1
      return
c
  200 gzeroar=0
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gaxisi(rp,xi,yi)
      real*8 rp(3),xi(3),yi(3)
c
c        local (automatic) variables
      real*8 sx
c *
c * Cartesian axes with z-axis along rp.
c *
c  Input: rp
c Output: xi, yi forming right-handed orthonormal system with rp
c
      sx=rp(1)**2+rp(3)**2
      if (sx.gt..5d0) then
        sx=sqrt(sx)
c        xi in direction y x rp (= x direction if rp is along z)
        xi(1)=rp(3)/sx
        xi(2)=0.d0
        xi(3)=-rp(1)/sx
      else
        sx=sqrt(rp(1)**2+rp(2)**2)
c        xi in direction rp x z
        xi(1)=rp(2)/sx
        xi(2)=-rp(1)/sx
        xi(3)=0.d0
      endif
c        yi in direction rp x xi (= y direction if rp is along z)
      yi(1)=xi(3)*rp(2)-xi(2)*rp(3)
      yi(2)=xi(1)*rp(3)-xi(3)*rp(1)
      yi(3)=xi(2)*rp(1)-xi(1)*rp(2)
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gaxisii(rpi,rp,xi,yi)
      real*8 rpi(3),rp(3),xi(3),yi(3)
c
c        local (automatic) variables
      real*8 ri,sik
c *
c * Cartesian axes with z-axis along rp, x-axis towards rpi.
c *
c  Input: rpi
c Output: xi, yi forming right-handed orthonormal system with rpi
c
c        set yi in direction rp x rpi
      yi(1)=rpi(3)*rp(2)-rpi(2)*rp(3)
      yi(2)=rpi(1)*rp(3)-rpi(3)*rp(1)
      yi(3)=rpi(2)*rp(1)-rpi(1)*rp(2)
c        project yi orthogonal to rp, to avoid problems near yi=0
      ri=yi(1)*rp(1)+yi(2)*rp(2)+yi(3)*rp(3)
      yi(1)=yi(1)-ri*rp(1)
      yi(2)=yi(2)-ri*rp(2)
      yi(3)=yi(3)-ri*rp(3)
      sik=yi(1)**2+yi(2)**2+yi(3)**2
      if (sik.gt.0.d0) then
c        sik = sin th(ik)
        sik=sqrt(sik)
        yi(1)=yi(1)/sik
        yi(2)=yi(2)/sik
        yi(3)=yi(3)/sik
c        rpi is same/opposite direction to rp: set yi along z x rp
      else
        ri=sqrt(rp(1)**2+rp(2)**2)
        if (ri.gt.0.d0) then
          yi(1)=-rp(2)/ri
          yi(2)=rp(1)/ri
          yi(3)=0.d0
c        if rp is also along z-axis, set yi along y-axis
        else
          yi(1)=0.d0
          yi(2)=1.d0
          yi(3)=0.d0
        endif
      endif
c        xi in direction yi x rp
      xi(1)=yi(2)*rp(3)-yi(3)*rp(2)
      xi(2)=yi(3)*rp(1)-yi(1)*rp(3)
      xi(3)=yi(1)*rp(2)-yi(2)*rp(1)
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gphij(rp,cm,np,i,rpi,scmi,cmi,xi,yi,big,ni,phi)
      integer np,i,scmi,ni
      real*8 rp(3,np),cm(np),rpi(3),cmi,xi(3),yi(3),big,phi(2,np)
c
c        local (automatic) variables
      integer j
      real*8 bi,bj,cmij,cmj,d,dc,xj,yj
c *
c * angles phi about z-axis rp(i) of intersection of i & j circles
c * phi = big means no intersection
c *
c  Input: rp, cm
c         np = number of caps
c         i
c         rpi = rp(i)
c         scmi = sign(cm(i))
c         cmi = abs(cm(i))
c         xi, yi form orthonormal axes with rp(i)
c         big
c Output: ni = number of intersections of i circle with j circles
c            = -1 if circle i is outside polygon
c            = -2 if area of polygon is zero
c         phi(1,j), phi(2,j) = azimuthal angle about rpi
c              of intersection of j circle with i circle;
c              zero azimuthal angle is in direction xi;
c
c        initialise phi to big, meaning no intersection
      do j=1,np
        phi(1,j)=big
        phi(2,j)=big
      enddo
c        set count of number of intersections with i circle to zero
      ni=0
c        find intersection of i circle with each j circle in turn
      do 150 j=1,np
c        skip self
        if (j.eq.i) goto 150
c        cm(j).ge.2 means include whole sphere, so no intersection
        if (cm(j).ge.2.d0) goto 150
c        cmij = 2 sin^2[th(ij)/2] = 1-cos th(ij)
        cmij=((rpi(1)-rp(1,j))**2+(rpi(2)-rp(2,j))**2
     *    +(rpi(3)-rp(3,j))**2)/2.d0
c        cmj = 1-cos th(j)
c        bj = cj-ci*cij
c        d = 1-ci^2-cj^2-cij^2+2*ci*cj*cij
c        dph = atan(sqrt(d)/bj) is angle from rp(j) to intersection
        cmj=abs(cm(j))
        bj=(cmi-cmj)+cmij*(1.d0-cmi)
        d=-(cmi-cmj)**2+cmij*(2.d0*((cmi+cmj)-cmi*cmj)-cmij)
        dc=cmi**2+cmj**2
c        positive d means i and j circles intersect
c        if d is zero within numerical round off, treat it as zero
        if (d+dc.gt.dc) then
          d=sqrt(d)
c        ph = atan(yj/xj) is angle from xi to rp(j)
          xj=xi(1)*rp(1,j)+xi(2)*rp(2,j)+xi(3)*rp(3,j)
          yj=yi(1)*rp(1,j)+yi(2)*rp(2,j)+yi(3)*rp(3,j)
c        order intersection angles so segment from first to second angle
c        is inside j circle
          if (cm(j).ge.0.d0) then
c        phi(1,j)=ph-dph , phi(2,j)=ph+dph
            phi(1,j)=atan2(yj*bj-xj*d,xj*bj+yj*d)
            phi(2,j)=atan2(yj*bj+xj*d,xj*bj-yj*d)
          elseif (cm(j).lt.0.d0) then
c        phi(1,j)=ph+dph , phi(2,j)=ph-dph
            phi(1,j)=atan2(yj*bj+xj*d,xj*bj-yj*d)
            phi(2,j)=atan2(yj*bj-xj*d,xj*bj+yj*d)
          endif
c        increment count of number of intersections
          ni=ni+2
c        zero d means i and j circles just touch;
c        negative d means i and j circles don't intersect
        elseif (d+dc.le.dc) then
c        bi = ci-cj*cij
          bi=(cmj-cmi)+cmij*(1.d0-cmj)
c        bi=0 means i and j circles coincide, implying also bj=0 and d=0
c        but test both bi and bj to guard against numerics
          if (bi.eq.0.d0.or.bj.eq.0.d0) then
c        null intersection of areas:
c        rp(i) and rp(j) point in same direction
            if (cmij.lt.1.d0) then
              if ((scmi.ge.0.and.cm(j).lt.0.d0)
     *          .or.(scmi.lt.0.and.cm(j).ge.0.d0)) goto 220
c        rp(i) and rp(j) point in opposite directions
            elseif (cmij.gt.1.d0) then
              if ((scmi.ge.0.and.cm(j).ge.0.d0)
     *          .or.(scmi.lt.0.and.cm(j).lt.0.d0)) goto 220
            endif
c        only do one of the two degenerate circles
            if (i.lt.j) goto 210
c        i circle is outside j circle
          elseif ((cm(j).ge.0.d0.and.bj.gt.0.d0)
     *      .or.(cm(j).lt.0.d0.and.bj.lt.0.d0)) then
c        j circle also outside i circle means null intersection area
            if ((scmi.ge.0.and.bi.gt.0.d0)
     *        .or.(scmi.lt.0.and.bi.lt.0.d0)) goto 220
c        skip i circle since it's entirely outside j circle
            goto 210
          endif
        endif
  150 continue
c        normal return
      return
c
c        circle lies outside area
  210 ni=-1
      return
c
c        area is zero
  220 ni=-2
      return
c
      end
c
c-----------------------------------------------------------------------	
      integer function gsegij(rp,cm,np,rpi,scmi,cmi,ni,tol,phi,iord,
     *  jml,jmu,jpl,jpu,jm,jp,km,kp,phm,php,ph,dph)
      integer np,scmi,ni,iord(2*np),jml,jmu,jpl,jpu,jm,jp,km,kp
      real*8 rp(3,np),cm(np),rpi(3),cmi,tol,phi(2,np),phm,php,ph,dph
c
c        parameters
      include 'pi.par'
      real*8 TWOPI
      parameter (TWOPI=2.d0*PI)
c        local (automatic) variables
      integer j,jj,k
      real*8 cmik,cmk,cpsi,cpsim,dphc,si,sk
c        local variables to be saved
      integer jl,ju
      save jl,ju
c *
c * Determine whether next segment of i circle
c * is an edge of the polygon.
c *
c  Input: rp
c         cm
c         rpi = rp(i)
c         scmi = sign(cm(i))
c         cmi = abs(cm(i))
c         np
c         ni = number of intersections of j circles with i circle
c         tol
c         phi = angles of j circles about i circle
c         iord = order of j circles about i circle
c         jpl = 0 on first call
c             = per output of previous call subsequently
c Output: jml, jmu = points jml to jmu are at lower point
c         jpl, jpu = points jpl to jpu are at upper point
c                    jml <= jmu < jpl <= jpu
c                    1 <= jmu <= ni
c                    jml may be <= 0, and jpl & jpu may be >= ni
c         jm = lower point to use
c         jp = upper point to use
c         km = circle at lower limit of segment
c         km = circle at upper limit of segment
c         phm = azimuthal angle at lower limit of segment
c         php = azimuthal angle at upper limit of segment
c         ph = azimuthal angle at centre of segment
c         dph = azimuthal length of segment between limits
c Return value: -1 = error (tol is too large)
c               0 = segment is not edge of polygon
c               1 = segment is edge of polygon
c               2 = done all segments
c
      if (jpl.ne.0) then
c        done
        if (jpl.eq.jl+ni) goto 210
      endif
c        sin th(i)
      si=sqrt(cmi*(2.d0-cmi))
c        dphc = azimuthal angle corresponding to great circle angle tol
      if (tol.gt.PI) goto 300
      dphc=sin(tol/2.d0)/si
      if (dphc.gt.1.d0) goto 300
      dphc=2.d0*asin(dphc)
c        first segment: jml <= 1 <= jmu < jpl <= jpu < jml+ni
      if (jpl.eq.0) then
c        lower point: jml to jmu are all at the same azimuth phi
        jml=1
        do j=jml,ni
          jmu=j
          jp=mod(j,ni)+1
          km=(iord(j)+1)/2
          kp=(iord(jp)+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km)
          php=phi(1+mod(iord(jp)+1,2),kp)
          dph=php-phm
          if (dph.lt.0.d0) dph=dph+TWOPI
          if (dph.gt.dphc) goto 110
        enddo
  110   continue
c        check if lower point goes lower than 1
        do j=ni,jmu+1,-1
          jp=mod(j,ni)+1
          km=(iord(j)+1)/2
          kp=(iord(jp)+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km)
          php=phi(1+mod(iord(jp)+1,2),kp)
          dph=php-phm
          if (dph.lt.0.d0) dph=dph+TWOPI
          if (dph.gt.dphc) goto 120
          jml=j
        enddo
  120   continue
        if (jml.gt.jmu) jml=jml-ni
c        segment should not include all points
        if (jmu-jml+1.ge.ni) then
          print *,'*** gsegij: segment',jml,' to',jmu,
     *      ' should not include all',ni,' intersections'
          goto 300
        endif
c        record range jl to ju of lower point of first segment
        jl=jml
        ju=jmu
c        upper point: jpl to jpu are all at the same azimuth phi
        jpl=jmu+1
        do j=jpl,jl-1+ni
          jpu=j
          jp=mod(j,ni)+1
          km=(iord(j)+1)/2
          kp=(iord(jp)+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km)
          php=phi(1+mod(iord(jp)+1,2),kp)
          dph=php-phm
          if (dph.lt.0.d0) dph=dph+TWOPI
          if (dph.gt.dphc) goto 130
        enddo
  130   continue
c        first segment: jml <= 1 <= jmu < jpl <= jpu < jl+ni
        if (jml.gt.1.or.jml.gt.jmu.or.jmu.lt.1.or.jmu.gt.ni
     *    .or.jmu.ge.jpl.or.jpl.gt.jpu.or.jpu.ge.jl+ni) then
          print *,'1st segment,',ni,' intersections:',jml,jmu,jpl,jpu
        endif
c        subsequent segments: jml <= jmu < jpl <= jpu < ju+ni
      elseif (jpl.ne.0) then
c        lower point: jml to jmu are all at the same azimuth phi
        jml=jpl
        do j=jml,jl-1+ni
          jmu=j
          jp=mod(j,ni)+1
          km=(iord(j)+1)/2
          kp=(iord(jp)+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(j)+1,2),km)
          php=phi(1+mod(iord(jp)+1,2),kp)
          dph=php-phm
          if (dph.lt.0.d0) dph=dph+TWOPI
          if (dph.gt.dphc) goto 140
        enddo
  140   continue
c        upper point: jpl to jpu are all at the same azimuth phi
        jpl=jmu+1
        do j=jpl,ni+ju
          jpu=j
          jm=mod(j-1+ni,ni)+1
          jp=mod(j,ni)+1
          km=(iord(jm)+1)/2
          kp=(iord(jp)+1)/2
c        lower and upper angles of segment
          phm=phi(1+mod(iord(jm)+1,2),km)
          php=phi(1+mod(iord(jp)+1,2),kp)
          dph=php-phm
          if (dph.lt.0.d0) dph=dph+TWOPI
          if (dph.gt.dphc) goto 150
        enddo
  150   continue
c        subsequent segments: 1 < jml <= jmu < jpl <= jpu <= ju+ni
        if (jml.le.1.or.jml.gt.jmu.or.jmu.lt.1.or.jmu.gt.ni
     *    .or.jmu.ge.jpl.or.jpl.gt.jpu.or.jpu.gt.ju+ni) then
          print *,'> 1st segment,',ni,' intersections:',jml,jmu,jpl,jpu
        endif
      endif
c        lower angle(s) must be lower limit of segments jml to jmu
      do j=jml,jmu
        jm=mod(j-1+ni,ni)+1
        if (mod(iord(jm),2).ne.1) goto 200
      enddo
c        upper angle(s) must be upper limit of segments jpl to jpu
      do j=jpl,jpu
        jp=mod(j-1+ni,ni)+1
        if (mod(iord(jp),2).ne.0) goto 200
      enddo
      km=(iord(jmu)+1)/2
      jp=mod(jpl-1+ni,ni)+1
      kp=(iord(jp)+1)/2
      phm=phi(1+mod(iord(jmu)+1,2),km)
      php=phi(1+mod(iord(jp)+1,2),kp)
c        check segment satisfies all conditions
      if (php.ge.phm) then
        do j=1,ni
          if (mod(iord(j),2).eq.0) then
            k=(iord(j)+1)/2
c        require order   k- ...ph- ph+... k+
            if (phi(2,k).gt.phi(1,k)) then
              if (phm.lt.phi(1,k).or.php.gt.phi(2,k)) goto 200
c        or   k+ k- ...ph- ph+...   or   - ...ph- ph+... k+ k-
            elseif (phi(2,k).lt.phi(1,k)) then
              if (phm.lt.phi(1,k).and.php.gt.phi(2,k)) goto 200
            endif
          endif
        enddo
      elseif (php.lt.phm) then
        do j=1,ni
          k=(iord(j)+1)/2
c        require order   ph+... k+ k- ...ph-
          if (phi(2,k).gt.phi(1,k)) goto 200
        enddo
      endif
c        point at lower limit is one with largest exterior angle psi
      if (jml.eq.jmu) then
        jm=mod(jml-1+ni,ni)+1
      else
        cpsim=1.
        do j=jml,jmu
          jj=mod(j-1+ni,ni)+1
          km=(iord(jj)+1)/2
          cmk=abs(cm(km))
          sk=sqrt(cmk*(2.d0-cmk))
c        cmik = 1-cos th(ik)
          cmik=((rpi(1)-rp(1,km))**2+(rpi(2)-rp(2,km))**2
     *      +(rpi(3)-rp(3,km))**2)/2.d0
c        cpsi = (cik-ci*ck)/(si*sk)
          cpsi=(cmi+cmk-cmi*cmk-cmik)/(si*sk)
          if ((scmi.ge.0.and.cm(km).lt.0.d0)
     *      .or.(scmi.le.0.and.cm(km).ge.0.d0)) cpsi=-cpsi
c        choose largest exterior angle psi (smallest cos psi)
          if (cpsi.lt.cpsim) then
            jm=jj
            cpsim=cpsi
          endif
        enddo
      endif
c        point at upper limit is one with largest exterior angle psi
      if (jpl.eq.jpu) then
        jp=mod(jpl-1+ni,ni)+1
      else
        cpsim=1.d0
        do j=jpl,jpu
          jj=mod(j-1+ni,ni)+1
          kp=(iord(jj)+1)/2
          cmk=abs(cm(kp))
          sk=sqrt(cmk*(2.d0-cmk))
c        cmik = 1-cos th(ik)
          cmik=((rpi(1)-rp(1,kp))**2+(rpi(2)-rp(2,kp))**2
     *      +(rpi(3)-rp(3,kp))**2)/2.d0
c        cpsi = (cik-ci*ck)/(si*sk)
          cpsi=(cmi+cmk-cmi*cmk-cmik)/(si*sk)
          if ((scmi.ge.0.and.cm(kp).lt.0.d0)
     *      .or.(scmi.lt.0.and.cm(kp).ge.0.d0)) cpsi=-cpsi
c        choose largest exterior angle psi (smallest cos psi)
          if (cpsi.lt.cpsim) then
            jp=jj
            cpsim=cpsi
          endif
        enddo
      endif
c        circles at upper and lower limits
      km=(iord(jm)+1)/2
      kp=(iord(jp)+1)/2
      phm=phi(1+mod(iord(jm)+1,2),km)
      php=phi(1+mod(iord(jp)+1,2),kp)
c        angular length, centre point of segment
      if (php.ge.phm) then
        dph=php-phm
        ph=(php+phm)/2.d0
      elseif (php.lt.phm) then
        ph=(php+phm)/2.d0
        if (ph.le.0.d0) then
          ph=ph+PI
          php=php+TWOPI
        elseif (ph.gt.0.d0) then
          ph=ph-PI
          phm=phm-TWOPI
        endif
        dph=php-phm
      endif
c      segment is edge of polygon
      gsegij=1
      return
c
c      segment is not edge of polygon
  200 gsegij=0
      return
c
c      done full circle
  210 gsegij=2
      return
c
c      tol is too large
  300 gsegij=-1
      return
c
      end
c
c-----------------------------------------------------------------------
      subroutine gvtrail(scmi,np,i,km,kp,vtrail,vik,nvmax,nv,ik,ikchk)
      integer scmi,np,i,km,kp,nvmax,vtrail(nvmax),vik(nvmax,2),nv,ik
      real*8 ikchk
c
c        externals
      integer ikrand,ik2i,ik2k,ik2ik
c        local (automatic) variables
      integer ikran,iv,k,l
c *
c * Record endpoints of edge.
c *
c        two end points of the edge, going right-handedly around edge
      do 250 l=1,2
c        end point is intersection of i circle with k circle
        if (l.eq.1.and.scmi.gt.0.or.l.eq.2.and.scmi.lt.0) then
          k=km
        else
          k=kp
        endif
c        from edge k to edge i right-handedly through vertex
        if (l.eq.1) then
          ik=ik2ik(np,k,i)
c        from edge i to edge k right-handedly through vertex
        else
          ik=ik2ik(np,i,k)
        endif
        if (nv.le.nvmax) then
          vik(nv,l)=ik
        endif
c        pseudo-random number from ik
        ikran=ikrand(ik)
c        new intersection
        if (i.lt.k) then
          ikchk=ikchk+ikran
c        intersection already met
        else
          ikchk=ikchk-ikran
          if (nv.le.nvmax) then
c        index of vertex already met
            if (l.eq.1) then
              do iv=1,nv
                if (vik(iv,2).eq.ik) then
c        vertex after iv is nv
                  vtrail(iv)=nv
                  goto 240
                endif
              enddo
            elseif (l.eq.2) then
              do iv=1,nv
                if (vik(iv,1).eq.ik) then
c        vertex after nv is iv
                  vtrail(nv)=iv
                  goto 240
                endif
              enddo
            endif
c        should not happen
            print *,'*** from gvtrail: no vertex',
     *        ' (',ik2i(np,ik),ik2k(np,ik),')',
     *        ' found at',i,', edge',km,kp,'???'
            print *,'vertices so far are:',
     *        (' (',ik2k(np,vik(iv,1)),ik2i(np,vik(iv,1)),')',
     *         ' (',ik2i(np,vik(iv,2)),ik2k(np,vik(iv,2)),')',
     *        iv=1,nv)
  240       continue
          endif
        endif
  250 continue
      return
      end
c
c-----------------------------------------------------------------------
      integer function ik2ik(np,i,k)
      integer np,i,k
c
      integer ik
c
c        gives ik = 1,2,...,np(np-1) for i,k = 0,1,...,np with i != k
      ik=i+np*k
      if (i.lt.k) ik=ik+1
      ik2ik=ik
      return
      end
c
c-----------------------------------------------------------------------
      integer function ik2i(np,ik)
      integer np,ik
c
      integer i,k
c
      i=mod(ik-1,np)
      k=(ik-1)/np
      if (i.ge.k) i=i+1
      ik2i=i
      return
      end
c
c-----------------------------------------------------------------------
      integer function ik2k(np,ik)
      integer np,ik
c
      integer k
c
      k=(ik-1)/np
      ik2k=k
      return
      end
c
c-----------------------------------------------------------------------
      subroutine gvord(np,nv,vtrail,vord,vik,ev,nev)
      integer np,nv,vtrail(nv),vord(nv),vik(nv),ev(nv),nev
c
c        externals
      integer ik2i,ik2k
c        local (automatic) variables
      integer i,iev,inter,j,k,nvm
c *
c * Order of vertices around polygon.
c *
c  Input: nv = number of vertices
c         vtrail = vertex trail
c Output: vord = vertex order
c         ev = end indices of each connected sequence of vertices
c         nev = number of connected sequences of vertices
c
      nev=0
      if (nv.gt.0) then
c       print *,'gvord: vertex trail:',(vtrail(i),i=1,nv)
c       print *,'gvord: vik_1 on input:',(vord(i),i=1,nv)
c       print *,'gvord: vik_2 on input:',(vik(i),i=1,nv)
        do i=1,nv
          vord(i)=0
        enddo
        j=0
        iev=0
c        do intersecting, then non-intersecting circles
        do inter=1,2
          do 320 nvm=1,nv
c        find first vertex not already traversed
            do i=1,nv
              if (vord(i).eq.nvm) goto 320
            enddo
            if (inter.eq.1) then
              if (ik2i(np,vik(nvm)).eq.0.or.ik2k(np,vik(nvm)).eq.0)
     *          goto 320
            else 
              if (ik2i(np,vik(nvm)).ne.0.and.ik2k(np,vik(nvm)).ne.0)
     *          goto 320
            endif
c        index of first vertex of this boundary
            j=j+1
            vord(j)=nvm
c        circulate
            k=j+1
            do j=k,nv
c        next vertex of this boundary
              vord(j)=vtrail(vord(j-1))
c        gone full circle
              if (vtrail(vord(j)).eq.nvm) then
c        record cumulative length
                iev=iev+1
                ev(iev)=j
                nev=nev+1
                if (j.eq.nv) goto 330
                goto 320
              endif
            enddo
  320     continue
        enddo
  330   continue
c       print *,'vertex order:',(vord(i),i=1,nv)
c        should not happen
        if (vtrail(vord(nv)).ne.nvm) then
          print *,'*** from gvord: last vertex',vord(nv),
     *      ' should connect to',nvm,'???'
        endif
      endif
      return
c
      end
c
c-----------------------------------------------------------------------
      real*8 function cmijf(rpi,rpj)
      real*8 rpi(3),rpj(3)
c
      cmijf=((rpi(1)-rpj(1))**2+(rpi(2)-rpj(2))**2
     *  +(rpi(3)-rpj(3))**2)/2.d0
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine vpermi(ix,n,iperm,iwk)
      integer n,iperm(n)
      integer ix(n),iwk(n)
c
c        local (automatic) variables
      integer i,j
c *
c * Permutes elements of integer array ix(n) by permutation iperm(n).
c *
      do i=1,n
        iwk(i)=ix(i)
      enddo
      do i=1,n
        j=iperm(i)
        ix(i)=iwk(j)
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      subroutine vpermd(x,n,iperm,wk)
      integer n,iperm(n)
      real*8 x(n),wk(n)
c
c        local (automatic) variables
      integer i,j
c *
c * Permutes elements of real*8 array x(n) by permutation iperm(n).
c *
      do i=1,n
        wk(i)=x(i)
      enddo
      do i=1,n
        j=iperm(i)
        x(i)=wk(j)
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      subroutine vpermdd(x,m,n,iperm,wk)
      integer m,n,iperm(n)
      real*8 x(m,n),wk(n)
c
c        local (automatic) variables
      integer i,j,k
c *
c * Permutes each column of real*8 array x(m,n) by permutation iperm(n).
c *
      do k=1,m
        do i=1,n
          wk(i)=x(k,i)
        enddo
        do i=1,n
          j=iperm(i)
          x(k,i)=wk(j)
        enddo
      enddo
      return
      end
c
