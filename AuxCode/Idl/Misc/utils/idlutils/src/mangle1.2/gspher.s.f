c-----------------------------------------------------------------------
      subroutine gspher(area,bound,vert,w,lmax1,im,nw,rp,cm,np,npc,ibv,
     *  iphi,tol,phw,iord,v,ldegen)
      integer lmax1,im,nw,np,npc,ibv,iphi
      logical ldegen
      real*8 area,bound(2),vert(2),w(im,nw),rp(3,np),cm(np),tol
c        work arrays (could be automatic if compiler supports it)
      integer iord(2*np)
      real*8 phw(2,np),v(lmax1)
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
c        local (automatic) variables
      integer i,icmmin,ik,iseg,j,jm,jml,jmu,jp,jpl,jpu,k,km,kp,l,
     *  nbd,ni,scmi
      logical warn,whole
      real*8 ci,cmi,cmik,cmmin,cmk,cpsi,cti,ctk,
     *  darea,dbound(2),dph,dvert(2),ikchk,
     *  p,ph,phi,phii,phm,php,psi,ri,rii,si,sk,spsi,sqrt4pi,xi(3),yi(3)
c *
c * Spherical transform of region W of sphere of unit radius bounded by
c *    1 - r.rp(i) < cm(i)  (if cm(i).ge.0)
c *    1 - r.rp(i) > -cm(i)  (if cm(i).lt.0)
c * for i=1,np where rp(i) are unit directions.
c * See AJSH notes Multfn C115.
c *
c * bound (boundary) and vert (vertex) are related to the coefficients
c * of the 1st & 2nd order terms of the power series expansion
c * in sin(th/2) of the correlation <WW> at angular separation th:
c *    <WW> = sum(lm) |W_lm|^2 P_l(cos th)
c *         = 2*pi*area - 4*bound*sin(th/2) + 2*vert*sin^2(th/2) + ...
c * In fact bound is just the total length of the boundary in radians,
c * and vert is the sum over vertex terms of 1-psi/tan(psi), where
c * psi is the exterior angle (=pi-interior angle) at the vertex.
c * The range of validity of the power series expansion in sin(th/2)
c * depends on the tortuousness of the boundary: the more convoluted the
c * boundary, and the sharper the vertex angles, the shorter the range.
c *
c * Oct 1992: next order term incorporated
c *    <WW> = 2*pi*area - 4*bound(1)*sin(th/2) + 2*vert(1)*sin^2(th/2)
c *           + [2/3*bound(2) + 8/9*vert(2)]*sin^3(th/2) + ...
c *
c * The option ibv modifies the evaluation of bound and vert, and
c * determines whether harmonics are added to or subtracted from input:
c * ibv = 0: standard option, for <WW>.
c *          Harmonics of W are ADDED to input harmonics.
c * ibv = 1: to take the cross-correlation between a region W1 and
c *          its intersection W12 with region W2.
c *          Here <W12 W1> = <W12^2> + <W12(W1-W12)>
c *          and this is what the code evaluates bound and vert for.
c *          Harmonics of W12 are ADDED to input harmonics.
c * ibv = 2: to take region W1 less its intersection W12 with region W2.
c *          Here <WW> = <(W1-W12)(W1-W12)>
c *                    = <W1^2> - <W12^2> - 2<W12(W1-W12)>
c *          The code evaluates bound and vert for
c *                               <W12^2> + 2<W12(W1-W12)>.
c *          Harmonics of W12 are SUBTRACTED from input harmonics.
c * ibv = 3: to take the union of regions W1 and W2.
c *          Here <WW> = <(W1+W2-W12)(W1+W2-W12)>
c *                    = <W1^2> + <W2^2> - <W12^2> + 2<(W1-W12)(W2-W12)>
c *          The code evaluates bound and vert for
c *                               <W12^2> - 2<(W1-W12)(W2-W12)>.
c *          Harmonics of W12 are SUBTRACTED from input harmonics.
c * The correct calling procedure in each case is exampled below.
c *
c * Complement of a region:
c *    area -> 4*pi-area,  w(1,1) -> sqrt(4*pi)-w(1,1)
c *    bound -> bound
c *    vert -> vert
c *    w -> -w   except monopole term w(1,1) as above
c *
c * To get intersection W12 of regions W1 & W2,
c * but evaluate bound & vert for cross-correlation between W12 and W1:
c *    put the np2 constraints of region W2 in 1 to npc
c *      & the np1 constraints of region W1 in npc+1 to np
c *    call gspher(area,bound,vert,w,parameters of W1,npc=np2,1,...)
c *
c * To get region W1 less its intersection with region W2:
c *    call gspher(area,bound,vert,w,parameters of W1,npc=0,0,...)
c *    put the np2 constraints of region W2 in 1 to npc
c *      & the np1 constraints of region W1 in npc+1 to np
c *    call gspher(darea,dbound,dvert,w,params of W1 & W2,npc=np2,2,...)
c *    area=area-darea
c *    bound=bound-dbound
c *    vert=vert-dvert
c *
c * To get the union of two regions W1 & W2:
c *    call gspher(area,bound,vert,w,parameters of W1,npc=0,0,...)
c *    call gspher(darea,dbound,dvert,w,parameters of W2,npc=0,0,...)
c *    area=area+darea
c *    bound=bound+dbound
c *    vert=dvert+dvert
c *    put the np2 constraints of region W2 in 1 to npc
c *      & the np1 constraints of region W1 in npc+1 to np
c *    call gspher(darea,dbound,dvert,w,params of W1 & W2,npc=np2,3,...)
c *    area=area-darea
c *    bound=bound-dbound
c *    vert=dvert-dvert
c *
c * Cautions:
c * (1) This subroutine underestimates the area (monopole harmonic)
c *     by 2*pi (sqrt(pi)) if 2*pi <= area < 4*pi
c *     and the area is bounded by more than one arc.
c * (2) This subroutine will usually work correctly when there are near
c *     multiple (.ge. 3) intersections of boundaries, but in rare
c *     instances it may fail.  If so, it should flag the failure with
c *     ldegen=.true.  This error condition should NOT be ignored.
c * (3) There is a mathematical ambiguity in the boundary term bound
c *     whenever two (or more) boundaries coincide.  To resolve the
c *     ambiguity correctly, perturb the boundaries slightly.
c *     Likewise there is a mathematical ambiguity in the vertex term
c *     vert whenever there are multiple (.ge. 3) intersections.
c *     Again, to resolve the ambiguity correctly, perturb the
c *     boundaries slightly.
c *     It should be noted that if these problems exist, then the series
c *     expansion of <WW>, whose coefficients involve bound and vert,
c *     probably breaks down already at tiny values of the separation
c *     angle th.
c *     
c  Input: lmax1 = lmax+1 where lmax is maximum desired l of transform.
c         im = 1 means compute only real part of harmonics;
c              2 means compute both real and imaginary parts.
c              Note harmonics are real if region possesses reflection
c              symmetry through plane defined by z-axis and direction
c              rp(iphi).
c         nw = [(lmax+1)*(lmax+2)]/2
c         rp(3,i),i=1,np = x, y, z coordinates of a set of unit
c              directions defining the region.
c              It is assumed without checking that
c              rp(1)^2 + rp(2)^2 + rp(3)^2 = 1 .
c         cm(i),i=1,np = set of 1-cos's defining the region.
c         np = number of directions.
c         npc = part number of directions, used in conjunction with ibv;
c               see above for more details.
c         ibv = 0 to 3 controls evaluation of bound & vert;
c               see above for more details.
c         iphi > 0 means compute harmonics in frame of reference where
c                  y-axis is along z x rp(iphi);
c              = 0 means use the input frame of reference.
c Output: area = area of region in steradians.
c         bound = length of boundary of region in radians if ibv=0,
c                 or as explained above if ibv>0.
c         vert = sum over vertices of 1-psi/tan(psi) if ibv=0,
c                where psi is exterior angle (=pi-interior angle)
c                at vertex, or as explained above if ibv>0.
c         ldegen = .true. signals an error: the code dealt incorrectly
c                  with a multiply intersecting boundary.
c Input/Output: w(i,lm) = spherical transform, dimensioned w(im,nw)
c            w(i,lm), i=1,im, lm=l*(l+1)/2+m+1, l=0,lmax, m=0,l;
c            w(1,lm) is real part, w(2,lm) is imaginary part (if im=2).
c            w(i,lm) is ADDED to input w(i,lm) if ibv=0 or 1,
c                    or SUBTRACTED from input w(i,lm) if ibv=2 or 3.
c            Note w(l,-m)=(-)**m*[Complex conjugate of w(l,m)], just as
c                 Y(l,-m)=(-)**m*[Complex conjugate of Y(l,m)].
c            tol
c Work arrays: phw and iord should be dimensioned at least 2*np.
c              v should be dimensioned at least lmax1.
c
c        set azimuthal angle of non-intersection to big
      data big /1.d6/
c        ok if area tests not too far outside [0,max]
      data arearnd /1.d-10/
c        warn about multiple intersection when dph < dphmin
c     data dphmin /1.d-8/
c
c        come here with enlarged tolerance
  100 continue
c        initialise error flag to no error
      ldegen=.false.
      warn=.false.
c        zero stuff
      area=0.d0
      bound(1)=0.d0
      bound(2)=0.d0
      vert(1)=0.d0
      vert(2)=0.d0
c     do j=1,nw
c       do i=1,im
c         w(i,j)=0.d0
c       enddo
c     enddo
c        check for zero area because one circle is null
      if (gzeroar(cm,np).eq.0) goto 410
c        no constraints at all will mean area is whole sphere
      whole=.true.
c        number of arc segments bounding area
      nbd=0
c        error check on evaluation of vertex terms
      ikchk=0.d0
c        area=sqrt(4pi)*monopole
      sqrt4pi=sqrt(4.d0*PI)
c        harmonics defined so point iphi is at zero azimuthal angle
      rii=0.d0
      if (iphi.ge.1) rii=sqrt(rp(1,iphi)**2+rp(2,iphi)**2)
      phii=0.d0
      if (rii.gt.0.d0) phii=atan2(rp(2,iphi),rp(1,iphi))
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
c        ci = cos th(i)
        ci=1.d0-cmi
c        si = sin th(i)
        si=sqrt(cmi*(2.d0-cmi))
c........ri, phi, rp(3,i) are cylindrical coordinates of rp(i)
        ri=sqrt(rp(1,i)**2+rp(2,i)**2)
        if (ri.eq.0.d0.or.i.eq.iphi) then
          phi=0.d0
        else
          phi=atan2(rp(2,i),rp(1,i))-phii
        endif
c........construct cartesian axes with z-axis along rp(i)
c        set yi in direction z x rp(i)
        if (ri.gt.0.d0) then
          yi(1)=-rp(2,i)/ri
          yi(2)=rp(1,i)/ri
          yi(3)=0.d0
c        if rp(i) is along z-axis, set yi in direction z x rp(iphi)
        elseif (rii.gt.0.d0) then
          yi(1)=-rp(2,iphi)/rii
          yi(2)=rp(1,iphi)/rii
          yi(3)=0.d0
c        if rp(iphi) is also along z-axis, set yi along y-axis
        elseif (ri.eq.0.d0.and.rii.eq.0.d0) then
          yi(1)=0.d0
          yi(2)=1.d0
          yi(3)=0.d0
        endif
c        xi in direction yi x rp(i)
        xi(1)=yi(2)*rp(3,i)-yi(3)*rp(2,i)
        xi(2)=yi(3)*rp(1,i)-yi(1)*rp(3,i)
        xi(3)=yi(1)*rp(2,i)-yi(2)*rp(1,i)
c........angles phi about z-axis rp(i) of intersection of i & j circles
        call gphij(rp,cm,np,i,rp(1,i),scmi,cmi,xi,yi,big,ni,phw)
c        i circle lies outside polygon
        if (ni.eq.-1) goto 280
c        area of polygon is zero
        if (ni.eq.-2) goto 410 
c........i circle has no intersections
        if (ni.eq.0) then
          nbd=nbd+1
          dph=TWOPI
          if (cm(i).lt.0.d0) dph=-dph
          ph=0.d0
c        increment area
          darea=cmi*dph
          area=area+darea
c        bound(1) term is length of boundary
          dbound(1)=si*abs(dph)
          dbound(2)=(1.d0/si-2.d0*si)*abs(dph)
c        standard
          if (ibv.eq.0
c        cross
     *      .or.(ibv.eq.1.and.i.gt.npc)
c        intersection
     *      .or.(ibv.eq.2.and.i.gt.npc)
c        union
     *      .or.(ibv.eq.3)) then
            bound(1)=bound(1)+dbound(1)
            bound(2)=bound(2)+dbound(2)
c        intersection
          elseif (ibv.eq.2.and.i.le.npc) then
            bound(1)=bound(1)-dbound(1)
            bound(2)=bound(2)-dbound(2)
          endif
c         print *,'at',i,': full circle area +=',darea,' =',area
c        increment spherical transform
          if (ibv.ge.2) dph=-dph
          call wlm(w,lmax1,im,nw,ri,phi,0,rp(3,i),ci,si,ph,dph,v)
c........i circle has intersections
        elseif (ni.gt.0) then
c        find ordering of intersection angles around i circle
          call findbot(phw,2*np,iord,ni)
c........contribution from each segment of i circle
          jpl=0
c        come here to do another segment
  200     continue
c........is segment edge of polygon?
            iseg=gsegij(rp,cm,np,rp(1,i),scmi,cmi,ni,tol,phw,iord,
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
c    *          '*** warning from gspher: near multiple intersection at'
c    *          ,i,': segment',km,kp,' dph=',dph
c             warn=.true.
c           endif
c........segment satisfies conditions
            nbd=nbd+1
c        contribution to area from the boundary segment
            if (cm(i).lt.0.d0) dph=-dph
            darea=cmi*dph-dph
            area=area+darea
c        bound(1) term is length of boundary
            dbound(1)=si*abs(dph)
            dbound(2)=(1.d0/si-2.d0*si)*abs(dph)
c        standard
            if (ibv.eq.0
c        cross
     *        .or.(ibv.eq.1.and.i.gt.npc)
c        intersection
     *        .or.(ibv.eq.2.and.i.gt.npc)
c        union
     *        .or.(ibv.eq.3)) then
              bound(1)=bound(1)+dbound(1)
              bound(2)=bound(2)+dbound(2)
c        intersection
            elseif (ibv.eq.2.and.i.le.npc) then
              bound(1)=bound(1)-dbound(1)
              bound(2)=bound(2)-dbound(2)
            endif
c           print *,'at',i,': edge',km,kp,
c    *       ' (',jm,' in',jml,jmu,',',jp,' in',jpl,jpu,' of',ni,')'
c           print *,'dph =',dph,' area +=',darea,' =',area
c        contribution to area and vert from end points of the segment
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
              if (psi.ne.0.d0) then
c        increment area
                area=area-psi
                spsi=sqrt(1.d0-cpsi**2)
                cti=(1.d0-cmi)/si
                if (cm(i).lt.0.d0) cti=-cti
                ctk=(1.d0-cmk)/sk
                if (cm(k).lt.0.d0) ctk=-ctk
c        standard
                if (ibv.eq.0
     *            .or.(ibv.eq.1.and.i.gt.npc.and.k.gt.npc)
     *            .or.(ibv.eq.2.and.i.gt.npc.and.k.gt.npc)
     *            .or.(ibv.eq.3.and.((i.gt.npc.and.k.gt.npc)
     *                           .or.(i.le.npc.and.k.le.npc)))) then
                  dvert(1)=1.d0-psi*cpsi/spsi
                  dvert(2)=spsi/(1.d0+cpsi)**2*(2.d0+cpsi)*(cti+ctk)
                  vert(1)=vert(1)+dvert(1)
                  vert(2)=vert(2)+dvert(2)
c        cross
                elseif (ibv.eq.1) then
                  if (i.le.npc.and.k.le.npc) then
                    continue
                  else
                    dvert(1)=PI/2.d0*cpsi/spsi
                    dvert(2)=spsi/(1.d0+cpsi)**2*(2.d0+cpsi)*(cti+ctk)
                    vert(1)=vert(1)-dvert(1)
                    vert(2)=vert(2)+dvert(2)/2.d0
                    dvert(2)=spsi/(1.d0-cpsi)**2*(2.d0-cpsi)*(cti-ctk)
                    if (i.gt.npc) then
                      vert(2)=vert(2)-dvert(2)/2.d0
                    elseif (k.gt.npc) then
                      vert(2)=vert(2)+dvert(2)/2.d0
                    endif
                  endif
c        intersection
                elseif (ibv.eq.2) then
                  if (i.le.npc.and.k.le.npc) then
                    dvert(1)=1.d0-psi*cpsi/spsi
                    dvert(2)=spsi/(1.d0+cpsi)**2*(2.d0+cpsi)*(cti+ctk)
                    vert(1)=vert(1)-dvert(1)
                    vert(2)=vert(2)-dvert(2)
                  else
                    dvert(1)=1.d0+(PI-psi)*cpsi/spsi
                    dvert(2)=spsi/(1.d0-cpsi)**2*(2.d0-cpsi)*(cti-ctk)
                    vert(1)=vert(1)-dvert(1)
                    if (i.gt.npc) then
                      vert(2)=vert(2)-dvert(2)
                    elseif (k.gt.npc) then
                      vert(2)=vert(2)+dvert(2)
                    endif
                  endif
c        union
                elseif (ibv.eq.3) then
                  dvert(1)=1.d0-psi*cpsi/spsi
                  dvert(2)=spsi/(1.d0+cpsi)**2*(2.d0+cpsi)*(cti+ctk)
                  vert(1)=vert(1)-dvert(1)
                  vert(2)=vert(2)+dvert(2)
                endif
c               print *,'     intersect',i,k,' area +=',-psi,' =',area
c        peculiar monopole term
                if (ibv.ge.2) psi=-psi
                w(1,1)=w(1,1)-psi/sqrt4pi
              endif
  240       continue
c        increment spherical transform
            if (ibv.ge.2) dph=-dph
            call wlm(w,lmax1,im,nw,ri,phi,0,rp(3,i),ci,si,ph,dph,v)
c        do another segment
          goto 200
        endif
  280 continue
c--------finish off
c        whole=.true. means area is whole sphere
      if (whole) then
        area=2.d0*TWOPI
        if (ibv.eq.0.or.ibv.eq.1) then
          w(1,1)=w(1,1)+2.d0*(TWOPI/sqrt4pi)
        elseif (ibv.eq.2.or.ibv.eq.3) then
          w(1,1)=w(1,1)-2.d0*(TWOPI/sqrt4pi)
        endif
        goto 410
c        just one boundary segment
      elseif (nbd.eq.1) then
        if (area.lt.0.d0) area=area+2.d0*TWOPI
        if (ibv.eq.0.or.ibv.eq.1) then
          w(1,1)=w(1,1)+TWOPI/sqrt4pi
        elseif (ibv.eq.2.or.ibv.eq.3) then
          w(1,1)=w(1,1)-TWOPI/sqrt4pi
        endif
        goto 410
      endif
c        add/subtract 2*pi's to area to ensure 0.le.area.lt.2*pi
      i=area/TWOPI
      area=area-i*TWOPI
      if (area.lt.0.d0) area=area+TWOPI
      j=w(1,1)/(TWOPI/sqrt4pi)
      w(1,1)=w(1,1)-j*(TWOPI/sqrt4pi)
      if (w(1,1).lt.0.d0) w(1,1)=w(1,1)+(TWOPI/sqrt4pi)
c     print *,'area adjusted by',-i,'  monopole by',-j
c        check on whether ik endpoints matched ki endpoints
      if (ikchk.ne.0.d0) then
c       print *,'*** from gspher: at tol =',tol,
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
c        check area is positive
      if (area.lt.0.d0) then
c        negative area should `never' happen
        print *,'*** from gspher: area=',area,' should be .ge. 0'
        warn=.true.
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
      if (p.gt.1.d0) then
c        check if discrepancy is from numerical roundoff
        if (abs(area-TWOPI).le.arearnd) then
          area=0.d0
          bound(1)=0.d0
          bound(2)=0.d0
          vert(1)=0.d0
          vert(2)=0.d0
          goto 410
        endif
        if (area.le.darea+arearnd) then
          area=darea
          goto 410
        endif
        print *,'*** from garea: area/area(',icmmin,')=',
     *    area,' /',darea,' =',p,'  should be .le. 1'
        ldegen=.true.
        warn=.true.
        goto 420
      endif
      if (warn) then
        write (*,'(a3,a20,4a24)')
     *    ' ','x','y','z','r','1-c'
        write (*,'(i3,5g24.16)')
     *    (j,(rp(i,j),i=1,3),sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *    cm(j),j=1,np)
      endif
  410 continue
      return
c
  420 print *,'*** from gspher: total failure'
      ldegen=.true.
      write (*,'(a3,a20,4a24)')
     *  ' ','x','y','z','r','1-c'
      write (*,'(i3,5g24.16)')
     *  (j,(rp(i,j),i=1,3),sqrt(rp(1,j)**2+rp(2,j)**2+rp(3,j)**2),
     *  cm(j),j=1,np)
      return
c
      end
c
