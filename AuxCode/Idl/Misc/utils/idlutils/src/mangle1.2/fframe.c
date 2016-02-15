/*------------------------------------------------------------------------------
*/
frame(fmt, azi, eli, azf, elf)
format *fmt;
double *azn, *eln, *azp;
{
    char *system;
    char fk4[] = "FK4", fk5[] = "FK5";
    int epoch1950 = 1950., epoch2000 = 2000.;

    /* initial and final frames identical */
    if (fmt->inframe == fmt->outframe && fmt->inepoch == fmt->outepoch) {
	azf = elf;
	azi = eli;
    } else {

c--------initialize
      if (init) then
c        ecliptic latitude of North Celestial Pole (1950 FK4)
        elp=felp(BEPOCH)
c        angles to transform between ecliptic and galactic coordinates
        call azell(RAG,DECG,L2P,RAEZ,elp,EAZP,azg,elg,l2z)
        init=.false.
      endif

/* 1950 frames -> */
	if (strcmp(fmt->inframe, 'eqB') == 0 || strcmp(fmt->inframe, eclipticB) == 0) {
	    system = fk4;
	    azf = azi;
	    elf = eli;
	    /* precess */
	    SLA_preces(system, &fmt->inepoch, &epoch1950, &azf, &elf)

	    /* RA & Dec B1950.0 FK4 -> */
	    if (strcmp(fmt->inframe, 'eqB') == 0) {
		if (strcmp(fmt->inframe, 'eclipticB') == 0) {
		    azn = RAEZ;
		    eln = felp_(&fmt->inepoch);
		    azp = EAZP;
		} else if (strcmp(fmt->inframe, 'galactic') == 0) {
		    azn = RAG;
		    eln = DECG;
		    azp = L2P;
		} else {
		    azi = 0.;	eli = 90. / RADIAN;
		    /* B1950.0 FK4 -> J2000.0 FK5 assuming zero proper motion in the FK5 frame */
		    SLA_fk45z_(&azi, &eli, &fmt->inepoch, &ra2k, &dec2k);
		    ra2k *= RADIAN;
		    dec2k *= RADIAN;
		}
	    /* Ecliptic 1950 -> */
	    } else if (strcmp(fmt->inframe, 'eclipticB') == 0) {
		if (strcmp(fmt->inframe, 'eqB') == 0) {
		    azn = EAZP;
		    eln = felp_(&fmt->inepoch);
		    azp = RAEZ;
		} else if (strcmp(fmt->inframe, 'galactic') == 0) {
		    elp = felp_(&fmt->inepoch);
		    azell_(&RAG, &DECG, &L2P, &RAEZ, &elp, &EAZP, &azn, &eln, &azp);
		} else {
		    azi = 0.;	eli = 90. / RADIAN;
		    date = SLA_epj2d(fmt->inepoch);
		    SLA_ecleq_(&azi, &eli, &date, &ra2k, &dec2k);
		    ra2k *= RADIAN;
		    dec2k *= RADIAN;
	    }
	    /* 2000 frames -> */
	    if (strcmp(fmt->outframe, 'eqJ') == 0) {
	        azf=ra2k
	        elf=dec2k
	    } else if (strcmp(fmt->outframe, 'eclipticJ') == 0) {
	        date = SLA_epj2d(JEPOCH)
	        call SLA_eqecl(ra2k/RADIAN,dec2k/RADIAN,date,azf,elf)
	        azf *= RADIAN
	        elf *= RADIAN
	    } else if (strcmp(fmt->outframe, 'SDSS') == 0) {
	        call azel(ra2k,dec2k,
		azn = RASDNP;
		eln = DECSDNP;
		azo = -ETANCP;
	    }

c--------2000 frames ->
      elseif (framei.eq.EQUATORIAL2K.or.framei.eq.ECLIPTIC2K
     *  .or.framei.eq.SDSS) then
c........initial frame -> RA & Dec 2000 FK5
        if (framei.eq.EQUATORIAL2K) then
          ra2k=azi
          dec2k=eli
        elseif (framei.eq.ECLIPTIC2K) then
          date=SLA_epj2d(JEPOCH)
          call SLA_ecleq(azi/RADIAN,eli/RADIAN,date,ra2k,dec2k)
          ra2k=ra2k*RADIAN
          dec2k=dec2k*RADIAN
        elseif (framei.eq.SDSS) then
          call azel(azi,eli,-ETANCP,DECSDNP,RASDNP,ra2k,dec2k)
        endif
c........RA & Dec 2000 FK5 -> final frame
        if (framef.eq.EQUATORIAL) then
          call SLA_fk54z(ra2k/RADIAN,dec2k/RADIAN,BEPOCH,azf,elf,dr,dd)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.ECLIPTIC) then
          date=SLA_epj2d(BEPOCH)
          call SLA_eqecl(ra2k/RADIAN,dec2k/RADIAN,date,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.GALACTIC) then
          call SLA_eqgal(ra2k/RADIAN,dec2k/RADIAN,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.EQUATORIAL2K) then
          azf=ra2k
          elf=dec2k
        elseif (framef.eq.ECLIPTIC2K) then
          date=SLA_epj2d(JEPOCH)
          call SLA_eqecl(ra2k/RADIAN,dec2k/RADIAN,date,azf,elf)
          azf=azf*RADIAN
          elf=elf*RADIAN
        elseif (framef.eq.SDSS) then
          call azel(ra2k,dec2k,RASDNP,DECSDNP,-ETANCP,azf,elf)
        endif

c--------Galactic frame ->
      elseif (framei.eq.GALACTIC) then
        if (framef.eq.EQUATORIAL) then
c         call SLA_ge50(azi/RADIAN,eli/RADIAN,azf,elf)
c         azf=azf*RADIAN
c         elf=elf*RADIAN
          call azel(azi,eli,L2P,DECG,RAG,azf,elf)
        elseif (framef.eq.ECLIPTIC) then
          call azel(azi,eli,l2z,elg,azg,azf,elf)
        else
          call SLA_galeq(azi/RADIAN,eli/RADIAN,ra2k,dec2k)
          if (framef.eq.EQUATORIAL2K) then
            azf=ra2k*RADIAN
            elf=dec2k*RADIAN
          elseif (framef.eq.ECLIPTIC2K) then
            date=SLA_epj2d(JEPOCH)
            call SLA_eqecl(ra2k,dec2k,date,azf,elf)
            azf=azf*RADIAN
            elf=elf*RADIAN
          elseif (framef.eq.SDSS) then
            ra2k=ra2k*RADIAN
            dec2k=dec2k*RADIAN
            call azel(ra2k,dec2k,RASDNP,DECSDNP,-ETANCP,azf,elf)
          endif
        endif

      endif
c--------put azf in interval [0,360)

      iaz=azf/360.d0
      if (azf.lt.0.d0) iaz=iaz-1
      azf=azf-iaz*360.d0

      return
      end
c
