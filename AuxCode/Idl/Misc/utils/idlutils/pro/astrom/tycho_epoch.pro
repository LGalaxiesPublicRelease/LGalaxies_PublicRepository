;------------------------------------------------------------------------------
;+
; NAME:
;   tycho_epoch
;
; PURPOSE:
;   Apply proper motion corrections to an epoch other than 2000
;   for the Tycho-2 catalog.
;
; CALLING SEQUENCE:
;   tycho_epoch, epoch, tycdat
;
; INPUTS:
;   epoch:       New epoch
;   tycdat:      Tycho-2 data structure
;
; OUTPUTS:
;   tycdat:      (Modified)
;
; COMMENTS:
;   The fields RAMDEG,DEMDEG are assumed to be the epoch 2000 positions
;   in degrees.  These positions are moved according to proper motions
;   given by PMRA,PMDE, which should be in milliarcsec/yr.
;
; BUGS:
;   I have not handled the case where proper motion will move stars
;   over the poles, i.e. to DEC > 90 deg. ???
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written D. Schlegel, 31 December 2002, Princeton
;-
;------------------------------------------------------------------------------
pro tycho_epoch, epoch, tycdat

   nyear = epoch - 2000.d0
   tycdat.RAmdeg = $
    tycdat.RAmdeg + nyear * tycdat.pmRA / (1000.D * 3600.D * cos(tycdat.DEmdeg))
   tycdat.DEmdeg = $
    tycdat.DEmdeg + nyear * tycdat.pmDE / (1000.D * 3600.D)

   ; Rectify RA to be in domain [0,360).
   i = where(tycdat.RAmdeg LT 0.0, ct)
   if (ct GT 0) then tycdat[i].RAmdeg = tycdat[i].RAmdeg + 360.D
   i = where(tycdat.RAmdeg GE 360.0, ct)
   if (ct GT 0) then tycdat[i].RAmdeg = tycdat[i].RAmdeg - 360.D

end
;------------------------------------------------------------------------------
