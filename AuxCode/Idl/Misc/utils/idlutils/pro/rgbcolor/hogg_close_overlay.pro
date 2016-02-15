;+
; NAME:
;  hogg_close_overlay
; PURPOSE:
;  close Z buffer and output contents to an image for overlaying
; INPUTS:
;  naxis        - [naxis1,naxis2] size of image
; OUTPUTS:
;  overlay      - image containing the overlay!
; BUGS:
;  - Code not written.
;  - Header not written.
;  - Ought to restore saved !P,!X,!Y variables.
;-
pro hogg_close_overlay, naxis,overlay
snapshot= tvrd()
tvlct, r,g,b,/get
overlay= (255-byte(rebin((r[snapshot]),naxis[0],naxis[1])))
return
end
