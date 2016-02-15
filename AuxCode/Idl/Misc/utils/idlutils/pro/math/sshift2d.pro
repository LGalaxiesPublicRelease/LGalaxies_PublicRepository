;+
; NAME:
;    sshift2d
; PURPOSE: (one line)
;    Shift a 2-D array using a damped sinc function for the fractional part.
; DESCRIPTION:
;
; CATEGORY:
;    Mathematical
; CALLING SEQUENCE:
;    result = sshift2d( array, shiftvec )
; INPUTS:
;    array    : Array to be shifted.
;    shiftvec : Vector of two elements: [ xshift, yshift ].
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;    The shifted array is returned as the function value.
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;    February, 1993:
;    Copied from "sincshift.pro" written by John Spencer, Lowell Observatory.
;    Very minor modifications by Doug Loucks, Lowell Observatory.
;-
;******************************************************************************
; Begin support functions.
; ------------------------------------------------------------------------------
; Function Padarray
; Returns original 2-D array padded by margin replication
; on all sides.  x margin = pad( 0 ), y margin=pad( 1 )
; ------------------------------------------------------------------------------
FUNCTION Padarray, array, pad
sizearr = SIZE( array )
xs = sizearr[ 1 ]
ys = sizearr[ 2 ]

parray = FLTARR( xs+2*pad[ 0 ], ys+2*pad[ 1 ] )
; Middle:
parray[ pad[ 0 ], pad[ 1 ] ] = array
; low-x side:
parray[ 0 : pad[0]-1, pad[1] : pad[1]+ys-1 ] = REBIN( array[0,*], pad[0], ys )
; high-x side:
parray[ pad[0]+xs : *, pad[1] : pad[1]+ys-1 ]= REBIN( array[xs-1,*], pad[0],ys )
; low-y side:
parray[ *, 0 : pad[1]-1 ] = REBIN( parray[ *, pad[1] ], xs+2*pad[0], pad[1] )
; high-y side:
parray[*, pad[1]+ys : * ]= REBIN( parray[ *, pad[1]+ys-1 ], xs+2*pad[0],pad[1])

RETURN, parray
END

; ------------------------------------------------------------------------------
; Function Shiftrep
; Shifts an image by integer numbers of pixels, replicating border pixels
; instead of wrapping around
; ------------------------------------------------------------------------------
FUNCTION Shiftrep, image, xshift, yshift
imsize = SIZE( image )
xsize = imsize[ 1 ]
ysize = imsize[ 2 ]

; Nint is the "nearest integer" function found in the Astron_misc
; library,
ix = Nint( xshift )
iy = Nint( yshift )

temp = SHIFT( image, ix, iy )
IF ix GT 0 THEN temp[ 0 : ix-1, * ] = REBIN( temp[ ix, * ], ix, ysize )
IF ix LT 0 THEN temp[xsize+ix : *, *] = REBIN(temp[ xsize+ix-1, * ], -ix, ysize)
IF iy GT 0 THEN temp[ *, 0 : iy-1 ] = REBIN( temp[ *, iy ], xsize, iy )
IF iy LT 0 THEN temp[*, ysize+iy : *] = REBIN(temp[ *, ysize+iy-1 ], xsize, -iy)
RETURN, temp
END
; End support functions.
;******************************************************************************

; ------------------------------------------------------------------------------
; Function sshift2d
; Shifts a 2-D array by shiftvec( 0 ) along x and shiftvec( 1 ) along y,
; replicating the margins.
; ------------------------------------------------------------------------------
FUNCTION sshift2d, array, shiftvec
sizearr = SIZE( array )
xs = sizearr[ 1 ]
ys = sizearr[ 2 ]

; Separate and do integer shift first:
ishift = Nint( shiftvec )
fshift = shiftvec - ishift
sarray = Shiftrep( array, ishift[ 0 ], ishift[ 1 ] )

; Return if there's no fractional shift:
IF fshift[ 0 ] EQ 0 AND fshift[ 1 ] EQ 0 THEN RETURN, sarray

dampfac = 3.25
sincrad = 10
sincsize = 2 * sincrad + 1

; Pad the array (replicating margins) in preparation for convolution:
sarray = Padarray( sarray, [sincrad,sincrad] )

; Generate the x and y sinc functions:
sinc = FLTARR( sincsize, 2 )
kernel = FLTARR( sincsize, sincsize )

FOR iindex = 0, 1 DO BEGIN
   IF fshift[ iindex ] NE 0 THEN BEGIN
      s = FINDGEN( 2*sincrad+1 ) - sincrad + fshift[ iindex ]
      sinc[ *, iindex ] = EXP( -( s/dampfac )^2 ) * SIN( !PI*s )/( !PI*s )
   ENDIF ELSE BEGIN
      sinc[ *, iindex ] = 0.0
      sinc[ sincrad, iindex ] = 1.0
   ENDELSE
ENDFOR

; Multiply the sinc functions to generate 2-D kernel
FOR i = 0, sincsize-1 DO BEGIN
   kernel[ *, i ] = sinc[ *, 0 ] * sinc[ i, 1 ]
ENDFOR

sarray = CONVOL( sarray, kernel, CENTER=1 )

;Trim back to original size:
sarray = sarray[ sincrad : sincrad+xs-1, sincrad : sincrad+ys-1 ]

RETURN, sarray
END
