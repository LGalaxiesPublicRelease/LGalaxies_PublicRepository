; $Id: hist_2d.pro,v 1.14 2004/01/21 15:54:54 scottm Exp $
;
; Copyright (c) 1992-2004, Research Systems, Inc.  All rights reserved.
;   Unauthorized reproduction prohibited.
;+
; NAME:
;   HIST_2D
;
; PURPOSE:
;   Return the density function (histogram) of two variables.
;
; CATEGORY:
;   Image processing, statistics, probability.
;
; CALLING SEQUENCE:
;   Result = hist_2d(V1, V2)
; INPUTS:
;   V1 and V2 = arrays containing the variables.  May be any non-complex
;       numeric type.
;
; Keyword Inputs:
;       MIN1:   MIN1 is the minimum V1 value to consider. If this
;               keyword is not specified, then if the smallest value of
;               V1 is greater than zero, then MIN1=0 is used, otherwise
;               the smallest value of V1 is used.
;
;       MIN2:   MIN2 is the minimum V2 value to consider. If this
;               keyword is not specified, then if the smallest value of
;               V2 is greater than zero, then MIN2=0 is used, otherwise
;               the smallest value of V2 is used.
;
;       MAX1:   MAX1 is the maximum V1 value to consider. If this
;               keyword is not specified, then V1 is searched for
;               its largest value.
;
;       MAX2    MAX2 is the maximum V2 value to consider. If this
;               keyword is not specified, then V2 is searched for
;               its largest value.
;
;       BIN1    The size of each bin in the V1 direction (column
;               width).  If this keyword is not specified, the
;               size is set to 1.
;
;       BIN2    The size of each bin in the V2 direction (row
;               height).  If this keyword is not specified, the
;               size is set to 1.
;
; OUTPUTS:
;   The two dimensional density function of the two variables,
;   a longword array of dimensions (m1, m2), where:
;       m1 = Floor((max1-min1)/bin1) + 1
;      and  m2 = Floor((max2-min2)/bin2) + 1
;   and Result(i,j) is equal to the number of sumultaneous occurences
;   of an element of V1 falling in the ith bin, with the same element
;   of V2 falling in the jth bin.
;
; RESTRICTIONS:
;   Not usable with complex or string data.
;
; PROCEDURE:
;   Creates a combines array from the two variables, equal to the
;   linear subscript in the resulting 2D histogram, then applies
;   the standard histogram function.
;
; EXAMPLE:
;   Return the 2D histogram of two byte images:
;       R = HIST_2D(image1, image2)
;
;   Return the 2D histogram made from two floating point images
;   with range of -1 to +1, and with 101 (= 2/.02 + 1) bins:
;       f1 = RANDOMN(seed, 256, 256)
;       f2 = RANDOMN(seed, 256, 256)
;       R = HIST_2D(f1, f2, MIN1=-1, MIN2=-1, MAX1=1, MAX2=1, $
;           BIN1=.02, BIN2=.02)
;       TVSCL, R
;
; MODIFICATION HISTORY:
;   Written by:
;   DMS, Sept, 1992     Written
;   DMS, Oct, 1995      Added MIN, MAX, BIN keywords following
;               suggestion of Kevin Trupie, GSC, NASA/GSFC.
;   CT, RSI, May 2001: Corrected MIN, MAX keywords so that the out-of-range
;               values are ignored rather than truncated to be within range.
;               Allow input arrays with negative values.
;-
function hist_2d, im1, im2, $
    MIN1 = min1in, MIN2 = min2in, $
    MAX1 = max1in, MAX2 = max2in, $
    BIN1 = b1in, BIN2 = b2in

    COMPILE_OPT idl2
    ON_ERROR, 2

    ;Find extents of arrays.
    im1max = MAX(im1, MIN=im1min)
    im2max = MAX(im2, MIN=im2min)

    ;Supply default values for keywords.
    min1 = (N_ELEMENTS(min1in) gt 0) ? min1in : (0 < im1min)
    max1 = (N_ELEMENTS(max1in) gt 0) ? max1in : im1max
    min2 = (N_ELEMENTS(min2in) gt 0) ? min2in : (0 < im2min)
    max2 = (N_ELEMENTS(max2in) gt 0) ? max2in : im2max
    b1 = (N_ELEMENTS(b1in) gt 0) ? b1in : 1L
    b2 = (N_ELEMENTS(b2in) gt 0) ? b2in : 1L

    ;Get # of bins for each
    im1bins = FLOOR((max1-min1) / b1) + 1L
    im2bins = FLOOR((max2-min2) / b2) + 1L

    if (im1bins le 0) then MESSAGE, 'Illegal bin size for V1.'
    if (im2bins le 0) then MESSAGE, 'Illegal bin size for V2.'

    noMinTruncation = (min1 eq im1min) and (min2 eq im2min)
    noMaxTruncation = (im1max le max1) and (im2max le max2)
    binSizeOne = (b1 eq 1) and (b2 eq 1)

    ; Combine im1 and im2 into a single array.
    ; Use im2 values as the "row" indices, use im1 values as the "columns".
    if (min1 eq 0) and (min2 eq 0) and $
        (noMinTruncation and noMaxTruncation and binSizeOne) then begin
        ; Fast case without scaling?
        h = im1bins*LONG(im2) + LONG(im1)
    endif else begin
        im1tmp = im1
        im2tmp = im2

        ; Only do the data conversions that are necessary.
        if (min1 ne 0) then im1tmp = TEMPORARY(im1tmp) - min1
        if (min2 ne 0) then im2tmp = TEMPORARY(im2tmp) - min2
        if (b1 ne 1) then im1tmp = TEMPORARY(im1tmp)/b1
        if (b2 ne 1) then im2tmp = TEMPORARY(im2tmp)/b2
        h = im1bins*LONG(TEMPORARY(im2tmp)) + LONG(TEMPORARY(im1tmp))

        ; Construct an array of out-of-range (0) and in-range (1) values.
        in_range = 1
        if (noMinTruncation eq 0) then $ ; set lt min to zero
            in_range = (im1 ge min1) and (im2 ge min2)
        if (noMaxTruncation eq 0) then $ ; set gt max to zero
            in_range = TEMPORARY(in_range) and (im1 le max1) and (im2 le max2)
        ; Set values that are out of range to -1
        h = (TEMPORARY(h) + 1L)*TEMPORARY(in_range) - 1L
    endelse

    h = HISTOGRAM(h, MIN=0, MAX=im1bins*im2bins-1)  ;Get the 1D histogram
    return, REFORM(h, im1bins, im2bins, /OVERWRITE) ;and make it 2D
end
