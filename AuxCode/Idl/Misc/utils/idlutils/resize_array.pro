;;
; 
; Copyright (C) 2006 Patricio Rojo
; 
; This program is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public License
; as published by the Free Software Foundation; either version 2
; of the License, or (at your option) any later version.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
;   MA  02110-1301, USA.
; 
;;

;+
; NAME:
;	RESIZE_ARRAY
;
; PURPOSE:
;	This function changes the size of an array maintaining the
;	original content.
;
; CATEGORY:
;	Pato's miscellaneous utilities.
;
; CALLING SEQUENCE:
;
;	RESIZE_ARRAY, Data, Newsize
;
; INPUTS:
;	Data:	   Data array of arbitrary dimensions and type.
;       Newsize:   New size of selected dimension.
;
; KEYWORD PARAMETERS:
;	DIMENSION: Dimension that want to be modified. Starts from 1.
;       NEWELEM:   Value of the new elements.
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	Algorithm description.
;
; EXAMPLE:
;	Describe example here
;
;               IDL> data = dblarr(2,3,5)
;		IDL> resize_array, data, 4, dim=2
;               IDL> help,data
;               ;DATA            DOUBLE    = Array[2, 4, 5]
;
; MODIFICATION HISTORY:
; 	Written by:	Pato Rojo, Cornell.  2006
;			pato@astro.cornell.edu
;-

pro resize_array, array, newsize, dimension=dimension, newelem=newelem

on_error, 2

if ~ keyword_set(dimension) then dim = 0 $
else dim = dimension-1

nd = size(array, /n_dim)
if dim gt nd-1 then message, string("Specified dimension (", dim+1, $
  ") is greater than array size (", nd, ")", format='(a,i0,a,i0,a)')

elemperdim = n_elements(array) / (size(array, /dim))[dim]
type = size(array, /type)

if ~ keyword_set(newelem) then begin
    if (type) eq 7 then newelem = '' $
    else newelem = 0
endif

n = newsize * elemperdim
nar = make_array(n, type=type) + newelem

odim = size(array, /dim)
ndim = odim
ndim[dim] = newsize
nar = reform(nar, ndim)

if n lt n_elements(array) then begin
    message, 'Reduction not implemented yet'
endif else begin
    case nd of
        1: nar[0:odim[0]-1]                                         = array
        2: nar[0:odim[0]-1, 0:odim[1]-1]                            = array
        3: nar[0:odim[0]-1, 0:odim[1]-1, 0:odim[2]-1]               = array
        4: nar[0:odim[0]-1, 0:odim[1]-1, 0:odim[2]-1, 0:odim[3]-1]  = array
        5: nar[0:odim[0]-1, 0:odim[1]-1, 0:odim[2]-1, 0:odim[3]-1, $
          0:odim[4]-1]                                              = array
        6: nar[0:odim[0]-1, 0:odim[1]-1, 0:odim[2]-1, 0:odim[3]-1, $
          0:odim[4]-1, 0:odim[5]-1]                                 = array
        7: nar[0:odim[0]-1, 0:odim[1]-1, 0:odim[2]-1, 0:odim[3]-1, $
          0:odim[4]-1, 0:odim[5]-1, 0:odim[6]-1]                    = array
        8: nar[0:odim[0]-1, 0:odim[1]-1, 0:odim[2]-1, 0:odim[3]-1, $
          0:odim[4]-1, 0:odim[5]-1, 0:odim[6]-1, 0:odim[7]-1]       = array
        else: message, string("Array can only have from 1 to 8 ", $
                              "dimensions (", nd, ")", format='(a,i0,a)')
    endcase
endelse

array = nar

end
