PRO select_w_event, event
;
;This procedure is the event handler for the XMENU widget below
COMMON select_w, val, exclusive

WIDGET_CONTROL, event.id, GET_VALUE = value, GET_UVALUE = i

;start:

; Get the selections
if (event.select EQ 1) then val = [val,i] $
                       else val = val[ where( val NE i) ]

if (value EQ 'DONE') or (exclusive) then begin  
              good  = where( val GE 0, nsel )
              if (nsel GT 0) THEN val = val[good] 
              widget_control, event.top, /DESTROY
 END
END

PRO select_w, items, iselected, comments, command_line, only_one, $
	Count = count, GROUP_LEADER=GROUP, selectin = selectin
;+
; NAME:
;	SELECT_W    
; PURPOSE:
;	Create a non-exclusive widget menu of items
; EXPLANATION:
;	More than one item may be selected or 'de-selected'.   
;	Normally called by SCREEN_SELECT
;
; CALLING SEQUENCE:
;	SELECT_W, items ,iselected, [ comments, command_line, only_one ]
;
; INPUTS:
;	items - string array giving list of items that can be
;		selected.
;
; OPTIONAL INPUTS:
;	comments - comments which can be requested for each item in
;		array selections.    NOT YET IMPLEMENTED
;	command_line - optional command line to be placed at the bottom
;		of the screen.  It is usually used to specify what the
;		user is selecting.
;	only_one - integer flag. If set to 1 then the user can only select
;		one item.  The routine returns immediately after the first
;		selection is made.
; OPTIONAL KEYWORD INPUT
;       SELECTIN - vector of items to be pre-selected upon input (not used for
;               only_one option)
;
; OUTPUT:
;	iselected - list of indices in selections giving the selected
;		items.
;
; OPTIONAL OUTPUT KEYWORD:
;       COUNT  - Integer scalar giving the number of items selected
; COMMON BLOCKS:
;	SELECT_W - Used to communicate with the SELECT_W_EVENT procedure 
;
; MODIFICATION HISTORY:
;	Written, K. Venkatakrishna & W. Landsman, Hughes/STX    January, 1992
;	Widgets made MODAL.  M. Greason, Hughes STX, 15 July 1992.
;       Changed handling of MODAL keyword for V5.0   W.Thompson  September 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Added selectin keyword  D. Lindler 01/12/99 
;-
;
 On_error,2
 common select_w, val, exclusive

 if N_params() LT 5 then exclusive = 0 else exclusive = only_one

 val = -1

 if N_params() LT 4 then command_line = $ 
' Select by pressing the left mouse button once; To de-select press twice; finally QUIT'
    
        MODAL = N_ELEMENTS(GROUP) GE 1
        base = WIDGET_BASE( TITLE = command_line, /COLUMN, MODAL=MODAL, $
                GROUP_LEADER=GROUP)
 if only_one then $
       XMENU, items, base, COLUMN=8  $
    else begin 
       donebut = WIDGET_BUTTON( base, VALUE = "DONE", UVALUE = -1) 
       XMENU, items, base, /NONEXCLUSIVE, COLUMN=8, buttons=buttons
       if n_elements(selectin) gt 0 then begin
                for i=0,n_elements(selectin)-1 do $
                        widget_control,buttons(selectin[i]),set_button=1
                val = [-1,selectin]
       endif
 endelse

; Realize the widgets:
 WIDGET_CONTROL, base, /REALIZE

; Hand off to the XMANAGER, i.e.,event-handler,:
  XMANAGER, 'select_w', base, GROUP_LEADER = GROUP
 if val[0] NE -1 then iselected = val
 count = N_elements( iselected)
 !ERR = count

 return
 end

