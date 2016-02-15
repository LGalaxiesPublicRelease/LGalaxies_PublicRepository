PRO SCREEN_SELECT, selections, iselected, comments, command_line, only_one, $
                COUNT  = count      
;+
; NAME:
;	SCREEN_SELECT
; PURPOSE:
;	Allow a user to make an interactive screen selection from a list
; EXPLANATION:
;	This procedure determines whether to use the dumb terminal version,  
;	or the widget version by examining the !D.NAME system variable.
;
; CALLING SEQUENCE:
;	screen_select, selections, iselected, comments, command_line, only_one
;
; INPUTS:
;	selections - string array giving list of items that can be
;		selected.
;
; OPTIONAL INPUTS:
;	comments - comments which can be requested for each item in
;		array selections.  It can be:
;			string array - same length as array selections.
;			null string - no comments available
;			scalar string - name of a procedure which will
;				return comments.  It will take selections
;				as its first argument and return comments
;				as its second argument.
;	command_line - optional command line to be placed at the bottom
;		of the screen.  It is usually used to specify what the
;		user is selecting.
;	only_one - integer flag. If set to 1 then the user can only select
;		one item.  The routine returns immediately after the first
;		selection is made.
;
; OUTPUTS:
;	iselected - list of indices in selections giving the selected
;		items.
;
; OPTIONAL OUTPUT KEYWORD:
;       COUNT - Integer scalar giving the number of selections made
;
; SIDE EFFECTS:
;	The obsolete system variable !err is set to the number of selections
;
; PROCEDURE:
;	The actual processing is farmed out to different procedures depending
;	on the terminal type.    
;
;	Widget Terminal   ==>  SELECT_W.PRO
;	VT100 Terminal  ==>    SELECT_O.PRO
; HISTORY:
;	Written by M. Greason, STX, May 1990.
;       Added widget support    W. Landsman           January, 1992
;	Remove X window but no widget option         November, 1994
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Added COUNT keyword, deprecate !ERR   W. Landsman   March 2000
;-
;--------------------------------------------------------------------------
;			Set defaults.
;
if N_params() LT 3 then comments = ''
if N_params() LT 4 then command_line = ''
if N_params() LT 5 then only_one = 0
;
;			If no selection array or output array has been 
;			supplied, set iselected to -1 and !err to 0.
;			Then quit.
;
 count = 0
 if N_params() LT 2 then begin
	!ERR = 0
	iselected = -1
 endif else begin
;
;			Determine which procedure to farm the work out
;			to depending upon the contents of !d.name and !D.flags
;
   if (!D.FLAGS and 65536) EQ 65536 then $       ;Widgets?

      select_w, selections,iselected, comments, command_line, only_one,  $
                Count = count else  $                                       
                                         ;Dumb Terminal?
      select_o, selections, iselected, comments, command_line, only_one, $
                Count = count

 endelse
;
 return
 end
