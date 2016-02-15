PRO SELECT_O, selections, iselected, comments, command_line, only_one, $
              COUNT = n_select
;+
; NAME:
;	SELECT_O
; PURPOSE:
;	Dumb-terminal routine to let a user interactively select from a list
; EXPLANATION: 
;	This is the non-widget version of SCREEN_SELECT
;
; CALLING SEQUENCE:
;	select_o, selections, iselected, comments, command_line, only_one, $
;                            [ COUNT = ]
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
; OUTPUTS:
;	iselected - list of indices in selections giving the selected
;		items.
; OPTIONAL OUTPUT KEYWORD:
;       COUNT - Integer scalar giving the number of selections
; SIDE EFFECTS:
;	The obsolete system variable !err is set to the number of selections
; HISTORY:
;	version 1, D. Lindler  April 88.
;	modified to IDL V2 (from screen_select).  M. Greason, May 1990.
;	changed name from screen_select_o         W. Landsman January 1993
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Added COUNT keyword, deprecate !ERR   W. Landsman   March 2000
;-
;--------------------------------------------------------------------------
;			set defaults
;
IF n_params(0) LT 3 THEN comments=''
IF n_params(0) LT 4 THEN command_line=''
IF n_params(0) LT 5 THEN only_one=0
;
; 			initilization
;
n_select=0			;number of selections made
n=n_elements(selections)
nchar=max(strlen(selections))
help_avail=0			;help available flag
ncom=0				;Length of comments
IF n_elements(comments) EQ n THEN BEGIN
	ncom=max(strlen(comments))
	help_avail=1
ENDIF
IF n_elements(comments) EQ 1 THEN BEGIN				;scalar string
	IF strlen(strtrim(comments)) GT 0 THEN help_avail=1	;function name
ENDIF
question=0			;user asked for help
;
; 			determine screen format
;
inpos = 0
spos = 0
totchr = (nchar+ncom*question+3) < 79 		;total characters required
nx = 79 / totchr				;number in x direction
ny = (n+nx-1)/nx				;total number in y direction
screen=strarr(nx,ny)				;fill display string array.
selected = replicate(0B,nx,ny)			;vector of selected values
k = 0
FOR j = 0, ny-1 DO BEGIN
	FOR i = 0, nx-1 DO BEGIN
		IF (k LT n) THEN BEGIN
			st=selections[k]
			IF (ncom GT 0) AND (question EQ 1) THEN $
				st=' '+st+' '+comments[k]
			screen[i,j]=st
		ENDIF ELSE screen[i,j] = "  "
		k = k + 1
      	ENDFOR
ENDFOR
nlines = 22 < ny			;number of screen lines to display
nscr = (ny + nlines - 1) / nlines	;number of screens to display.
xpos = indgen(nx) * totchr + 1
ypos = indgen(nlines) + 1
;
; 			screen format init.
;
if !VERSION.OS EQ "vms" then cr = 13B else cr = 10B	;keystrokes
up = 128B
down = 130B
left = 129B
right = 131B
ix=0				;current position on screen
iy=0
iscr = 0
scr_curpos, 0, 0			;get scr_curpos compiled.
key = read_key(0)			;get read_key compiled.
refresh:
scr_attrib, 0       			;clear attributes
scr_other, '[?25l'			;disable visible cursor
scr_erase, 5         			;clear screen
;
; 			print initial contents of the screen
;
FOR j = 0, nlines-1 DO BEGIN
	k = j + spos
	FOR i = 0, nx-1 DO BEGIN 
		IF k LT ny THEN BEGIN
 		   	IF selected[i,k] THEN scr_attrib, 0, 1
			scr_curpos, ypos[j], xpos[i]
			print, screen[i,k]
			scr_attrib, 0
		ENDIF
	ENDFOR
ENDFOR
;
; 			print help lines in the text window.
;
message='Use arrow keys to move.  '
IF only_one THEN message=message+'<cr> or <space bar> to select' $
	    ELSE message=message+'<space bar> to select.  <cr> when done'
IF help_avail THEN message=message+'  ? for info.'
scr_attrib, 0, 1
scr_curpos, 22, 1
print, message
print, command_line
scr_attrib, 0
;
; 			loop until <cr>
;
key= 0B
WHILE key NE cr DO BEGIN
;
; 			high light current location
;
	scr_attrib, 0, 4
	scr_curpos, ypos[iy], xpos[ix]
	print, screen[ix,iy+spos]
	scr_attrib, 0
;
; 			process next key
;
	key = read_key(1)
;
; 			arrow key processing
;
	IF (key EQ up) OR (key EQ down) OR (key EQ right) OR (key EQ left) $
	    THEN BEGIN
	    scroll=0				;scroll flag
;
;				unhighlight the previous selection.
;
	    scr_curpos, ypos[iy], xpos[ix]	;remove cursor attrib.
	    IF selected[inpos] THEN BEGIN	;if selected, change attrib.
		scr_attrib, 0, 1
            ENDIF
	    print, screen[ix,iy+spos]
	    scr_attrib, 0
;
;				decode arrow key.
;
	    CASE key OF
		up: BEGIN
			IF iy GT 0 THEN BEGIN
				iy = iy - 1
				inpos = inpos - nx
		    	ENDIF ELSE BEGIN
				IF iscr GT 0 THEN BEGIN
				    iscr = iscr - 1
				    scroll = 1
				ENDIF
			ENDELSE
		    END
		down: BEGIN
			IF (iy LT (nlines-1)) AND ((iy+spos) LT (ny-1)) $
			THEN BEGIN
				inpos=inpos+nx
				iy=iy+1
			ENDIF ELSE BEGIN
				IF iscr LT (nscr - 1) THEN BEGIN
				    iscr = iscr + 1
				    scroll = 1
				ENDIF
			ENDELSE
		     END
		right: BEGIN
 			 IF ix LT (nx-1) THEN BEGIN
				ix=ix+1
				inpos=inpos+1
			 ENDIF
		       END
		left : BEGIN
			 IF ix GT 0 THEN BEGIN
				ix=ix-1
				inpos=inpos-1
			 ENDIF
		       END
	    ENDCASE
	    WHILE inpos GE n DO BEGIN		;prevent passing end-of-list
		inpos=inpos-1
		ix=ix-1
	    ENDWHILE
;
; 			do we need to scroll ?
;
		IF scroll THEN BEGIN
			iy = 0
			spos = iscr * nlines
			inpos = (nx * spos) + ix
			goto, refresh
		ENDIF
	ENDIF
;
; 			process other keys
;
	IF (only_one EQ 1) AND (key EQ cr) THEN key=' '	;select with cr also
	   if string(key) EQ ' ' THEN BEGIN
		    IF (NOT selected[inpos]) THEN BEGIN
			selected[inpos]=1B
			n_select=n_select+1
                        IF n_select EQ 1 THEN iselected = lonarr(1) + inpos $
                                         ELSE iselected = [iselected,inpos]
			IF only_one THEN BEGIN		;got our one selection?
                                iselected = iselected[0]
				goto,done
			ENDIF
		    ENDIF
		  ENDIF
	CASE strupcase(key) OF
	    cr  : goto,done
	    'R' : BEGIN
			selected[inpos]=0B
			n_select = (n_select - 1) > 0
		  END
	    '?' : BEGIN
		   IF (help_avail) THEN BEGIN
			IF (ncom EQ 0) THEN BEGIN  ;go get help text
			    scr_erase, 5
			    print, 'PLEASE WAIT....'
			    istat=execute(comments+',selections,comments')
			    ncom=strlen(comments[0])
			    print, 'FINISHED WITH HELP'
			ENDIF
			question=1
			goto,refresh
		   ENDIF
	          END
	    ELSE :
	ENDCASE
ENDWHILE
;
;			Finished.  Set !err to the number of items selected.
;
done:
scr_other, '[?25h'			;enable visible cursor
!err = n_select
scr_erase, 5
;
RETURN
END
