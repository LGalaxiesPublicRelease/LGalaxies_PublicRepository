pro RESET_RDPLOT
;
;   If the user crashes out of the RDPLOT program, they can call this procedure
; to reset the graphics device functions to default values.
;
device, /CURSOR_CROSSHAIR, SET_GRAPHICS_FUNCTION=3, BYPASS_TRANSLATION=0
return
end



pro RDPLOT, x, y, WaitFlag, DATA=Data, DEVICE=Device, NORMAL=Normal, $
   NOWAIT=NoWait, WAIT=Wait, DOWN=Down, CHANGE=Change, Err=Err, $
   PRINT=Print, XTITLE=XTitle,YTITLE=YTitle, XVALUES=XValues,YVALUES=YValues, $
   FULLCURSOR=FullCursor, NOCLIP=NoClip, LINESTYLE=Linestyle, THICK=Thick, $
   COLOR=Color, CROSS=Cross
   
;*******************************************************************************
;+
; NAME:
;   RDPLOT
;
; PURPOSE:
;   Like CURSOR but with a full-screen cursor and continuous readout option
; EXPLANATION:
;   This program is designed to essentially mimic the IDL CURSOR command,
;   but with the additional options of continuously printing out the data
;   values of the cursor's position, and using a full-screen cursor rather 
;   than a small cross cursor.  The Full screen cursor uses OPLOT and 
;   X-windows graphics masking to emulate the cursor.
;      One difference is that IF the PRINT keyword is set but the DOWN, WAIT,
;   or CHANGE keywords are not set, then the leftmost mouse button will 
;   print a "newline" line-feed, but not exit.
;
; CALLING SEQUENCE:
;   RDPLOT, [X, Y, WaitFlag], [/DATA, /DEVICE, /NORMAL,
;      /NOWAIT, /WAIT, /DOWN, /CHANGE, ERR=,
;      PRINT=, XTITLE=, YTITLE=, XVALUES=, YVALUES=,
;      /FULLCURSOR, /NOCLIP, LINESTYLE=, THICK=, COLOR=, /CROSS]
;
; REQUIRED INPUTS:
;   None.
;
; OPTIONAL INPUTS: 
;   WAITFLAG = Uses the same table as the intrinsic CURSOR command, But note
;	that unlike the CURSOR command, there is no UP keyword.
;		WaitFlag=0 sets the NOWAIT keyword
;		WaitFlag=1 sets the WAIT keyword {default}
;		WaitFlag=2 sets the CHANGE keyword
;		WaitFlag=3 sets the DOWN keyword
;
; OPTIONAL OUTPUTS:
;    X - a named variable to receive the final cursor X position, scalar
;    Y - a named variable to receive the final cursor Y position, scalar
; OPTIONAL KEYWORD INPUT PARAMETERS:
;   /DATA = Data coordinates are displayed and returned.
;   /DEVICE = device coordinates are displayed and returned.
;   /NORMAL = normal coordinates are displayed and returned.
;          Default is to use DATA coordinates if available (see notes).
;   /NOWAIT = if non-zero the routine will immediately return the cursor's
;      present position.
;   WAIT = if non-zero will wait for a mouse key click before returning.  If
;      cursor key is already down, then procedure immediately exits.
;   DOWN = equivalent to WAIT *except* that if the mouse key is already down
;      when the procedure is called, the procedure will wait until the mouse
;      key is clicked down again.
;   CHANGE = returns when the mouse is moved OR a key is clicked up or down.
;   PRINT = if non-zero will continuously print out (at the terminal) the data 
;      values of the cursor's position.  If PRINT>1, program will printout a 
;      brief header describing the mouse button functions.  However, note that 
;      the button functions are overridden if any of the DOWN, WAIT, mouse
;      or CHANGE values are non-zero.
;   XTITLE = label used to describe the values of the abscissa if PRINT>0.
;   YTITLE = label used to describe the values of the ordinate if PRINT>0.
;   XVALUES = a vector corresponding to the values to be printed when the
;	PRINT keyword is set.  This allows the user the option of printing
;	out other values rather than the default X coordinate position of
;	the cursor.  E.g., if XVALUES is a string vector of dates such as
;	['May 1', 'May 2', ...], then those dates will be printed rather than
;	the X value of the cursor's position: if X=1 then 'May 2' would be
;	printed, etc.  This requires that the values of the X coordinate read
;	by the cursor must be positive (can't access negative elements).
;       If XVALUES=-1, then NO values for X will be printed.
;   YVALUES = analagous to the XVALUES keyword.
;   FULLCURSOR = if non-zero default cursor is blanked out and full-screen 
;      (or full plot window, depending on the value of NOCLIP) lines are
;      drawn; their intersecton is centered on the cursor position.
;   NOCLIP = if non-zero will make a full-screen cursor, otherwise it will
;      default to the value in !P.NOCLIP.
;   LINESTYLE = style of line that makes the full-screen cursor.
;   THICK = thickness of the line that makes the full-screen cursor.
;   COLOR = color of the full-screen cursor.
;   CROSS = if non-zero will show the regular cross AND full screen cursors.
;
; OPTIONAL KEYWORD OUTPUT PARAMETER:
;   ERR = returns the most recent value of the !mouse.button value.
;
; NOTES:
;   Note that this procedure does not allow the "UP" keyword/flag...which 
;   doesn't seem to work too well in the origianl CURSOR version anyway.
;
;   If a data coordinate system has not been established, then RDPLOT will
;   create one identical to the device coordinate system.   Note that this
;   kluge is required even if the user specified /NORMAL coordinates, since
;   RDPLOT makes use of the OPLOT procedure.  This new data coordinate system
;   is effectively "erased" (!X.CRange and !Y.CRange are both set to zero)
;   upon exit of the routine so as to not change the plot status from the
;   user's point of view.
;
;   Only tested on X-windows systems.  If this program is interrupted, the
;   graphics function might be left in a non-standard state; in that case,
;   run the program RESET_RDPLOT to return the standard graphics functions,
;   or type the command:   DEVICE, /CURSOR_CROSS, SET_GRAPHICS=3, BYPASS=0
;
; BUGS:
;   It is assumed that the current background of the plot is correctly
;   defined by the value in !P.Background.  Otherwise, the color of the
;   long cursor probably will not be correct.  Sometimes the color doesn't
;   work anyway, and I'm not sure why.
;
;   There may be some cases (e.g., when THICK>1 and NOCLIP=0) when the
;   full-screen cursor is not correctly erased, leaving "ghost images" on the
;   plot.  It just seems that the screen updates get slow or the positions
;   ambiguous with a thick line and the cursor off the plot.
;
; PROCEDURE:
;   Basically is a bells-n-whistles version of the CURSOR procedure.  All
;   the details are covered in the above discussion of the keywords.
;
; EXAMPLE (a silly, but informative one):
;   Months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', $
;             'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
;   plot, indgen(12), xrange=[-5, 15]
;   rdplot, /FULL, /PRINT, XTITLE='Month: ', YTITLE='Y-value per month = ', $
;      xvalues=Months
;
; MODIFICATION HISTORY:
;   Written (originally named CURFULL) by J.Wm.Parker  1993 Nov 22 
;   Created data coordinates if not already present, W. Landsman Nov. 93
;   Added continuous printout of data values, COLOR and FULLCURSOR keywords
;      (so that default is that it acts just like the cursor command).
;      Changed name from CURFULL to RDPLOT.   J.Wm.Parker  1994 Apr 20
;   Modified (with some translation table assistance from the IDL support 
;      group) to correctly plot the crosshair with the desired IDL 
;      color using the device's translation table to determine the XOR 
;      function and using the BYPASS function.  Added the RESET_RDPLOT
;      procedure to cleanup crashes that might occur while running
;      RDPLOT.  Other minor changes/bug fixes.  J.Wm.Parker  1994 May 21
;   Modified DOWN, WAIT, CHANGE functions to behave more similar to the
;      generic CURSOR procedure.   J.Wm.Parker  1995 April 24
;   Added XVALUES, YVALUES keywords and cleanup.   J.Wm.Parker  1995 April 24
;   Convert to IDL V5.0,  W. Landsman    July 1998
;   Change !D.NCOLORS to !D.TABLE_SIZE for 24 bit displays W. Landsman May 2000
;   Skip translation table for TrueColor visuals   W. Landsman  March 2001
;-
;*******************************************************************************
On_error,2

;;;
;   If the device does not support windows, then this program can not be used.
;

if ((!D.Flags and 256) ne 256) then message, $
  'ERROR - Current graphics device ' + !D.NAME + ' does not support windows'

;;;
;   Keywords, keywords.
;
if (N_Params() eq 3) then begin
   case WaitFlag of
      0 : NoWait = 1
      1 : Wait = 1
      2 : Change = 1
      3 : Down = 1
      else : Wait = 1
   endcase
endif

NoWait = keyword_set(NoWait)
Wait = keyword_set(Wait)
Down = keyword_set(Down) or Wait
Change = keyword_set(Change)
FullCursor = keyword_set(FullCursor)


;;;
;   If plotting coordinates are not already established, and the NORMAL keyword
; is not set, then use device coordinates.
;   Note that even if this procedure was called with the DATA keyword set, that
; the DEVICE keyword will always take precedence over the DATA keyword in the
; cursor command.  However, if the NORMAL and DEVICE keywords are both set,
; then very strange values are returned.
;
UndefinedPlot = ((!X.CRange[0] eq 0) and (!X.CRange[1] eq 0))
if UndefinedPlot then plot, [0,!D.X_Size], [0,!D.Y_Size], /NODATA, $
   XSTYLE=5, YSTYLE=5, XMARGIN=[0,0], YMARGIN=[0,0], /NOERASE


;;;
;   Initialize the !mouse.button variable.  The value of !mouse.button 
; corresponds to the BYTE  value of the buttons on the mouse from left to right,
; lowest bit first.  So, the left button gives !mouse.button = 1, next button 
; gives !mouse.button = 2, then 4.
;  Read in the cursor with no wait.  If the user does not want to wait, or if 
; the DOWN or WAIT keywords are set AND the mouse key is depressed, then we're
; done (I hate GOTO's, but it is appropriate here).
;
!mouse.button = 0
cursor, X, Y, /NOWAIT, DATA=Data, DEVICE=Device, NORMAL=Normal
if (keyword_set(NoWait) or (Wait and (!mouse.button gt 0))) then $
            goto, LABEL_DONE


;;;
;   PRINTOUT SETUP SECTION ==================================================
;;;

;;;
;   Is the PRINT keyword set?  Then we have a lot of things to set up.  First,
; set up carriage return and line feed variables for the formatted printout,
; and define the titles for the printed values.
;
if keyword_set(Print) then begin 
   CR = string("15b)
   LF = string("12b)
   if not(keyword_set(XTitle)) then XTitle = "X = "
   if not(keyword_set(YTitle)) then YTitle = "Y = "
   Blanks  = "                    "

;;;
;   Now, if the XValues and/or YValues keywords are set, then deal with them.
; Also, we may want to suppress the printing of the X or Y values (e.g.,
; XValues=-1 or YValues=-1 sets the XhowX and ShowY variables).
;
   ShowX = 1
   XVfmt = "(A13)"
   UseXV = keyword_set(XValues)
   if UseXV then begin
      XVSt = string(XValues)
      XVtop = n_elements(XValues) - 1
      XVfmt = "(A" + strtrim(max(strlen(XVst))+3,2) + ")"
      if ((XVtop eq 0) and (strtrim(XVSt[0],2) eq '-1')) then ShowX = 0
   endif else XVfmt = "(A13)"
   if not(ShowX) then XTitle = ''


   ShowY = 1
   UseYV = keyword_set(YValues)
   if UseYV then begin
      YVSt = string(YValues)
      YVtop = n_elements(YValues) - 1
      YVfmt = "(A" + strtrim(max(strlen(YVst)),2) + ")"
      if ((YVtop eq 0) and (strtrim(YVSt[0],2) eq '-1')) then ShowY = 0
   endif else YVfmt = "(A13)"
   if not(ShowY) then YTitle = ''

;;;
;   If Print>1, then printout the informative header, which will vary depending
; on the values of the DOWN and CHANGE keywords.
;
   if (Print gt 1) then begin
      print, ' '
      if Change then begin
         print, " Hit any mouse button or move the mouse to exit."
      endif else begin
         if Down then begin
            print, " Hit any mouse button to exit."
         endif else begin
            print, ' Mouse Button:    LEFT         MIDDLE        RIGHT'
            print, ' Result Action:   New Line     Exit          Exit'
         endelse
      endelse
      print, ' '
   endif

endif else Print = 0



;;;
;   FULL-SCREEN CURSOR SETUP SECTION =======================================
;;;

;;;;
; If using the full-screen cursor:
;   Determine the data range for the full screen.
;   Blank out the regular cross cursor if the CROSS keyword is not set.
;   Set up the linestyle, thickness, clipping, and color parameters for the 
; oplot commands.
;   Set up the graphics to be XOR with the overplotted crosshair, and figure
; out the color to use for plotting the crosshair {details below}.
;
if FullCursor then begin
   Yfull = convert_coord([0.0,1.0], [0.0,1.0], /NORMAL, /TO_DATA)
   Xfull = Yfull[0,*]
   Yfull = Yfull[1,*]

   device, GET_GRAPHICS=OldGraphics, SET_GRAPHICS=6
   if not(keyword_set(Cross)) then device, CURSOR_IMAGE=intarr(16)

   if not(keyword_set(Linestyle)) then Linestyle = 0
   if not(keyword_set(Thick)) then Thick = 1
   NoClip = keyword_set(NoClip)

;;;
;   I think the best way to make the fullscreen cursor work is to use the XOR
; graphics function - overplotting a line will XOR with the data already on
; the screen, then overplotting the same line again will XOR again, effectively
; erasing the line and returning the device to it original state/appearance.
;    But first, let me present a quick primer on plotting colors in IDL and the 
; related color tables and translation table:
;   Normally, when a color N (a number between 0 and 255 which refers to a
; particular color in the currently loaded IDL color table) is used in one of
; the plotting or tv commands, the value that is actually sent to the display is
; the value in the N-th bin of the translation table.  E.g., if the background
; color is 0, then the actual (device) color value of the background is the
; value in the zeroth bin of the translation table.  Similarly, if the user
; wants to plot the color defined by number 147 in the IDL color table, the
; actual (device) color value of that color is the value in the 147th bin
; of the translation table.
;  So in the following example, let's pretend we have the following situation:
;   IDL> PRINT, !D.N_Colors
;            222
;   IDL> PRINT, !P.Background
;              0
;   IDL> DEVICE, TRANSLATION=TTab
;   IDL> PRINT, TTab(0)
;             34
;   IDL> PRINT, TTab(147)
;            181
;   When we set DEVICE,SET_GRAPHICS=6, and do an overplot, it performs an XOR
; function between the overplot's translated color value and the background's
; translated color value.
;   If we want the resulting color to be the IDL color 147, then we have to 
; overplot with the color whose translated color value XOR'ed with the 
; background's translated color value (34) will equal 181, which is the 
; translated color value of the desired IDL color 147.
;
; Symbolically:
; *  TTab(Desired Color) = TTab(OPLOT color) XOR TTab(Background)
; *  OPLOT Color = where( TTab eq (TTab(Desired Color) XOR TTab(Background)) )
;
; Numerically {using the above example}:
; *  OPLOT Color = where( TTab eq (TTab(147) XOR TTab(0)) )
; *  OPLOT Color = where( TTab eq (181 XOR 34) )
; *  OPLOT Color = where( TTab eq 151 )
;
;   Fine.
;   HOWEVER...since the translation table often does NOT contain the full range
; of possible numbers (e.g., 0 to 255), the result of the XOR function between 
; the background and the oplot color may be a value that does NOT appear in the 
; translation table.  This is particularly a problem for colors near the bottom
; of the translation table where the result of the XOR function may be less than
; the lowest value in TTab.
;   To fix this problem, I bypass the translation table, and directly send the
; device color (e.g., the value 151 in the above example) to the OPLOT command.
;   There is still some bug here - sometimes the color still isn't right.  I'll
; have to talk to the IDL support people about this {as soon as our support
; license is renewed!}
;
   if ((size(Color))[1] eq 0) then Color = !D.Table_size - 1   ;  if undefined
   device, get_visual_name=visualName
   if visualName NE 'TrueColor' then begin
    device, TRANSLATION=TTab, BYPASS_TRANSLATION=1
    DevColor = TTab[Color < (!D.Table_size - 1)]
    DevBack  = TTab[!P.Background]
    OColor = DevColor xor DevBack
   endif else Ocolor = color
endif 



;;;
;   FINALLY...THE PLOT READING SECTION  ====================================
;;;

;;;
;   If the cursor is beyond the boundaries of the window (device coordinates of
; X=-1 and Y=-1), then wait until the cursor is moved into the window.
;
cursor, X, Y, /NOWAIT, /DEVICE
if ((X lt 0) or (Y lt 0)) then cursor, X, Y, /CHANGE


;;;
;   Begin the loop that will repeat until a button is clicked (or a change if
; that is what the user wanted).   Err0 is used to keep track if the procedure
; was entered with a key already down, then it will be non-zero until that
; key has been released, at which point it will be permanantly set to zero.
;   Wait for a change (movement or key click).  Delete the old lines, and
; if we don't exit the loop, repeat and draw new lines.
;
cursor, X, Y, /NOWAIT, DATA=Data, DEVICE=Device, NORMAL=Normal
Err0 = !mouse.button


repeat begin	; here we go!

;;;
;   If doing a full-screen cursor, overplot two full-screen lines intersecting 
; at that position.
;
   if FullCursor then begin
      XY = convert_coord(X,Y, DATA=Data,DEVICE=Device,NORMAL=Normal, /TO_DATA)
      Xdata = XY[0] * [1.0,1.0]
      Ydata = XY[1] * [1.0,1.0]
      oplot, Xdata,Yfull, LINE=Linestyle,THICK=Thick,NOCLIP=NoClip,COLOR=OColor
      oplot, Xfull,Ydata, LINE=Linestyle,THICK=Thick,NOCLIP=NoClip,COLOR=OColor
   endif

;;;
;   If printing out data values, do so.
;   !mouse.button=1 is the signal for a new line.
;
   if (Print gt 0) then begin

      if ShowX then begin
         if UseXV then Xst = XVSt[(X+0.5) > 0 < XVtop] else Xst = strtrim(X,2)
         XSt = XTitle + string(Xst + Blanks, FORMAT=XVfmt)
      endif else Xst = ''
      if ShowY then begin
         if UseYV then Yst = YVSt[(Y+0.5) > 0 < YVtop] else Yst = strtrim(Y,2)
         YSt = YTitle + string(Yst + Blanks, FORMAT=YVfmt)
      endif else Yst = ''

      print, CR, Xst, Yst, format='($,3A)'

      if ((!mouse.button eq 1) and not(Down or Change)) then begin  ;  new line?
         print, LF, format="($,a)"
         while (!mouse.button gt 0) do begin  ; if button is held down, don't print
            wait, 0.1
            cursor, X, Y, /NOWAIT
         endwhile
      endif
   endif

   Err0 = Err0 < !mouse.button

;;;
;   Check to see that the cursor's current position is really the last measured 
; position (the mouse could have moved during a delay in the last section).  If
; so, then go on.  If not, then wait for some change in the mouse's status 
; before going on.
;  In either case, once we are going on, then if doing a full-screen cursor, 
; overplot the previous lines {the XOR graphics function will return the plot
; to its original appearance}.  Repeat until exit signal.
;
   cursor, XX, YY, /NOWAIT, DATA=Data, DEVICE=Device, NORMAL=Normal
   if ((XX eq X) and (YY eq Y)) then $
      cursor, XX, YY, /CHANGE, DATA=Data, DEVICE=Device, NORMAL=Normal

   if FullCursor then begin
      oplot, Xdata,Yfull, LINE=Linestyle,THICK=Thick,NOCLIP=NoClip,COLOR=OColor
      oplot, Xfull,Ydata, LINE=Linestyle,THICK=Thick,NOCLIP=NoClip,COLOR=OColor
   endif

;;;
;   Load the new XX and YY values into the X and Y variables.  Return or exit.
;
   X = XX
   Y = YY
   Err = !mouse.button

   if Down then ExitFlag = (!mouse.button gt 0) and (Err0 eq 0) else $
                ExitFlag = (!mouse.button gt 1) or Change

endrep until ExitFlag

if (Print gt 0) then print, ""



LABEL_DONE:

;;;
;  Done!  Go back to the default Graphics and cursor in case they were changed.
;  Also erase the plot ranges if they originally were not defined.
;
if FullCursor then device,/CURSOR_CROSSHAIR,SET_GRAPHICS=OldGraphics,Bypass=0

if UndefinedPlot then begin
   !X.CRange = 0
   !Y.CRange = 0
endif

return
end   ;   RDPLOT   by   Joel Parker   18 May 94
