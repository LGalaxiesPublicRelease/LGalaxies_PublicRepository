;+
; NAME:
;    ERROR_MESSAGE
;
; PURPOSE:
;
;    The purpose of this function  is to have a device-independent
;    error messaging function. The error message is reported
;    to the user by using DIALOG_MESSAGE if widgets are
;    supported and MESSAGE otherwise.
;
;    In general, the ERROR_MESSAGE function is not called directly.
;    Rather, it is used in a CATCH error handler. Errors are thrown
;    to ERROR_MESSAGE with the MESSAGE command. A typical CATCH error
;    handler is shown below.
;
;       Catch, theError
;       IF theError NE 0 THEN BEGIN
;          Catch, /Cancel
;          ok = Error_Message(/Traceback, /Error)
;          RETURN
;       ENDIF
;
;    Error messages would get into the ERROR_MESSAGE function by
;    throwing an error with the MESSAGE command, like this:
;
;       IF test NE 1 THEN Message, 'The test failed.'
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;    Utility.
;
; CALLING SEQUENCE:
;
;    ok = Error_Message(the_Error_Message)
;
; INPUTS:
;
;    the_Error_Message: This is a string argument containing the error
;       message you want reported. If undefined, this variable is set
;       to the string in the !Error_State.Msg system variable.
;
; KEYWORDS:
;
;    ERROR: Set this keyword to cause Dialog_Message to use the ERROR
;       reporting dialog. Note that a bug in IDL causes the ERROR dialog
;       to be used whether this keyword is set to 0 or 1!
;
;    INFORMATIONAL: Set this keyword to cause Dialog_Message to use the
;       INFORMATION dialog instead of the WARNING dialog. Note that a bug
;       in IDL causes the ERROR dialog to be used if this keyword is set to 0!
;
;    TITLE: Set this keyword to the title of the DIALOG_MESSAGE window. By
;       default the keyword is set to 'System Error' unless !ERROR_STATE.NAME
;       equals "IDL_M_USER_ERR", in which case it is set to "Trapped Error'.
;
;    TRACEBACK: Setting this keyword results in an error traceback
;       being printed to standard output with the PRINT command. Set to
;       1 (ON) by default. Use TRACEBACK=0 to turn this functionality off.
;
; OUTPUTS:
;
;    Currently the only output from the function is the string "OK".
;
; RESTRICTIONS:
;
;    The WARNING Dialog_Message dialog is used by default.
;
; EXAMPLE:
;
;    To handle an undefined variable error:
;
;    IF N_Elements(variable) EQ 0 THEN $
;       ok = Error_Message('Variable is undefined', /Traceback)
;
; MODIFICATION HISTORY:
;
;    Written by: David W. Fanning, 27 April 1999.
;    Added the calling routine's name in the message and NoName keyword. 31 Jan 2000. DWF.
;    Added _Extra keyword. 10 February 2000. DWF.
;    Forgot to add _Extra everywhere. Fixed for MAIN errors. 8 AUG 2000. DWF.
;    Adding call routine's name to Traceback Report. 8 AUG 2000. DWF.
;    Added ERROR, INFORMATIONAL, and TITLE keywords. 19 SEP 2002. DWF.
;    Removed the requirement that you use the NONAME keyword with the MESSAGE
;      command when generating user-trapped errors. 19 SEP 2002. DWF.
;    Added distinctions between trapped errors (errors generated with the
;      MESSAGE command) and IDL system errors. Note that if you call ERROR_MESSAGE
;      directly, then the state of the !ERROR_STATE.NAME variable is set
;      to the *last* error generated. It is better to access ERROR_MESSAGE
;      indirectly in a Catch error handler from the MESSAGE command. 19 SEP 2002. DWF.
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 1999-2002 Fanning Software Consulting
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


FUNCTION ERROR_MESSAGE, theMessage, Error=error, Informational=information, $
   Traceback=traceback, NoName=noname, Title=title, _Extra=extra

On_Error, 2

   ; Check for presence and type of message.

IF N_Elements(theMessage) EQ 0 THEN theMessage = !Error_State.Msg
s = Size(theMessage)
messageType = s[s[0]+1]
IF messageType NE 7 THEN BEGIN
   Message, "The message parameter must be a string.", _Extra=extra
ENDIF

   ; Get the call stack and the calling routine's name.

Help, Calls=callStack
IF Float(!Version.Release) GE 5.2 THEN $
   callingRoutine = (StrSplit(StrCompress(callStack[1])," ", /Extract))[0] ELSE $
   callingRoutine = (Str_Sep(StrCompress(callStack[1])," "))[0]

   ; Are widgets supported?

widgetsSupported = ((!D.Flags AND 65536L) NE 0)
IF widgetsSupported THEN BEGIN

      ; If this is an error produced with the MESSAGE command, it is a trapped
      ; error and will have the name "IDL_M_USER_ERR".

   IF !ERROR_STATE.NAME EQ "IDL_M_USER_ERR" THEN BEGIN

      IF N_Elements(title) EQ 0 THEN title = 'Trapped Error'

         ; If the message has the name of the calling routine in it,
         ; it should be stripped out. Can you find a colon in the string?

      colon = StrPos(theMessage, ":")
      IF colon NE -1 THEN BEGIN

            ; Extract the text up to the colon. Is this the same as
            ; the callingRoutine? If so, strip it.

         IF StrMid(theMessage, 0, colon) EQ callingRoutine THEN $
            theMessage = StrMid(theMessage, colon+1)

      ENDIF

         ; Add the calling routine's name, unless NONAME is set.

      IF Keyword_Set(noname) THEN BEGIN
         answer = Dialog_Message(theMessage, Title=title, _Extra=extra, $
            Error=error, Information=information)
      ENDIF ELSE BEGIN
         answer = Dialog_Message(StrUpCase(callingRoutine) + ": " + $
            theMessage, Title=title, _Extra=extra, $
            Error=error, Information=information)
      ENDELSE

   ENDIF ELSE BEGIN

         ; Otherwise, this is an IDL system error.

      IF N_Elements(title) EQ 0 THEN title = 'System Error'

      IF StrUpCase(callingRoutine) EQ "$MAIN$" THEN $
         answer = Dialog_Message(theMessage, _Extra=extra, Title=title, $
            Error=error, Information=information) ELSE $
      IF Keyword_Set(noname) THEN BEGIN
         answer = Dialog_Message(theMessage, _Extra=extra, Title=title, $
            Error=error, Information=information)
      ENDIF ELSE BEGIN
         answer = Dialog_Message(StrUpCase(callingRoutine) + "--> " + $
            theMessage, _Extra=extra, Title=title, $
            Error=error, Information=information)
      ENDELSE
   ENDELSE
ENDIF ELSE BEGIN
      Message, theMessage, /Continue, /NoPrint, /NoName, /NoPrefix, _Extra=extra
      Print, '%' + callingRoutine + ': ' + theMessage
      answer = 'OK'
ENDELSE

   ; Provide traceback information if requested.

IF Keyword_Set(traceback) THEN BEGIN
   Help, /Last_Message, Output=traceback
   Print,''
   Print, 'Traceback Report from ' + StrUpCase(callingRoutine) + ':'
   Print, ''
   numlines = N_Elements(traceback)
   IF numlines EQ 1 THEN BEGIN
      Print, "     " + traceback[0]
   ENDIF ELSE BEGIN
      Print, "     " + traceback[0]
      Print, ""
      FOR j=1,numlines-1 DO Print, "     " + traceback[j]
   ENDELSE
ENDIF

RETURN, answer
END ; ----------------------------------------------------------------------------

