;+
; NAME: 
;  IDL_VALIDNAME()
;
; PURPOSE: 
;   Modify a string if necessary, so that it can used as a IDL variable name.
;   
; EXPLANATION:
;   Duplicates the intrinsic V6.0 function IDL_VALIDNAME() with the 
;       CONVERT_ALL keyword, for  pre-V6.0 compatibility
;   IDL_VALIDNAME performs the following:
;   (1) All non-alphanumeric characters (except '_' and '$') are converted
;       to underscores
;   (2) If the first character is not alphabetic or a '_' or a '!' then an 
;       underscore is preprended to the string
;   (3) If the string is an IDL reserved word (e.g. 'endif') then an 
;       underscore is preprended to the string.
;
; CALLING SEQUENCE: 
;    result = IDL_VALIDNAME( str, /Convert_all ) 
;
; INPUT: 
;      str - Scalar string to be converted to a valid variable name. 
;
; OUTPUT:
;      result - the input string modified, if necessary, according to the
;           3 rules above, so that it can be used as a valid IDL variable.
; OPTIONAL INPUT KEYWORDS:
;    /CONVERT_ALL -- Does nothing, but ensures compatibility with V6.0
;          intrinsic call IDL_VALIDNAME(/CONVERT_ALL)
; PROCEDURES USED:
;     None.
;
; EXAMPLES:
;      (1) IDL> print,idl_validname('switch',/convert_all)
;          _switch           ;IDL reserved keyword
;
;      (2) IDL> print,idl_validname('2mass',/convert_all)
;          _2mass            ;1st character must be alphabetic
;
;      (3) IDL> print,idl_validname('date$',/convert_all)
;          date$             ;no change, input string already valid
; MODIFICATION HISTORY: 
;         Written by W. Landsman  SSAI  October 2003
;
;-
    function idl_validname, st, convert_all = convert_all

    reserved_names = ['AND','BEGIN','BREAK','CASE','COMMON','COMPILE_OPT','DO',$
     'ELSE','END','ENDCASE', 'ENDELSE','ENDFOR','ENDIF','ENDREP','ENDSWITCH', $
     'ENDWHILE','EQ','FOR','FORWARD_FUNCTION','FUNCTION','GE','GOTO','GT','IF',$
     'INHERITS','LE','LT','MOD','NE','NOT','OF','ON_IOERROR', $
     'OR','PRO','REPEAT','SWITCH','THEN','UNTIL','WHILE','XOR' ]
  
    reserved = where(strupcase(st) EQ reserved_names, Nreserv)
    if Nreserv GT 0 then st = '_' + st else begin
     
    ; 
    ; Replace invalid characters with underscores.    Only valid special
    ; characters are  _ (95b) and $(36b).   We assume ASCII ordering.
    ;
    bst = byte(st)
  
    if N_elements(bst) GT 1 then begin   ;More than 1 character?
      b = bst[1:*]
      btest =   ((b GE 97b  and b LE 122b)  or  $    ;small case letter
                 (b GE 65b  and b LE 90b)  or  $    ;upper case letter
                 (b GE 48b  and b LE 57b)  or  $    ;numeric
                 (b eq 95b) or (b eq 36b)  ) 
   
      btest = [1b, btest]
      g  = where(btest EQ 0, Ng)
      if Ng GT 0 then bst[g] =  95b           ;Underscore
    endif

    ; First character must be alphabetic or a underscore or a "!" 
 
    c = bst[0] 
    st = string(bst)

; If first character is numeric or a dollar sign then preprend a '_'

    btest =    ((c GE 97b  and c LE 122b)  or  $    ;small case letter
                 (c GE 65b  and c LE 90b)  or  $    ;upper case letter
                  (c eq 95b) or (c eq 33b)  ) 

    if btest EQ 0b then st = '_' + st 
    endelse
 
    return,st
    end
