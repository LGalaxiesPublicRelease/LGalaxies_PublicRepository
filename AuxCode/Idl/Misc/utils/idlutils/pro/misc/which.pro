;=================================================================

function which_find_routine, proname, _REF_EXTRA=_extra
; LOOKS FOR A MATCH BETWEEN ROUTINE INFORMATION AND
; AN IDL MODULE NAME...
; CLEVERLY, COMPILATION OF WHICH GUARANTEES THAT THERE WILL ALWAYS
; BE AT LEAST ONE PROCEDURE (WHICH) AND FUNCTION (WHICH_FIND_ROUTINE)...
compile_opt idl2, hidden
return, strmatch(routine_info(_EXTRA=_extra), proname, /FOLD_CASE)
end; which_find_routine

;=================================================================

pro which, proname
;+
; NAME:
;       WHICH
;
; PURPOSE: 
;       To search for any file in the IDL !path that contains the
;       user-supplied IDL routine (procedure or function) name.  Also
;       returns compilation status of each routine (in IDL lingo,
;       whether or not the routine is "resolved".)
;
; CALLING SEQUENCE:
;       WHICH, file
;
; INPUTS:
;       FILE - file name to search for.  The suffix .pro will be
;              appended if not included.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS: 
;       The IDL !path is searched for file names that are simply the
;       module (in IDL documentation, "module" and "routine" are used
;       interchangeably) name with a ".pro" suffix appended to them.
;       A module stored inside a file whose name is different than the
;       module name (followed by a ".pro") will not be found UNLESS
;       that module happens to be the currently-resolved module!
;       E.g., if the module "pro test_proc" lives in a file named
;       "dumb_name.pro", then it will not be found:
;
;       IDL> which, 'test_proc'
;       Module TEST_PROC Not Compiled.
;       % WHICH: test_proc.pro not found on IDL !path.
;
;       unless it happens to be resolved:
;
;       IDL> .run dumb_name
;       % Compiled module: TEST_PROC.
;       IDL> which, 'test_proc'
;       Currently-Compiled Module TEST_PROC in File:
;       /hvc/robishaw/dumb_name.pro
;
;       However, this is terrible programming style and sooner or
;       later, if you hide generically-named modules in
;       inappropriately-named files, bad things will (deservedly)
;       happen to you.
;
;       The routine further assumes that a file named "dumb_name.pro"
;       actually contains a module named "dumb_name"!  If it doesn't,
;       then you are a bad programmer and should seek professional
;       counseling.
; 
; PROCEDURES CALLED:
;       STRSPLIT(), WHICH_FIND_ROUTINE()
;
; EXAMPLES:
;       You haven't yet resolved (compiled) the routine (module)
;       DEFROI.  Let's look for it anyway:
;
;         IDL> which, 'defroi
;         Module DEFROI Not Compiled.
;
;         Other Files Containing Module DEFROI in IDL !path:
;         /usr/local/rsi/idl/lib/defroi.pro
;
;       For some reason you have two modules with the same name.
;       (This can occur in libraries of IDL routines such as the
;       Goddard IDL Astronomy User's Library; an updated version of a
;       routine is stored in a special directory while the old version
;       is stored in its original directory.) Let's see which version
;       of the module ADSTRING we are currently using:
;
;         IDL> which, 'adstring.pro'
;         Currently-Compiled Module ADSTRING in File:
;         /hvc/robishaw/idl/goddard/pro/v5.4+/adstring.pro
;
;         Other Files Containing Module ADSTRING in IDL !path:
;         /hvc/robishaw/idl/goddard/pro/astro/adstring.pro
;
; NOTES:
;       First, all currently-compiled procedures and functions are searched.
;       Then the remainder of the IDL !path is searched.
;
; MODIFICATION HISTORY:
;   30 May 2003  Written by Tim Robishaw, Berkeley
;   17 Feb 2004  Fixed oddity where user tries to call a function as
;                if it were a procedure, thus listing the module in both
;                the Compiled Functions and Compiled Procedures list.
;-

on_error, 2
resolve_routine, 'strsplit', /IS_FUN, /NO_RECOMPILE

if (N_params() lt 1) then begin
    message, 'syntax: which, proname (suffix .pro assumed)', /INFO
    return
endif

; IF .PRO SUFFIX INCLUDED, DROP IT...
proname = strtrim(proname,2)
if strmatch(proname,'*.pro', /FOLD_CASE) $
  then proname = strmid(proname,0,strlen(proname)-4)

; SEARCH THE CURRENTLY-COMPILED PROCEDURES AND FUNCTIONS FIRST...
pindx = where(which_find_routine(proname),presolved)
findx = where(which_find_routine(proname,/FUNCTIONS),fresolved)

; IF PROCEDURE OR FUNCTION WAS FOUND, IS IT UNRESOLVED...
punresolved = total(which_find_routine(proname,/UNRESOLVED))
funresolved = total(which_find_routine(proname,/UNRESOLVED,/FUNCTIONS))

if (presolved and not punresolved) OR $
   (fresolved and not funresolved) then begin

    ; THE PROCEDURE OR FUNCTION WAS FOUND...
    resolved_routine = (presolved AND not fresolved) ? $
      (routine_info(/SOURCE))[pindx].PATH : $
      (routine_info(/SOURCE,/FUNCTIONS))[findx].PATH

    print, 'Currently-Compiled Module '+strupcase(proname)+' in File:'
    print, resolved_routine, format='(A,%"\N")'

endif $
else print, strupcase(proname), format='("Module ",A," Not Compiled.",%"\N")'

; EXTRACT THE !PATH INTO A STRING ARRAY...
path = strsplit(!path, ':', /EXTRACT)

; GET RID OF "." IF USER INCLUDES THIS IN PATH...
path = path[where(path ne '.')]

; SEARCH CURRENT DIRECTORY, EVEN IF NOT IN IDL PATH...
cd, CURRENT=current
if (total(strmatch(path,current)) eq 0) then path = [current,path]

; ADD THE FILENAME TO EACH PATH DIRECTORY...
filenames = path + '/' + proname + '.pro'

; DOES ANY SUCH FILE EXIST IN THE CURRENT PATH...
file_exists = where(file_test(filenames), N_exists)

; IF THERE IS NO SUCH FILE THEN SPLIT...
if (N_exists eq 0) then begin
    if (N_elements(resolved_routine) eq 0) then $
        message, proname + '.pro not found on IDL !path.', /INFO
    return
endif

; PULL OUT ALL THE FILES THAT EXIST...
filenames = filenames[file_exists]

; TAKE RESOLVED ROUTINE OUT OF THE LIST...
if (N_elements(resolved_routine) gt 0) then begin

    ; GET THE INDICES OF THE UNRESOLVED ROUTINES...
    file_exists = where(strmatch(filenames,resolved_routine) eq 0, N_exists)

    ; WAS THE RESOLVED ROUTINE THE ONLY ONE...
    if (N_exists eq 0) then return

    filenames = filenames[file_exists]
endif

; PRINT THE REMAINING ROUTINES...
print, 'Other Files Containing Module '+strupcase(proname)+' in IDL !path:'
print, transpose(filenames)
print

end; which
