;+
; NAME:
;   mjd2datelist
;
; PURPOSE:
;   Construct a list of MJDs and date strings spanning a range of MJDs
;   (useful for plot limits).
;
; CALLING SEQUENCE:
;   mjd2datelist, mjstart, [ mjend, step=, mjdlist=, datelist= ]
;
; INPUTS:
;   mjstart    - Starting modified Julian date to span.
;
; OPTIONAL INPUTS:
;   mjend      - Ending modified Julian date to span; if not set, then
;                only the date string for MJSTART is returned.
;   step       - Step in either 'year', '6month', 'month', or 'day';
;                default to 'year', or 'day' if MJEND not set.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mjlist     - List of modified Julian dates (MJDs)
;   datelist   - List of dates in the form DD-MMM-YYYY, where MMM is
;                the first three letters of the month name.
;
; COMMENTS:
;   This routine returns a list of MJDs and date strings spaced by
;   the amount specified by STEP that span the range [MJSTART,MJEND].
;   If using STEP='year', the output list will be on the first date
;   of each year and [MJSTART,MJEND] will fall internal to that list.
;   If using STEP='month', the output list will be on each 01-Jan
;   and 01-Jul.  If using STEP='month', the output list will be on
;   the first of each month.
;
; EXAMPLES:
;   Construct a list of all the Jan 1st dates that span the dates
;   of the SDSS spectroscopic survey:
;     IDL> mjd2datelist, 51433, 52356, mjdlist=mjdlist, datelist=datelist
;     IDL> print, mjdlist
;          51179.500   51544.500   51910.500   52275.500   52640.500
;     IDL> print, datelist
;          01-Jan-1999 01-Jan-2000 01-Jan-2001 01-Jan-2002 01-Jan-2003
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   26-Mar-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro mjd2datelist, mjstart, mjend, step=step, $
 mjdlist=mjdlist, datelist=datelist

   if (n_params() LT 1) then begin
      doc_library, 'mjd2datelist'
      return
   endif
   if (NOT keyword_set(step)) then begin
      if (keyword_set(mjend)) then step = 'year' $
       else step = 'day'
   endif

   monthname = ['','Jan','Feb','Mar','Apr','May','Jun', $
    'Jul','Aug','Sep','Oct','Nov','Dec']
   numdays = [0,31,28,31,30,31,30,31,31,30,31,30,31]

   offset = 2400000.5D

   jd1 = offset + mjstart
   caldat, jd1, month1, day1, year1

   ; Force back to the first day of the month, 6-month, or year
   case strupcase(step) of
   'YEAR' : begin
      month1 = 1
      day1 = 1
      end
   '6MONTH' : begin
      if (month1 GE 6) then month1 = 7 $
       else month1 = 1
      day1 = 1
      end
   'MONTH' : begin
      day1 = 1
      end
   'DAY' : begin
      end
   endcase
   mjdlist = julday(month1, day1, year1) - offset
   jd1 = mjdlist + offset

   datelist = string(day1, monthname[month1], year1, $
    format='(i2.2,"-",a3,"-",i4.4)')
   if (NOT keyword_set(mjend)) then return

   while (max(mjdlist) LT mjend) do begin
      case strupcase(step) of
      'YEAR' : begin
         year1 = year1 + 1
         end
      '6MONTH' : begin
         month1 = month1 + 6
         if (month1 EQ 13) then begin
            month1 = 1
            year1 = year1 + 1
         endif
         end
      'MONTH' : begin
         month1 = month1 + 1
         if (month1 EQ 13) then begin
            month1 = 1
            year1 = year1 + 1
         endif
         end
      'DAY' : begin
         ; Check for leap year (if February)
         if (month1 EQ 2) then begin
            qleap = ((year1 MOD 4) EQ 0) $
             AND ( ((year1 MOD 100) NE 0) OR ((year1 MOD 400) EQ 0) )
            numdays[month1] = 28 + qleap
         endif
         day1 = day1 + 1
         if (day1 GT numdays[month1]) then begin
            day1 = 1
            month1 = month1 + 1
            if (month1 EQ 13) then begin
               month1 = 1
               year1 = year1 + 1
            endif
         endif
         end
      else: message, 'Unknown value for STEP'
      endcase

      jd1 = julday(month1, day1, year1)

      mjdlist = [mjdlist, jd1 - offset]
      datelist = [datelist, $
       string(day1, monthname[month1], year1, $
       format='(i2.2,"-",a3,"-",i4.4)') ]
   endwhile

end
;------------------------------------------------------------------------------
