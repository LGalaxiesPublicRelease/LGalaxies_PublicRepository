function date_conv,date,type
;+
; NAME:
;     DATE_CONV
; PURPOSE:
;     Procedure to perform conversion of dates to one of three possible formats.
;
; EXPLANATION:
;     The following date formats are allowed
;
;       format 1: real*8 scalar encoded as:
;               year*1000 + day + hour/24. + min/24./60 + sec/24./60/60
;               where day is the day of year (1 to 366)
;       format 2: Vector encoded as:
;               date[0] = year (eg. 1987)
;               date[1] = day of year (1 to 366)
;               date[2] = hour
;               date[3] = minute
;               date[4] = second
;       format 3: string (ascii text) encoded as
;               DD-MON-YEAR HH:MM:SS.SS
;               (eg.  14-JUL-1987 15:25:44.23)
;            OR
;               YYYY-MM-DD HH:MM:SS.SS  (ISO standard)
;               (eg.  1987-07-14 15:25:44.23 or 1987-07-14T15:25:44.23)
;	            
;       format 4: three element vector giving spacecraft time words
;       from a Hubble Space Telescope (HST) telemetry packet.
;
; CALLING SEQUENCE
;       results = DATE_CONV( DATE, TYPE )
;
; INPUTS:
;       DATE - input date in one of the three possible formats.
;       TYPE - type of output format desired.  If not supplied then
;               format 3 (real*8 scalar) is used.
;                       valid values:
;                       'REAL'  - format 1
;                       'VECTOR' - format 2
;                       'STRING' - format 3
;			'FITS' - YYYY-MM-DDTHH:MM:SS.SS'
;               TYPE can be abbreviated to the single character strings 'R',
;               'V', 'S' and 'F'.
;               Nobody wants to convert TO spacecraft time (I hope!)
; OUTPUTS:
;       The converted date is returned as the function value.
;
; NOTES:
;      Prior to Oct 1998, the returned real*8 date (format 1) was given as
;      (year-1900)*1000 + day + hour/24. + min/24./60 + sec/24./60/60 
;      This output is ambiguous with respect to the year 2000.   Note that the
;      current version of DATE_CONV() may not be backwards compatible with 
;      versions prior to Oct 1998.
;
; HISTORY:
;      version 1  D. Lindler  July, 1987
;      adapted for IDL version 2  J. Isensee  May, 1990
;      Made year 2000 compliant; allow ISO format input  jls/acc Oct 1998
;      DJL/ACC Jan 1998, Modified to work with dates such as 6-JAN-1996 where
;		day of month has only one digit.
;      DJL, Nov. 2000, Added input/output format YYYY-MM-DDTHH:MM:SS.SS
;-
;-------------------------------------------------------------
;
; data declaration
;
days = [0,31,28,31,30,31,30,31,31,30,31,30,31]
months = ['   ','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT',$
        'NOV','DEC']
;
; set default type if not supplied
;
if n_params(0) lt 2 then type = 'REAL'
;
; Determine type of input supplied
;
s = size(date) & ndim = s[0] & datatype = s[ndim+1]
if ndim gt 0 then begin                 ;vector?
        if ndim gt 1 then goto,notvalid
        if (s[1] ne 5) and (s[1] ne 3) then goto,notvalid
        if (s[1] eq 5) then form = 2 else form = 4
   end else begin                       ;scalar input
        if datatype eq 0 then goto,notvalid
        if datatype eq 7 then form = 3 $        ;string
                         else form = 1  ;numeric scalar
end
;
;      -----------------------------------
;
;*** convert input to year,day,hour,minute,second
;
;      -----------------------------------
case form of

        1: begin                                        ;real scalar
                idate = long(date)
                year = long(idate/1000)
;
; if year is only 2 digits, assume 1900
;
                if year lt 100 then begin
                   print,'DATE_CONV: Warning: Year specified is only 2 digits, assuming 19xx'
                   year=1900+year
                   idate=1900000+idate
                   date=1900000.+date
                end
;
                day = idate - year*1000
                fdate = date-idate
                fdate = fdate*24.
                hour = fix(fdate)
                fdate = (fdate-hour)*60.0
                minute = fix(fdate)
                sec = float((fdate-minute)*60.0)
           end

        2: begin                                        ;vector
                year = fix(date[0])
                if year lt 100 then begin
                   print,'DATE_CONV: Warning: Year specified is only 2 digits, assuming 19xx'
                   year=1900+year
                end
;
; if year is only 2 digits, assume 1900
;
                if year lt 100 then begin
                   print,'DATE_CONV: Warning: Year specified is only 2 digits, assuming 19xx'
                   year=1900+year
                end
;
                day = fix(date[1])
                hour = fix(date[2])
                minute = fix(date[3])
                sec = float(date[4])
           end

        3: begin                                        ;string
                temp = date
;
; check for old type of date, DD-MMM-YYYY
;
                if strpos(temp,'-') le 2 then begin
                  day_of_month = fix(gettok(temp,'-'))
                  month_name = gettok(temp,'-')
                  year = fix(gettok(temp,' '))
                  hour = fix(gettok(temp,':'))
                  minute = fix(gettok(temp,':'))
                  sec = float(strtrim(strmid(temp,0,5)))
;
; determine month number from month name
;
                  month_name = strupcase(month_name)
                  for mon = 1,12 do begin
                        if month_name eq months[mon] then goto,found
                  end
                  print,'DATE_CONV -- invalid month name specified'
                  retall
;
; check for new type of date, ISO: YYYY-MM-DD
;
                end else if strpos(temp,'-') eq 4 then begin
                  year = fix(gettok(temp,'-'))
                  month_name = gettok(temp,'-')
                  mon=month_name
                  day_of_month=gettok(temp,' ')
                  if strlen(temp) eq 0 then begin
                        dtmp=gettok(day_of_month,'T')
                        temp=day_of_month
                        day_of_month=dtmp
                  end
                  day_of_month=fix(day_of_month)
                  hour = fix(gettok(temp,':'))
                  minute = fix(gettok(temp,':'))
                  sec = float(strtrim(strmid(temp,0,5)))
                end else goto, notvalid
              found:
;
; if year is only 2 digits, assume 1900
;
                if year lt 100 then begin
                   print,'DATE_CONV: Warning: Year specified is only 2 digits, assuming 19xx'
                   year=1900+year
                end
;
;
;            convert to day of year from month/day_of_month
;
;            correction for leap years
;
;               if (fix(year) mod 4) eq 0 then days(2) = 29     ;add one to february
                lpyr = ((year mod 4) eq 0) and ((year mod 100) ne 0) $
                        or ((year mod 400) eq 0)
                if lpyr eq 1 then days[2] = 29 ; if leap year, add day to Feb.
;
;
;            compute day of year
;
                  day = fix(total(days[0:mon-1])+day_of_month)
           end

        4 : begin                       ;spacecraft time
                SC = DOUBLE(date)
                SC = SC + (SC LT 0.0)*65536.    ;Get rid of neg. numbers 
;
;            Determine total number of secs since midnight, JAN. 1, 1979
;
                SECS = SC[2]/64 + SC[1]*1024 + SC[0]*1024*65536.
                SECS = SECS/8192.0D0            ;Convert from spacecraft units 
;
;            Determine number of years 
;
                MINS = SECS/60.
                HOURS = MINS/60.
                TOTDAYS = HOURS/24.
                YEARS = TOTDAYS/365.
                YEARS = FIX(YEARS)
;
;            Compute number of leap years past 
;
                LEAPYEARS = (YEARS+2)/4
;
;           Compute day of year 
;
                DAY = FIX(TOTDAYS-YEARS*365.-LEAPYEARS)
;
;           Correct for case of being right at end of leapyear
;
                IF DAY LT 0 THEN BEGIN
                  DAY = DAY+366
                  LEAPYEARS = LEAPYEARS-1
                  YEARS = YEARS-1
                END
;
;            COMPUTE HOUR OF DAY
;
                TOTDAYS = YEARS*365.+DAY+LEAPYEARS
                HOUR = FIX(HOURS - 24*TOTDAYS)
                TOTHOURS = TOTDAYS*24+HOUR
;
;            COMPUTE MINUTE
;
                MINUTE = FIX(MINS-TOTHOURS*60)
                TOTMIN = TOTHOURS*60+MINUTE
;
;            COMPUTE SEC
;
                SEC = SECS-TOTMIN*60
;
;            COMPUTE ACTUAL YEAR
;
                YEAR = YEARS+79
;
; if year is only 2 digits, assume 1900
;
                if year lt 100 then begin
                   print,'DATE_CONV: Warning: Year specified is only 2 digits, assuming 19xx'
                   year=1900+year
                end
;
;
;            START DAY AT ONE AND NOT ZERO
;
                DAY=DAY+1
           END
ENDCASE
;
;            correction for leap years
;
        if form ne 3 then begin         ;Was it already done?
           lpyr = ((year mod 4) eq 0) and ((year mod 100) ne 0) $
                or ((year mod 400) eq 0)
           if lpyr eq 1 then days[2] = 29 ; if leap year, add day to Feb.
        end
;
;            find month which day occurs
;
        day_of_month = day
        month_num = 1
        while day_of_month gt days[month_num] do begin
               day_of_month = day_of_month - days[month_num]
               month_num = month_num+1
        end
;
;            check for valid day
;
        if (day lt 1) or (day gt total(days)) then begin
           print,'DATE_CONV -- There are only',total(days),' in year ',year
           retall
        end
;           ---------------------------------------
;
;   *****       Now convert to output format
;
;           ---------------------------------------
;
; is type a string
;
s = size(type)
if (s[0] ne 0) or (s[1] ne 7) then begin
        print,'DATE_CONV- Output type specification must be a string'
        retall
end
;
case strmid(strupcase(type),0,1) of

        'V' : begin                             ;vector output
                out = fltarr(5)
                out[0] = year
                out[1] = day
                out[2] = hour
                out[3] = minute
                out[4] = sec
             end
 
        'R' : begin                             ;floating point scalar
;               if year gt 1900 then year = year-1900
                out = sec/24.0d0/60./60. + minute/24.0d0/60. + hour/24.0d0 $
                        +  day + year*1000d0
              end

        'S' : begin                             ;string output 

                month_name = months[month_num]
;
;            encode into ascii_date
;
                out = string(day_of_month,'(i2)') +'-'+ month_name +'-' + $
                        string(year,'(i4)') + ' '+ $
                        string(hour,'(i2)') +':'+ $
                        strmid(string(minute+100,'(i3)'),1,2) + ':'+ $
                        strmid(string(sec+100,'(f6.2)'),1,5)
           end
	'F' : begin
		out = string(year,'(i4)')+'-'+string(month_num,'(I2)')+'-'+ $
			string(day_of_month,'(i2)')+'T'+string(hour,'(i2)') +':'+ $
                        string(minute,'(i2.2)') + ':'+ $
                        strmid(string(sec+100,'(f6.2)'),1,5)
	      end
        else: begin                     ;invalid type specified
                print,'DATE_CONV-- Invalid output type specified'
                print,' It must be ''REAL'', ''STRING'', or ''VECTOR'''
                retall
              end
endcase
return,out
;
; invalid input date error section
;
notvalid:
print,'DATE_CONV -- invalid input date specified'
retall
end
