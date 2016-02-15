;+
; NAME:
;   hogg_unquoted_regex
; PURPOSE:
;   return the regex which matches the first occurence of the given
;     regex not inside quotemarks
; INPUT:
;   regex      - naked regular expression
; OPTIONAL INPUT:
;   quotemark  - thing to use as the quotation mark, default to '"'
; REVISION HISTORY:
;   2002-10-11  written - Hogg
;-
function hogg_unquoted_regex, regex,quotemark=quotemark
if NOT keyword_set(quotemark) then quotemark='"'
return, '('+regex+')((.*(\".*\")+.*$)|([^\"]*$))'
end
