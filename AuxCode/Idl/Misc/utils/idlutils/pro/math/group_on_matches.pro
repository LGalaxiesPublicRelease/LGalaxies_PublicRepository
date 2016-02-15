;+
; NAME:
;   group_on_matches
; PURPOSE:
;   given a list of objects and matches between them, find the groups
; CALLING SEQUENCE:
;   group_on_matches, matches [, first=, next=, mult=, in=]
; INPUTS:
;   matches - [N] array of pointers to list of which objects are matched
; OUTPUTS:
;   first - [Ngroup] first member of each group
;   mult - [Ngroup] multiplicity of each group
;   in - [N] what group an object is in 
;   next - [N] next member of group an object is in (-1 if last)
; REVISION HISTORY:
;   2005-Oct-3  Written by Mike Blanton, NYU
;----------------------------------------------------------------------
pro group_on_matches, inmatches, first=first, next=next, mult=mult, in=in

first=lonarr(n_elements(inmatches))-1L
next=lonarr(n_elements(inmatches))-1L
mult=lonarr(n_elements(inmatches))
in=lonarr(n_elements(inmatches))-1L
mapgroup=lonarr(n_elements(inmatches))-1L

;; fix up matches to make sense
matches=ptrarr(n_elements(inmatches))
for i=0L, n_elements(matches)-1L do begin
    if(NOT keyword_set(inmatches[i])) then $
      matches[i]=ptr_new([i]) $
    else $
      matches[i]=ptr_new(*(inmatches[i]))
endfor 
for i=0L, n_elements(matches)-1L do begin
    imatch=*(matches[i])
    
    ;; make sure it matches self
    imatch=[imatch, i]
    
    ;; make sure matches are reciprocal
    for j=0L, n_elements(imatch)-1L do begin
        imatch2=*(matches[imatch[j]])
        imatch2=[imatch2, i]
        ptr_free, matches[imatch[j]]
        matches[imatch[j]]=ptr_new(imatch2)
    endfor

    ptr_free, matches[i]
    matches[i]=ptr_new(imatch)
endfor    
for i=0L, n_elements(matches)-1L do begin
    imatch=*(matches[i])
    isort=sort(imatch)
    iuniq=uniq(imatch[isort])
    imatch=imatch[isort[iuniq]]
    ptr_free, matches[i]
    matches[i]=ptr_new(imatch)
endfor

igroup=0L
for i=0L, n_elements(matches)-1L do begin
    imatch=*(matches[i])
    
    ;; find the earliest group anybody is attached to
    iearly=where(in[imatch] ge 0, nearly)
    minearly=igroup
    for j=0L, nearly-1L do begin
        checkearly=in[imatch[iearly[j]]]
        while(mapgroup[checkearly] ne checkearly) do $
          checkearly=mapgroup[checkearly]
        if(checkearly lt minearly) then $
          minearly=checkearly
    endfor

    if(minearly eq igroup) then begin
        ;; all matches are brand new, so we've found a new group
        mapgroup[igroup]=igroup
        in[imatch]=igroup
        igroup=igroup+1L
    endif else begin
        ;; there is an earlier group, so we want to assign all the
        ;; current objects to that group, and make sure that mapgroups
        ;; for all the other groups map to the earliest one too
        iearly=where(in[imatch] ge 0, nearly)
        for j=0L, nearly-1L do begin
            checkearly=in[imatch[iearly[j]]]
            while(mapgroup[checkearly] ne checkearly) do begin
                tmpearly=mapgroup[checkearly]
                mapgroup[checkearly]=minearly
                checkearly=tmpearly
            endwhile
            mapgroup[checkearly]=minearly
        endfor
        in[imatch]=minearly
    endelse
endfor

ngroups=0L
for i=0L, n_elements(mapgroup)-1L do begin
    if(mapgroup[i] ne -1) then begin
        if(mapgroup[i] eq i) then begin
            mapgroup[i]=ngroups
            ngroups=ngroups+1L
        endif else begin
            mapgroup[i]=mapgroup[mapgroup[i]]
        endelse
    endif
endfor

in=mapgroup[in]

npoints=n_elements(matches)
mult=lonarr(ngroups)
next=lonarr(npoints)-1L
first=lonarr(npoints)-1L
for i = 0l, npoints-1l do mult[in[i]]=mult[in[i]]+1l
for i = npoints-1l, 0l, -1l do begin 
    next[i]=first[in[i]]
    first[in[i]]=i
end

for i=0L, npoints-1L do $
  ptr_free, matches[i]

end
