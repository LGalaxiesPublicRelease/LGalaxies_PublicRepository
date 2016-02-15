;+
; NAME:
;   group_indx
; PURPOSE:
;   given group indices, yield multiplicity, plus an index linked list
; CALLING SEQUENCE:
;   group_indx, ingroup, multgroup=multgroup, firstgroup=firstgroup, $
;      nextgroup=nextgroup
; INPUTS:
;   ingroup - [N] group # of each element
; OUTPUTS:
;   multgroup - [Ngroup] multiplicity of each group
;   firstgroup - [Ngroup] first member of each group
;   nextgroup - [N] for each member, the next member of its group
;               (-1) if no more
; COMMENTS:
; BUGS:
; REVISION HISTORY:
;   2003-03-05 written by Michael Blanton (NYU)
;-
pro group_indx, ingroup, multgroup=multgroup, firstgroup=firstgroup, $
                nextgroup=nextgroup

npoints=n_elements(ingroup)
ngroups=max(ingroup)+1L

multgroup=lonarr(ngroups)
nextgroup=lonarr(npoints)-1L
firstgroup=lonarr(ngroups)-1L
for i = 0l, npoints-1l do multgroup[ingroup[i]]=multgroup[ingroup[i]]+1l
for i = npoints-1l, 0l, -1l do begin 
    nextgroup[i]=firstgroup[ingroup[i]]
    firstgroup[ingroup[i]]=i
end

end
