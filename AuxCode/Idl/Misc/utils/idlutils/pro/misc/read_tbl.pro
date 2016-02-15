;+
; NAME:
;   read_tbl
; PURPOSE:
;   Read .tbl/.hdr format files output by postgres (like those of 2MASS)
; CALLING SEQUENCE:
;   table= read_tbl(filebase [, chunksize=])
; INPUTS:
;   filebase - filebase for filebase.hdr and filebase.tbl files
; OPTIONAL INPUTS:
;   chunksize - size of chunks to read in at a time before converting [1000]
; COMMENTS:
;   Returns a structure in "table"
;   Not very extensively tested
; BUGS:
;   I don't know if I've included ALL possible types in type2var
; REVISION HISTORY:
;   2003-01-29 written by Michael Blanton (NYU)
;-
function type2var, type, val

if(strtrim(type,2) eq 'double') then $
  return,double(val) $
else if(strtrim(type,2) eq 'int') then $
  return,long(val) $
else if(strtrim(type,2) eq 'float') then $
  return,float(val) $
else if(strtrim(type,2) eq 'char') then $
  return,string(val) $
else $
  message, 'No registered type '+type

end
;
function read_tbl,filebase,chunksize=chunksize

if(n_params() ne 1) then begin
    doc_library,'read_tbl'
    return,-1
endif

; defaults
if(NOT keyword_set(chunksize)) then chunksize=1000

; build names
hdrfile=filebase+'.hdr'
tblfile=filebase+'.tbl'

; read in and interpret header
line=''
openr,unit,hdrfile,/get_lun
readf,unit,line
names=strtrim(strsplit(line,'|',/extract),2)
ncols=n_elements(names)
limits=strsplit(line,'|')
lengths=lonarr(ncols)
lengths[0:ncols-2]=limits[1:ncols-1]-limits[0:ncols-2]+1
lengths[ncols-1]=strlen(line)-limits[ncols-1]+1
readf,unit,line
types=strtrim(strsplit(line,'|',/extract),2)
free_lun,unit

; now build structure
for i=0L, n_elements(names)-1L do begin
    dummy=type2var(types[i],'')
    if(n_tags(instr1) eq 0) then $
      instr1=create_struct(names[i],dummy) $
    else $
      instr1=create_struct(instr1,names[i],dummy)
endfor

; create table 
nelem=numlines(tblfile)
instr=replicate(instr1,nelem)

; now read in table
chunk=strarr(ncols,chunksize)
openr,unit,tblfile,/get_lun
nchunks=long(ceil(double(nelem)/double(chunksize)))
for i=0L, nchunks-1L do begin
    startchunk=i*chunksize
    endchunk=(((i+1L)*chunksize)<nelem)-1L
    ninchunk=endchunk-startchunk+1L
    splog,'startchunk= '+string(startchunk)+'; endchunk= '+string(endchunk)
    for j=0L, ninchunk-1L do begin
        readf,unit,line
        chunk[*,j]=strmid(line,limits,lengths)
    endfor
    for j=0L,ncols-1L do begin
        instr[startchunk:endchunk].(j)= $
          transpose(type2var(types[j],chunk[j,0L:ninchunk-1L]))
    endfor
endfor
free_lun,unit

return,instr

end
