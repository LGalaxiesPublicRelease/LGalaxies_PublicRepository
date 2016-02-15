;+
; NAME:
;   ppmtompeg
;
; PURPOSE:
;   Wrapper for ppmtompeg, the open-source UNIX mpeg writer
;
; CALLING SEQUENCE:
;   ppmtompeg, bytcube, [ output, tmpdir= ]
; 
; INPUTS:
;   bytcube  - Byte array [X,Y,NFRAME]
;   output   - Output file name; default 'idl.mpeg'
;
; KEYWORDS:
;   tmpdir   - temp directory to do dirty work in; default to '/tmp/mpeg1234/'
;
; OUTPUTS:
;
; COMMENTS:
;   A single mpeg file is written.
;   Temporary files are removed after the MPEG file is made.
;
; EXAMPLES:
;   ppmtompeg, bytscl(imagecube, min=800, max=1600), 'movie.mpeg'
;
; PROCEDURES CALLED:
;   rmfile
;
; INTERNAL SUPPORT PROCEDURES:
;   ppmtompeg_parameters
;
; REVISION HISTORY:
;   2002-Mar-29  Douglas Finkbeiner, Princeton University
;                  dfink@astro.princeton.edu
;
;----------------------------------------------------------------------

; Write parameter file
pro ppmtompeg_parameters, parname, input_dir=input_dir, $
             input_name=input_name, $
             yuv_size=yuv_size, frame_rate=frame_rate, output=output

  if NOT keyword_set(input_dir) then input_dir = '.'
  if NOT keyword_set(input_name) then message, 'input_name required!'
  if NOT keyword_set(yuv_size) then  yuv_size = '320x240'
  if NOT keyword_set(frame_rate) then frame_rate = 30
  if NOT keyword_set(output) then output = 'idl.mpeg'
  if NOT keyword_set(parname) then message, 'parname required!'

  openw, wlun, parname, /get_lun

  printf, wlun, '# Written by IDL ', systime()
  printf, wlun, 'PATTERN         IBBPBBPBBPBBPBB'
  printf, wlun, 'OUTPUT          '+output
  printf, wlun, 'GOP_SIZE        15'
  printf, wlun, 'SLICES_PER_FRAME        1'
  printf, wlun, 'BASE_FILE_FORMAT        PPM'
  printf, wlun, 'YUV_SIZE        '+yuv_size
  printf, wlun, 'INPUT_CONVERT   *'
  printf, wlun, 'INPUT_DIR       '+input_dir
  printf, wlun, 'INPUT'
  printf, wlun, input_name
  printf, wlun, 'END_INPUT'
  printf, wlun, '# motion vector search parameters'
  printf, wlun, '# MAD or MSE -- must be upper case'
  printf, wlun, 'ERROR           MAD'
  printf, wlun, '# FULL or HALF -- must be upper case'
  printf, wlun, 'PIXEL           FULL'
  printf, wlun, '# means +/- this many pixels'
  printf, wlun, 'RANGE           8'
  printf, wlun, ''
  printf, wlun, 'PSEARCH_ALG     EXHAUSTIVE'
  printf, wlun, 'BSEARCH_ALG     SIMPLE'
  printf, wlun, ''
  printf, wlun, '#FRAME_RATE     '+frame_rate
  printf, wlun, 'FORCE_ENCODE_LAST_FRAME'
  printf, wlun, ''
  printf, wlun, 'IQSCALE         1'
  printf, wlun, 'PQSCALE         6'
  printf, wlun, 'BQSCALE         6'
  printf, wlun, 'REFERENCE_FRAME ORIGINAL'

  free_lun, wlun
  return
end


; Make the mpeg
pro ppmtompeg, bytcube, output, tmpdir=tmpdir

; defaults
  if NOT keyword_set(output) then output = 'idl.mpeg'

; temporary directory
  if NOT keyword_set(tmpdir) then tmpdir = '/tmp/mpeg1234/'
  spawn, 'mkdir -p '+tmpdir

; size string
  sz = size(bytcube, /dim)
  yuv_size = strcompress(string(sz[0], sz[1], format='(I,"x",I)'), /rem)

; write ppm files
  nframe = sz[2]
  fname = 'f-'+string(lindgen(nframe), format='(I6.6)')+'.ppm'
  for i=0L, nframe-1 do begin 
     write_ppm, tmpdir+fname[i], bytcube[*, *, i]
  endfor

; generate param file for ppmtompeg
  input_dir = tmpdir
  input_name = 'f-*.ppm    [000000-'+string(nframe-1, format='(I6.6)')+']'
  frame_rate = '30'
  parname = tmpdir+'ppmtompeg.param'
  ppmtompeg_parameters, parname, input_dir=input_dir, $
    input_name=input_name, $
    yuv_size=yuv_size, frame_rate=frame_rate, output=output

; generate mpeg
  spawn, 'ppmtompeg -quiet '+parname

; Remove temporary files
  if (strmid(tmpdir, 0, 5) EQ '/tmp/') then begin
    spawn, '\rm -r '+tmpdir
  endif else begin
    rmfile, parname
    rmfile, tmpdir+fname
  endelse

  return
end
