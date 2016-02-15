; rewrite USNO-B data as FITS files for easier access (for Ivezic)

pro usnob2fits_1, path, subdir, fname, outpath, hash

  rec_len = 80L
  
  dirstr = string(subdir, format='(I3.3)')
  catfile = path+'/'+dirstr+'/'+fname+'.cat'
  outdir  = outpath+'/'+dirstr
  if not file_test(outdir, /dir) then spawn, 'mkdir -p '+outdir
  outfile = outdir+'/'+fname+'.fit'

  openr, readlun, catfile, /get_lun, /swap_if_big_endian
  nbyte = (fstat(readlun)).size
  nstars = nbyte/rec_len
  data = ulonarr(rec_len/4, nstars)
  readu, readlun, data
  free_lun, readlun

  usnostruct = usnob10_extract(data)
  usnostruct.fldepoch = hash[usnostruct.fldid]

; -------- write temporary file
  tempname = '/tmp/'+fname+'.fit'
  mwrfits, usnostruct, tempname, /create
  spawn, 'gzip -vf '+tempname+' ; \mv '+tempname+'.gz '+outdir

  return
end



pro usnob2fits

; -------- set paths
  usno_dir = getenv('USNO_DIR')
  print, 'Using input directory ', usno_dir
  dataroot = getenv('PHOTO_DATA')+'/'
  if dataroot eq '/' then message, 'you need to set PHOTO_DATA'
  path = usno_dir
  outpath = dataroot+'usnob/fits/'

; -------- set up epoch hash
  epochfile = path+'/USNO-B-epochs.fit'
  ep = mrdfits(epochfile, 1)
  hash = fltarr(10000)
  hash[ep.field] = ep.epoch

  for subdir = 0, 179 do begin 
     for i=0, 9 do begin 
        fname = 'b'+string(subdir*10+i, format='(I4.4)')
        print, fname
        usnob2fits_1, path, subdir, fname, outpath, hash
     endfor
  endfor 
  return
end
