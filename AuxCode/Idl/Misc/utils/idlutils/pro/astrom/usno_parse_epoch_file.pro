pro usno_parse_epoch_file
  
  path = '/u/dss/data/usnob/'
  file = path+'allEpochs.dat'
  if not file_test(file) then message, 'cannot find file'
  spawn, 'wc -l '+file, res
  len = long((strsplit(res, /ex))[0])
  data = dblarr(12, len)
  openr, rlun, file, /get_lun
  readf, rlun, data
  free_lun, rlun
  plate = data[0, *]
  epoch = data[1:10, *]
  sbox = lindgen(10, len) mod 10
  pbox = lindgen(10, len)  /  10
  field = 1+pbox+sbox*1000L
  w = where(epoch ge 0, ngood)
  
  str = {field: 0L, epoch:0.d}
  ep = replicate(str, ngood)
  ep.field = field[w]
  ep.epoch = epoch[w]

  mwrfits, ep, 'USNO-B-epochs.fit', /create

  return
end
