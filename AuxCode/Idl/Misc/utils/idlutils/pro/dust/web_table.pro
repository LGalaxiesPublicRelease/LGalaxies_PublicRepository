; Make r_V tables for SFD dustmap web site (for Peregrine)
; 2004-Mar-24  D. Finkbeiner

pro web_table

  sdss_rv_table, source='Galaxy', zsource=0.0, fname='Galaxy-z=0.0.txt'
  sdss_rv_table, source='Galaxy', zsource=0.1, fname='Galaxy-z=0.1.txt'
  sdss_rv_table, source='Galaxy', zsource=0.2, fname='Galaxy-z=0.2.txt'
  sdss_rv_table, source='Galaxy', zsource=0.3, fname='Galaxy-z=0.3.txt'
  sdss_rv_table, source='Galaxy', zsource=0.4, fname='Galaxy-z=0.4.txt'
  sdss_rv_table, source='Galaxy', zsource=0.5, fname='Galaxy-z=0.5.txt'

  sdss_rv_table, source='Fstar', zsource=0.0, fname='Fstar-z=0.txt'

  sdss_rv_table, source='QSO', zsource=0.0, fname='QSO-z=0.0.txt'
  sdss_rv_table, source='QSO', zsource=0.5, fname='QSO-z=0.5.txt'
  sdss_rv_table, source='QSO', zsource=1.0, fname='QSO-z=1.0.txt'
  sdss_rv_table, source='QSO', zsource=1.5, fname='QSO-z=1.5.txt'
  sdss_rv_table, source='QSO', zsource=2.0, fname='QSO-z=2.0.txt'
  sdss_rv_table, source='QSO', zsource=2.5, fname='QSO-z=2.5.txt'
  sdss_rv_table, source='QSO', zsource=3.0, fname='QSO-z=3.0.txt'
  sdss_rv_table, source='QSO', zsource=3.5, fname='QSO-z=3.5.txt'
  sdss_rv_table, source='QSO', zsource=4.0, fname='QSO-z=4.0.txt'
  sdss_rv_table, source='QSO', zsource=4.5, fname='QSO-z=4.5.txt'
  sdss_rv_table, source='QSO', zsource=5.0, fname='QSO-z=5.0.txt'
  sdss_rv_table, source='QSO', zsource=5.5, fname='QSO-z=5.5.txt'
  sdss_rv_table, source='QSO', zsource=6.0, fname='QSO-z=6.0.txt'
  sdss_rv_table, source='QSO', zsource=6.5, fname='QSO-z=6.5.txt'





  return
end
