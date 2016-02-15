;+
; NAME:
;   atv_joe_writeps
;
; PURPOSE:
;   Write a PostScript file of the current ATV display w/out using a widget
;
; CALLING SEQUENCE:
;   atv_joe_writeps, filename, [ _EXTRA ]
;
; INPUTS:
;   filename   - Name of PostScript file
;
; OPTIONAL INPUTS:
;   aspect     - retain aspect ratio in .ps file
;   _EXTRA     - Optional keywords for DEVICE
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine allows one to send the current contents of the ATV window
;   to a PostScript file without using a widget.  This makes it convenient
;   for using ATV in a batch mode to make plots (though the ATV window will
;   still pop up on your terminal).
;
;   Note that there are a number of defaults set in the call to DEVICE
;   that can be over-written in the call to this routine.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   08-May-2003  Written by J. Hennawi and D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro atv_joe_writeps, filename, aspect=aspect, _EXTRA=KeywordsForTV $
                     , CLOBBER = CLOBBER

   common atv_state
   common atv_images
   common atv_color

   widget_control, /hourglass

   view_min = round(state.centerpix - $
    (0.5 * state.draw_window_size / state.zoom_factor))
   view_max = round(view_min + state.draw_window_size / state.zoom_factor)

   xsize = (state.draw_window_size[0] / state.zoom_factor) > $
    (view_max[0] - view_min[0] + 1)
   ysize = (state.draw_window_size[1] / state.zoom_factor) > $
    (view_max[1] - view_min[1] + 1)

   viewaspect = float(ysize) / float(xsize)
   fname = strcompress(state.current_dir + 'atv.ps', /remove_all)

   tvlct, rr, gg, bb, 8, /get

   if (NOT keyword_set(filename)) then return
   forminfo = create_struct( $
    'BITS_PER_PIXEL', 8, $
    'COLOR'         , 1, $
    'ENCAPSULATED'  , 1, $
    'FILENAME'      , filename, $
    'FONT_SIZE'     , 12, $
    'INCHES'        , 1, $
    'ISOLATIN1'     , 0, $
    'PREVIEW'       , 0, $
    'TT_FONT'       , 0, $
    'XOFFSET'       , 1.00, $
    'XSIZE'         , 6.50, $
    'YOFFSET'       , 2.25, $
    'YSIZE'         , 6.50, $
    'PORTRAIT'      , 1, $
    'LANDSCAPE'     , 0, $
    'HELVETICA'     , 1, $
    'BOLD'          , 0, $
    'BOOK'          , 0, $
    'DEMI'          , 0, $
    'ITALIC'        , 0, $
    'LIGHT'         , 0, $
    'MEDIUM'        , 0, $
    'NARROW'        , 0, $
    'OBLIQUE'       , 0 )
   if (keyword_set(KeywordsForTV)) then begin
      tags_form = tag_names(forminfo)
      tags_in = tag_names(KeywordsForTV)
      for i=0L, n_elements(tags_in)-1 do begin
         j = (where(tags_in[i] EQ tags_form, ct))[0]
         if (ct GT 0) then forminfo.(j) = KeywordsForTV.(i)
      endfor
   endif
   if(keyword_set(aspect)) then forminfo.ysize=forminfo.xsize*viewaspect

   tvlct, rr, gg, bb, 8

   tmp_result = findfile(forminfo.filename, count = nfiles)

   result = ''
   if (nfiles GT 0) AND NOT KEYWORD_SET(CLOBBER) then begin
      mesg = strarr(2)
      mesg[0] = 'Overwrite existing file:'
      tmp_string = strmid(forminfo.filename, $
       strpos(forminfo.filename, '/', /reverse_search) + 1)
      mesg[1] = strcompress(tmp_string + '?', /remove_all)
      result = dialog_message(mesg, /default_no, $
       dialog_parent=state.base_id, /question)                 
   endif

   if (strupcase(result) EQ 'NO') then return
    
   widget_control, /hourglass

   screen_device = !d.name

   ; In 8-bit mode, the screen color table will have fewer than 256
   ; colors.  Stretch out the existing color table to 256 colors for the
   ; postscript plot.

   set_plot, 'ps'

   device, _EXTRA=forminfo

   tvlct, rr, gg, bb, 8, /get

   rn = congrid(rr, 248)
   gn = congrid(gg, 248)
   bn = congrid(bb, 248)

   tvlct, temporary(rn), temporary(gn), temporary(bn), 8

   ; Make a full-resolution version of the display image, accounting for
   ; scalable pixels in the postscript output

   newdisplay = bytarr(xsize, ysize)

   startpos = abs(round(state.offset) < 0)

   view_min = (0 > view_min < (state.image_size - 1)) 
   view_max = (0 > view_max < (state.image_size - 1)) 

   dimage = bytscl(scaled_image[view_min[0]:view_max[0], $
    view_min[1]:view_max[1]], $
    top=247, min=8, max=(!d.table_size-1)) + 8

   newdisplay[startpos[0], startpos[1]] = temporary(dimage)

   ; If there's blank space around the image border, keep it black
   tv, newdisplay
   atv_plotall

   if (state.frame EQ 1) then begin    ; put frame around image
      plot, [0], [0], /nodata, position=[0,0,1,1], $
       xrange=[0,1], yrange=[0,1], xstyle=5, ystyle=5, /noerase
      boxx = [0,0,1,1,0,0]
      boxy = [0,1,1,0,0,1]
      oplot, boxx, boxy, color=0, thick=state.framethick
   endif

   tvlct, temporary(rr), temporary(gg), temporary(bb), 8


   device, /close
   set_plot, screen_device

   return
end
;------------------------------------------------------------------------------
