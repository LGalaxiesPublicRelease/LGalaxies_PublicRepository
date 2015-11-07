pro start_plot,thick=thick,size=size
; Defines plotting symbol/character size and thickness
; Also loads a simple colour table for line graphics

; Check arguments and set defaults
if (n_elements(thick) eq 0) then begin
    thick=1.
    if !d.name eq 'X' then thick=1
    if !d.name eq 'PS' then thick=3
endif
if (n_elements(size) eq 0) then size=1.4

; Load simple colour table for line graphics
red =   [0,1,1,0,0,0,1,1,0.5,0,0,0,0.5,0.5]
green = [0,1,0,1,0,1,0,1,0,0.5,0,0.5,0,0.5]
blue =  [0,1,0,0,1,1,1,0,0,0,0.5,0.5,0.5,0]
tvlct, 255*red,255*green,255*blue
; white background; black foreground
!p.background=1
!p.color=0

; Redfine system variables
!p.thick=thick
!p.charthick=thick
!x.thick=thick
!y.thick=thick
!z.thick=thick
!p.symsize=size
!p.charsize=size

end

