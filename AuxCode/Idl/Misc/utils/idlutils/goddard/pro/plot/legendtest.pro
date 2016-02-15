
;+
; NAME:
;	LEGENDTEST
; PURPOSE:
;	Demo program to show capabilities of  the legend procedure.
; CALLING SEQUENCE:
;	legendtest
; INPUTS:
;	none
; OPTIONAL INPUTS:
;	none
; KEYWORDS:
;	none
; OUTPUTS:
;	legends of note
; COMMON BLOCKS:
;	none
; SIDE EFFECTS:
;	Sets !20 font to symbol if PostScript and !p.font=0.
; RESTRICTIONS:
;	With the vectorfont test, you'll get different results for PostScript
;	depending on the value of !p.font.
; MODIFICATION HISTORY:
;	write, 27 Aug 92, F.K.Knight (knight@ll.mit.edu)
;	add test of /left,/right,/top,/bottom keywords, 21 June 93, FKK
;	update based on recent changes to legend, 7 Feb 94, FKK
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
pro legendtest
if (!d.name eq 'PS') and (!p.font eq 0) then device,/Symbol,font_index=20
items = ['diamond','asterisk','square']
explanation = ['The legend procedure annotates plots---' $
  ,'  either using text alone,' $
  ,'  or text with plot symbols, lines, and special characters.' $
  ,'The following are some examples.' $
  ,'Hit return to continue.']
psym = [4,2,6]
lineitems = ['solid','dotted','DASHED']
linestyle = [0,1,2]
citems = 'color '+strtrim(string(indgen(8)),2)
colors = 15*indgen(8)+50
z =	['legend,explanation,charsize=1.5' $
	,'legend,items,psym=[4,2,6]' $
	,'plot,findgen(10) & legend,items,psym=[4,2,6] & legend,items,psym=[4,2,6],/bottom,/right' $
	,'legend,lineitems,linestyle=linestyle,/right,/bottom' $
	,'legend,items,psym=psym,/horizontal,chars=1.5	; horizontal format' $
	,'legend,[items,lineitems],psym=[psym,0,0,0],line=[0,0,0,linestyle],/center,box=0		; sans border' $
	,'legend,items,psym=psym,margin=1,spacing=2,char=2,delimiter="=",/top,/center; delimiter & larger margin' $
	,'legend,lineitems,line=linestyle,pos=[.3,.5],/norm,char=3,number=4	; position of legend' $
	,'legend,items,psym=-psym,number=2,line=linestyle,/right; plot two symbols, not one' $
	,'legend,citems,/fill,psym=8+intarr(8),colors=colors,char=2; 8 filled squares' $
	,'legend,[citems[0:4],lineitems],/fill,psym=[8+intarr(5),0*psym],line=[intarr(5),linestyle],colors=colors,char=2,text=colors' $
	,"legend,['Absurd','Sun Lover','Lucky Lady','Fishtail Palm'],vector=['ab!9r!3','!9nu!3','!9Wf!3','!9cN!20K!3'],charsize=2,/pos,psp=3"$
	]
prompt = 'Hit return to continue:'
for i = 0,n_elements(z) - 1 do begin
  erase
  stat = execute(z[i])
  xyouts,.01,.15,'COMMAND TO MAKE LEGEND:',charsize=1.7,/norm
  xyouts,.01,.05,z[i],/norm,charsize=1.2
  print,'Command: ',z[i]
  print,prompt,format='($,a)'
  a = get_kbrd(1)
  print
  endfor
;stop
erase
!p.charsize=2
c1_items = ['Plus','Asterisk','Period','Diamond','Triangle','Square','X']
c1_psym = indgen(7)+1
c2_items = ['Solid','Dotted','Dashed','Dash Dot','Dash Dot Dot Dot','Long Dashes']
c2_line = indgen(6)
legend,c1_items,psym=c1_psym,corners=c1,box=0
legend,c2_items,line=c2_line,corners=c2,box=0,pos=[c1[2],c1[3]],/norm
c = [c1[0]<c2[0],c1[1]<c2[1],c1[2]>c2[2],c1[3]>c2[3]]
plots,[c[0],c[0],c[2],c[2],c[0]],[c[1],c[3],c[3],c[1],c[1]],/norm
!p.charsize=0
xyouts,.01,.05,'Multiple columns---type "legend,/help" for details.',/norm,charsize=1.2
return
end

