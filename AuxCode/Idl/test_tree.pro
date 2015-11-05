;=========================================================================
;
pro test_tree
;  procedure to read in L-galaxies galaxy tree data
;  and test the tree pointers
;
;-------------------------------------------------------------------------

compile_opt idl2

; Define which files you want to read in
DirName = '../output/'
@'../output/galaxytree_template.pro'
FileName = 'SA_galtree_new'
ModelName = DirName + FileName
;FirstFile = 0
;LastFile = 511
FirstFile = 287
LastFile = 287
NumSnap=64

read_lgaltree,ModelName,Template,GalStruct $
         ,FirstFile=FirstFile,LastFile=LastFile ;,/swap_endian

;-------------------------------------------------------------------------------

GalID=GalStruct.GalID
Type=GalStruct.Type
First=Galstruct.FirstProgGal
Next=Galstruct.NextProgGal
Last=Galstruct.LastProgGal
FoFGal=GalStruct.FoFCentralGal
Descendant=GalStruct.DescendantGal
SnapNum=Galstruct.SnapNum
MMSubID=Galstruct.MMSubID
Mass=Galstruct.StellarMass
ICM=Galstruct.ICM

ngal=n_elements(Galstruct)
FoFGalList=lonarr(ngal)

;------------------------------------------------------------------------------

; The following code gives the same list of mis-matches whether I use
; FoFCentralGal or MMSubID.
; However the mismatches all come from galaxies that have missed a
; snap.

;; ifg=-1
;; for i=0L,ngal-1 do begin
;;    j=First[i]
;;    if j gt -1 then begin
;;       ;fg=FoFGal[j]
;;       fg=MMSubID[j]
;;       ifg=ifg+1
;;       FoFGalList[ifg]=fg
;;       if Descendant[j] ne i then message,'Descendant ne GalID'
;;       while Next[j] gt -1 do begin
;;          j=Next[j]
;;          ;if FoFGal[j] ne fg then print,'FoFGal ne fg, i=',long(i),$
;;          if MMSubID[j] ne fg then print,'MMSubID ne fg, i=',long(i),$
;;              ', j=',long(j),', Delta_Snap=',SnapNum[i]-SnapNum[j]
;;          if Descendant[j] ne i then message,'Descendant ne GalID'
;;       endwhile
;;    endif
;; endfor

; Check all FoF central gals are different
; - actually, don't expect them to be because we can have many
;   "First progenitor" galaxies in a FoF group - the tree structure is
;   running over galaxies, not FoF halos.
;FoFGalList=FoFGalList[0:ifg]
;FoFGalList=FoFGalList[sort(FoFGalList)]
;diff=FoFGalList[1:ifg-1]-FoFGalList[0:ifg-2]
;if (where(diff eq 0))[0] ne -1 then message,'Repeated FoFCentralGal'

;; ; Check that FoFgal is equal to galaxy id iff type eq 0 - it is
;; if (where(FoFGal eq lindgen(ngal) and type ne 0))[0] ne -1 then $
;;    print,'Galaxy is a FoFGal and not type 0'
;; if (where(FoFGal ne lindgen(ngal) and type eq 0))[0] ne -1 then $
;;    print,'Galaxy is type 0 and not a FoFGal'

;; ; Check that First galaxy in a FoF group is always the one with the
;; ; lowest id. - it is not
;; if (where(FoFGal gt lindgen(ngal)))[0] ne -1 then $
;;    print,'The first galaxy in a FoF group is not always the central galaxy'

;; ; Check algorithm for creating FoFchains
;; FoFChain=intarr(ngal)-1
;; for i=0,ngal-1 do begin
;;    if Type[i] ne 0 then begin
;;       FoFChain[i]=FoFChain[FoFGal[i]]
;;       FoFChain[FoFGal[i]]=i
;;    endif
;; endfor
;; ; Loop over chains to see if it works - it does
;; Test=intarr(ngal)
;; for i=0,ngal-1 do begin
;;    if Type[i] eq 0 then begin
;;       Test[i]+=1
;;       j=FoFChain[i]
;;       while j ne -1 do begin
;;          Test[j]+=1
;;          j=FoFChain[j]
;;       endwhile
;;    endif
;; endfor
;; if (where(Test ne 1))[0] ne -1 then $
;;    print,'FoFChain test failed'

;; ; Check to see whether the Type 0 is always the most massive galaxy in
;; ; a chain - it is not
;; Countgal=0
;; Counthalo=0
;; for i=0,ngal-1 do begin
;;    if Type[i] eq 0 then begin
;;       j=FoFChain[i]
;;       while j ne -1 do begin
;;          if Mass[j] gt Mass[i] then Countgal+=1
;;          if Galstruct[j].Mvir gt Galstruct[i].Mvir then Counthalo+=1
;;          j=FoFChain[j]
;;       endwhile
;;    endif
;; endfor
;; if Countgal gt 0 then print,'Type 0 galaxy not most massive ',Countgal,' times'
;; if Counthalo gt 0 then print,'Type 0 halo not most massive ',Counthalo,' times'

;; ; Try to reproduce Chris' selection scheme - it does reproduce it
;; ; The missing galaxies are all progenitors of orphaned galaxies.
;; Test=intarr(ngal)
;; for i=0,ngal-1 do begin
;;    if Type[i] eq 0 and SnapNum[i] eq NumSnap-1 then begin
;;       j=i
;;       while j ne -1 do begin
;;          Test[j:Last[j]]+=1
;;          j=FoFChain[j]
;;       endwhile
;;    endif
;; endfor
;; print,'Chris'' method found ',n_elements(where(Test eq 0)),' no times'
;; print,'Chris'' method found ',n_elements(where(Test eq 1)),' once'
;; print,'Chris'' method found ',n_elements(where(Test gt 1)),' twice or more'

; Test to see if a galaxy can ever be less massive than the sum of its
; progenitors - it can, because of the transfer of material to icm.
; Note that there will be cases where both stars are formed and
; disruption occurs; in that case Massdiff will be zero even though
; stars have formed.
Massdiff=Mass
for i=0,ngal-1 do begin
   j=First[i]
   while j ne -1 do begin
      Massdiff[i] = Massdiff[i]-Mass[j]
      j=Next[j]
   endwhile
;   if Massdiff[i] lt -1e-5*Mass[i] then $
;      print,i,type[i],snapnum[i],-Massdiff[i],Mass[i],-Massdiff[i]/Mass[i]
   if Massdiff[i] lt 0. then Massdiff[i]=0.
endfor

;; ; Test to see if a galaxy +icm can ever be less massive than the sum of its
;; ; progenitors - it can, only when a galaxy changes from type 1 to type 2.
;; Massdiff=Mass+ICM
;; for i=0,ngal-1 do begin
;;    j=First[i]
;;    while j ne -1 do begin
;;       Massdiff[i] = Massdiff[i]-Mass[j]-ICM[j]
;;       j=Next[j]
;;    endwhile
;;    if Massdiff[i] lt -1e-5*(Mass[i]+ICM[i]) then $
;;       print,i,type[i],snapnum[i],-Massdiff[i],(Mass[i]+ICM[i]),$
;;             -Massdiff[i]/(Mass[i]+ICM[i])
;; endfor

; Count growth of mass over time and compare to changes in stellar
; mass.
Masstot=dblarr(NumSnap)
dMass=dblarr(NumSnap)
Masslost=dblarr(NumSnap) ; Galaxies that go missing on this snap for ever
Massextra=dblarr(NumSnap) ; Galaxies that are present in this snap but not the next
Massmissing=dblarr(NumSnap) ; Galaxies that are missing from this snap
for i=0,ngal-1 do begin
   snap=SnapNum[i]
   Masstot[snap]+=(Mass[i]+ICM[i])
   dMass[snap]+=Massdiff[i]
   if Descendant[i] eq -1 and snap lt NumSnap-1 then $
      Masslost[snap+1]+=(Mass[i]+ICM[i])
   if Descendant[i] ne -1 then $
      if SnapNum[Descendant[i]] ne snap+1 then $
         Massextra[snap]+=(Mass[i]+ICM[i])
   if First[i] ne -1 then $
      if SnapNum[i] ne SnapNum[First[i]]+1 then $
         Massmissing[snap-1]+=(Mass[First[i]]+ICM[First[i]])
endfor
print,'        Snap,           Mass,     DeltaMass,          dMass,       Corrected'
for i=1,NumSnap-1 do $
   print,i,Masstot[i],Masstot[i]-Masstot[i-1],dMass[i],dMass[i]-Massextra[i-1]+Massmissing[i-1]

stop

end
