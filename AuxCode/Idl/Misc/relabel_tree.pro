;=========================================================================
;
pro relabel_tree
;  procedure to read in L-galaxies galaxy tree data
;  and relabel the tree pointers
;
;  Note: to do all files simultaneously, rather than one at a time as
;  below, would need to change TID to include the file number also (as
;  trees are not unique between files).  This is easily done as
;  TID=Galstruct.TreeRootID/1000000 mod 1000000
;  FID=Galstruct.TreeRootID/1000000000000
;  FTID=Galstruct.TreeRootID/1000000
;
;-------------------------------------------------------------------------

compile_opt idl2

; Define which files you want to read in
DirName = '../output/'
@'../Output/galaxytree_template.pro'
FileName = 'SA_galtree'
ModelName = DirName + FileName
;FirstFile = 0
;LastFile = 511
FirstFile = 287
LastFile = 287

;-------------------------------------------------------------------------------

; Loop over files

for File=FirstFile,LastFile do begin

read_lgaltree,ModelName,Template,GalStruct $
         ,StructSize=StructSize $
         ,FirstFile=File,LastFile=File ;,/swap_endian

ngal=n_elements(Galstruct)

OldPos=lonarr(ngal) ; Old ordering
BushID=lonarr(ngal)
GalID=lon64arr(ngal)
FirstProgGal=lon64arr(ngal)
NextProgGAl=lon64arr(ngal)
LastProgGal=lon64arr(ngal)
DescendantGal=lon64arr(ngal)
MainLeafID=lon64arr(ngal)
TreeRootID=lon64arr(ngal)
FoFCentralGal=lon64arr(ngal)

;------------------------------------------------------------------------------

; All galaxies in a tree have the same value of treerootid and they
; are depth-first ordered (I think) within a tree.

; The tree information is stored in decimal places 7-12
; NO IT IS NOT.  This looks like it might be bush info.
; But anyway, the important thing is that galaxy IDs run
; consecutively upwards from these IDs, then jump to the next one.
TID=Galstruct.TreeRootID/1000000 mod 1000000
; i_gal is new galaxy order
; i_tree is new BUSH label
; off is offest that needs to be applied to old ids
; Note that these offsets should not be applied to -1 entries
; - this is corrected below
i_gal=-1LL
; Loop over BUSHES
for i_tree=min(TID),max(TID) do begin
; galaxies in BUSH
   index=where(TID eq i_tree)
; set root of BUSH to nextavailable galaxy ID
   off=i_gal+1LL-((Galstruct[index[0]].TreeRootID/1000000)*1000000)
; relabel
   for l=0,n_elements(index)-1 do begin
      i_gal=i_gal+1
      OldPos[i_gal]=index[l]
      BushID[i_gal]=i_tree
      GalID[i_gal]=galstruct[index[l]].GalID+off
      FirstProgGal[i_gal]=galstruct[index[l]].FirstProgGal+off
      NextProgGal[i_gal]=galstruct[index[l]].NextProgGal+off
      LastProgGal[i_gal]=galstruct[index[l]].LastProgGal+off
      DescendantGal[i_gal]=galstruct[index[l]].DescendantGal+off
      MainLeafID[i_gal]=galstruct[index[l]].MainLeafID+off
      TreeRootID[i_gal]=galstruct[index[l]].TreeRootID+off
      FoFCentralGal[i_gal]=galstruct[index[l]].FoFCentralGal+off
   endfor
   if i_tree mod 100 eq 0 then $
      print,'Tree=',i_tree,', totngal=',i_gal
endfor
print,'Done: Tree=',i_tree,', totngal=',i_gal

; Reset -1 entries
GalID=GalID > (-1)
FirstProgGal=FirstProgGal > (-1)
NextProgGal=NextProgGal > (-1)
LastProgGal=LastProgGal > (-1)
DescendantGal=DescendantGal > (-1)
MainLeafID=MainLeafID > (-1)
TreeRootID=TreeRootID > (-1)
FoFCentralGal=FoFCentralGal > (-1)
print,'-1s reset'

; Now lets re-order the TREES into depth-first order
; I think that this is as simple as ordering in GalID order
index=sort(GalID)
OldPos=OldPos[index]
BushID=BushID[index]
GalID=GalID[index]
FirstProgGal=FirstProgGal[index]
NextProgGal=NextProgGal[index]
LastProgGal=LastProgGal[index]
DescendantGal=DescendantGal[index]
MainLeafID=MainLeafID[index]
TreeRootID=TreeRootID[index]
FoFCentralGal=FoFCentralGal[index]
print,'Re-ordered'

; Create a new galaxy structure with the sorted information and new
; galaxy labelling.
; For ease of use of galaxy template, do not alter entries
GalStructNew=GalStruct[OldPos]
GalstructNew.FirstProgGal=FirstProgGal
GalstructNew.NextProgGal=NextProgGal
GalstructNew.LastProgGal=LastProgGal
GalstructNew.DescendantGal=DescendantGal
GalstructNew.MainLeafID=MainLeafID
GalstructNew.TreeRootID=TreeRootID
GalstructNew.FofCentralGal=FofCentralGal
; Extract the File number
GalstructNew.FileTreeNr=GalstructNew.FileTreeNr/1000000000000LL
print,'New structure created'

; Write out 
write_lgaltree,ModelName+'new',GalStructNew,StructSize;,swap_endian=swap_endian
print,'Written out new data'

; End loop over files
endfor

end
