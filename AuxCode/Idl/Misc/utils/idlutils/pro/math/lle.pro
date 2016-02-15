;+
; NAME:
;    lle
; PURPOSE: (one line)
;    Perform local linear embedding
; DESCRIPTION:
;    Uses Sam Roweis's local linear embedding technique to reduce the 
;    dimensionality of a data set.
; CATEGORY:
;    Mathematical
; CALLING SEQUENCE:
;    lle, data, ncoords, coords, weights=weights
; INPUTS:
;    data - [p,N] data to be reduced
;    ncoords - number of output dimensions desired
; OUTPUTS: 
;    coords - [ncoords,N] embedding coordinates
; OPTIONAL OUTPUTS PARAMETERS:
;    weights - reconstruction weights
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; BUGS:
;    Not completed yet, do not use
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;    2003-05-14 - Written by Michael Blanton (NYU)
;-
pro lle_find_neighbors, data, nneighbors, neighbors

ndim=(size(data))[0]
if(ndim eq 1) then $
  p=1 $
else $
  p=(size(data,/dimensions))[0]
ndata=n_elements(data)/p

neighbors=lonarr(nneighbors,ndata)

for i=0L, ndata-1L do begin
    dist2=total((data-data[*,i]#replicate(1.,ndata))^2,1)
    sortdist2=sort(dist2)
;   assumes i is one of first nneighbor neighbors
    n=0
    splog,i
    for j=0L, nneighbors do begin
        if(i ne sortdist2[j]) then begin
            neighbors[n,i]=sortdist2[j]
            n=n+1
        endif
    endfor
endfor

end
;
pro lle_reconstruct_weights, data, neighbors, weights

ndim=(size(data))[0]
if(ndim eq 1) then $
  p=1 $
else $
  p=(size(data,/dimensions))[0]
ndata=n_elements(data)/p
nneighbors=n_elements(neighbors)/ndata

weights=dblarr(nneighbors,ndata)
for i=0L, ndata-1L do begin
    zz=data[*,neighbors[*,i]]
    for j=0L, p-1L do $
      zz[j,*]=zz[j,*]-data[j,i]
    covar=transpose(zz)#zz
    covar=covar+1.e-3*trace(covar)*identity(nneighbors)
    ww=invert(covar)#replicate(1.,nneighbors)
    weights[*,i]=ww/total(ww,/double)
endfor

end
;                                
pro lle_embedding_coords, neighbors, weights, ncoords, coords

ndim=(size(neighbors))[0]
if(ndim eq 1) then $
  nneighbors=1 $
else $
  nneighbors=(size(neighbors,/dimensions))[0]
ndata=n_elements(neighbors)/nneighbors

; create the sparse matrix
splog,'   --- creating sparse matrix'
count=(nneighbors+1L)*ndata
imwsm={n:ndata, $
       m:ndata, $
       vals:dblarr(count), $
       nindx:lonarr(count), $
       mindx:lonarr(count)}
for i=0L, ndata-1L do begin
    imwsm.nindx[i*(nneighbors+1L):(i+1L)*(nneighbors+1L)-1L]=i
    imwsm.mindx[i*(nneighbors+1L)]=i
    imwsm.mindx[i*(nneighbors+1L)+1:(i+1L)*(nneighbors+1L)-1L]= $
      neighbors[*,i]
    imwsm.vals[i*(nneighbors+1L)]=1.
    imwsm.vals[i*(nneighbors+1L)+1:(i+1L)*(nneighbors+1L)-1L]= $
      -weights[*,i]
endfor
isort=sort(imwsm.mindx*imwsm.n+imwsm.nindx)
imwsm.nindx=imwsm.nindx[isort]
imwsm.mindx=imwsm.mindx[isort]
imwsm.vals=imwsm.vals[isort]

; square it
splog,'   --- squaring sparse matrix'
msm=lle_sm_mult(lle_sm_transpose(imwsm),imwsm)

; get last few eigenvectors
splog,'   --- getting last few evecs'

; 1. get first eigenvalue
splog,'        --- getting first evec'
em_pca_sparse_matrix,msm,1,eigenvec,maxiter=1000,tol=1.e-6,/verbose
esm=lle_sm_transpose(lle_sm(eigenvec))
evalhigh=(lle_sm_full(lle_sm_mult(msm,esm)))[0]/eigenvec[0]
evalsafe=evalhigh*1.05

; WORKS FINE UP TO THIS --- gets eigenvectors WRONG
; 2. now get the flipped spectrum matrix
;    (take advantage of the fact that we *know* the diagonal is
;    filled)
splog,'        --- flipping matrix'
rsm=msm
rsm.vals=-rsm.vals
indx=where(rsm.nindx eq rsm.mindx)
rsm.vals[indx]=rsm.vals[indx]+evalsafe

; 3. get the first few eigenvectors of the flipped matrix
;    these are the SAME as the last few eigenvectors of the 
;    original matrix
splog,'        --- getting last few evecs'
em_pca_sparse_matrix,rsm,ncoords+1,eigenvec,maxiter=1000,tol=1.e-6,/verbose
coords=eigenvec[*,1:ncoords]

end
;                                ;
pro lle, data, ncoords, coords, weights=weights

if(NOT keyword_set(nneighbors)) then nneighbors=40

; check args
if (n_params() lt 1) then begin
    print, 'Syntax - lle, data, ncoords, coords [, weights= ]'
    return
endif

; check dimensions
ndim=(size(data))[0]
if(ndim eq 1) then $
  p=1 $
else $
  p=(size(data,/dimensions))[0]
n=n_elements(data)/p
if(ncoords gt p) then begin
    splog, 'ncoords must be less than or equal to p!'
    splog, 'ncoords= '+string(ncoords)
    splog, 'p= '+string(p)
    return
endif

splog,'finding neighbors'
lle_find_neighbors, data, nneighbors, neighbors

splog,'reconstructing weights'
lle_reconstruct_weights, data, neighbors, weights

splog,'calculating embedding coords'
lle_embedding_coords, neighbors, weights, ncoords, coords

return

end
