;+
; NAME:
;   VWPCA_GAPPY
;
; AUTHOR:
;   Vivienne Wild <wild@iap.fr>
;
; PURPOSE:
;   Robust Principal Component Analysis projection
;   
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;   PCs = VWPCA_GAPPY (Flux, Error, Evecs, Mean, [COV=,
;                /CONSERVE_MEMORY,/RECONSTRUCT])
;
; INPUTS:
;   Flux  = The NORMALISED data array (npix,ngal)
;   Error = The 1-sigma error array (npix,ngal) =0 for bad pixel
;   Evecs = The eigenvector array (npix,nvec)
;   Mean  = The mean array (npix)
;
; OUTPUTS:
;   PCs   = Principal Component Amplitudes (nvec,ngal)
;
; OPTIONAL OUTPUTS:
;   cov   = Covariance Matrix (used to calculate errors) 
; 
; KEYWORD PARAMETERS:
;   RECONSTRUCT = fill in gaps in flux array with PCA reconstruction
;                 and recalculate mean array (useful
;                 for performing iterative calculation of eigenspectra)
;   CONSERVE_MEMORY  = perform calculations conserving memory instead
;                      of IDL optimised matrix calculations (slower) 
; 
; ADDITIONAL PROGRAMS REQUIRED: 
;   IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;
; REFERENCES:
;
;   [1] Connolly & Szalay (1999, AJ, 117, 2052)
;   http://www.journals.uchicago.edu/AJ/journal/issues/v117n5/980466/980466.html
;
; MODIFICATION HISTORY:
;
;   2005 First implementation in IDL (V.Wild)
;   2007-03-13  Rewritten and properly commented (V.Wild)
;   2007-05-13  Corrected calculation of errors (V.Wild)
;-
;****************************************************************************************;
;  Copyright (c) 2005, Vivienne Wild                                                     ;
;                                                                                        ;
;  Permission to use, copy, modify, and/or distribute this software for any              ;
;  purpose with or without fee is hereby granted, provided that the above                ;
;  copyright notice and this permission notice appear in all copies.                     ;
;                                                                                        ;
;  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES              ;
;  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF                      ;
;  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR               ;
;  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES                ;
;  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN                 ;
;  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF               ;
;  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.                        ;
;****************************************************************************************;



FUNCTION VWPCA_GAPPY, flux, error, espec, mean, cov=cov, reconstruct=reconstruct, conserve_memory=slow,silent=silent

;;-- Return to caller on error.
On_Error, 2

;;-- All inputs must be given

if N_ELEMENTS(flux) eq 0 or N_ELEMENTS(error) eq 0 or $
  N_ELEMENTS(espec) eq 0 or N_ELEMENTS(mean) eq 0 then $
  message, 'USAGE: PCs = VWPCA_GAPPY(flux,error,espec,mean,[cov=,recon=])'

tmp = (SIZE(espec,/dim))        ;number of eigenvectors
if n_elements(tmp) eq 2 then nrecon = tmp[1] else nrecon=1 
nbin = (SIZE(espec,/dim))[0]    ;number of data points
tmp = (SIZE(flux,/dim))         ;number of observations to project
if n_elements(tmp) eq 2 then ngal = tmp[1] else ngal=1

;;-- Check for dimensionality mismatches

if (SIZE(flux,/dim))[0] ne nbin then message,'VWPCA_GAPPY: "flux" must have same dimension as eigenvectors'
if (SIZE(error,/dim))[0] ne nbin then message,'VWPCA_GAPPY: "error" must have same dimension as eigenvectors'
if (SIZE(mean,/dim))[0] ne nbin then message,'VWPCA_GAPPY: "mean" must have same dimension as eigenvectors'


;;------------------------------------------------------------------
;;*** project each galaxy in turn
;;------------------------------------------------------------------

pcs = dblarr(nrecon,ngal)
if ARG_PRESENT(cov) then cov = dblarr(nrecon,nrecon,ngal)

for j=0L,ngal-1 do begin

    if NOT(KEYWORD_SET(silent)) then if j/100 eq j/100. then $
      print, 'VWPCA_GAPPY: processing spectrum: ',j

;;-- Calculate weighting array from 1-sigma error array
;; If all bins have error=0 continue to next spectrum

    weight = dblarr(nbin)
    ind = where(error[*,j] ne 0.) 
    if ind[0] ne -1 then weight[ind] = 1./(double(error[ind,j])^2) else begin
        if NOT(KEYWORD_SET(silent)) then $
          print,'VWPCA_GAPPY: error array problem in spectrum (setting pcs=0)',j
        continue
    endelse

    ind = where(finite(weight) eq 0)
    if ind[0] ne -1 then begin
        if NOT(KEYWORD_SET(silent)) then $
          print,'VWPCA_GAPPY: error array problem in spectrum (setting pcs=0)',j
        continue
    endif

;;-- Subtract the mean from the data

    data = double(flux[*,j] - mean)

;;-- Calculate the weighted eigenvectors, multiplied by the eigenvectors (eq. 4-5 [1])

    if nrecon gt 1 then begin
       
        if KEYWORD_SET(slow) then begin
            for l=0,nbin-1 do $
              M = M + weight[l]*MATRIX_MULTIPLY(espec[l,*],espec[l,*], /atranspose)
        endif else begin
            M = MATRIX_MULTIPLY(double(espec)*rebin(weight,nbin,nrecon),$
                                double(espec), /atranspose)
        endelse

;;--- Calculate the weighted data array, multiplied by the eigenvectors (eq. 4-5 [1])
        
        F = MATRIX_MULTIPLY(weight*data,double(espec),/atranspose)
        
        
;;-- Solve for Principal Component Amplitudes (eq. 5 [1])        
        
        Minv = invert(M,status,/double)
        
;;-- If status = 1 or 2 there's a definite problem

        if status ne 0 then begin
            if NOT(KEYWORD_SET(silent)) then $
              print, 'VWPCA_GAPPY: problem with matrix inversion (setting pcs=0)', j,status            
            continue            ;pcs will be 0.0 for this file
        endif
        
        pcs[*,j] = reform(MATRIX_MULTIPLY(F,Minv))


;;-- Calculate covariance matrix (eq. 6 [1])
;; NB - there seems to be a mistake in [1] here. Using just Minv gives very
;;      consistent answers c.f. ordinary propogation of errors

    if ARG_PRESENT(cov) then cov[*,*,j] = Minv

;;-- If we only have one eigenvector

    endif else begin
        M = total(weight*espec*espec)
        F = total(weight*data*espec)
        pcs[0,j] = F/M
        if ARG_PRESENT(cov) then cov[0,0,j] = total((1./weight)*espec*espec)
    endelse
        


;;-- If reconstruction of flux array required, fill in regions with
;;   weight=0 with PCA reconstruction  

    if KEYWORD_SET(reconstruct) then begin

        bad_pix = WHERE(weight eq 0.,count)
        if count eq 0 then continue
        
        reconstruct = TOTAL(REBIN(TRANSPOSE(pcs[*,j]),N_ELEMENTS(bad_pix),nrecon)*espec[bad_pix,*],2)
        reconstruct += mean[bad_pix]
        flux[bad_pix,j] = reconstruct

    endif

endfor

;;-- Calculate new mean if flux array has been updated

if KEYWORD_SET(reconstruct) then mean = TOTAL(flux,2,/double)/double(nbin)

return,pcs

END
