;+
; NAME:
;   VWPCA_NORMGAPPY
;
; AUTHOR:
;   Vivienne Wild <wild@iap.fr>
;   Gerard Lemson <lemson@mpe.mpg.de>
;
; PURPOSE:
;   Robust Principal Component Analysis projection, including
;   robust estimation of normalisation of input data
;   
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;   PCs = VWPCA_NORMGAPPY (Data, Error, Evecs, Mean, [/RECONSTRUCT,
;                    /CONSERVE_MEMORY, Norm=, COV=])
;
; INPUTS:
;   Data  = The 2-d data matrix (nbin,nobj), arbitrary normalisation allowed
;   Error = The 1-sigma error matrix (nbin,nobj)
;   Evecs = The eigenvector matrix (nbin,nvec)
;   Mean  = The mean vector (nbin). This is recalculated if
;           /reconstruct set.
;
; KEYWORD PARAMETERS:
;   /RECONSTRUCT = fill in gaps in data matrix with PCA reconstruction
;                 and recalculate mean array (useful
;                 for performing iterative calculation of eigenspectra)
;   /CONSERVE_MEMORY  = perform calculations conserving memory instead
;                      of IDL optimised matrix calculations (slower) 
;
; OUTPUTS:
;   PCs   = Principal Component Amplitudes (nvec,nobj)
;
; OPTIONAL OUTPUTS:
;   Norm  = normalisation applied to data vectors (nobj)
;   Cov   = covariance matrix of PCs
; 
; ADDITIONAL PROGRAMS REQUIRED: 
;   IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;
; REFERENCES:
;
;   [1] Connolly & Szalay (1999, AJ, 117, 2052)
;   http://www.journals.uchicago.edu/AJ/journal/issues/v117n5/980466/980466.html
;   [2] Lemson, "Normalized gappy PCA projection"
;
; MODIFICATION HISTORY:
;
;   2007-03-16  First implementation in IDL (V. Wild)
;   2008-12-11  Tidy up for public release (V. Wild) 
;-
;****************************************************************************************;
;  Copyright (c) 2007, Gerard Lemson, Vivienne Wild                                                     ;
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

FUNCTION VWPCA_NORMGAPPY, data, error, espec, mean, norm=norm, cov=cov,reconstruct=reconstruct, conserve_memory=slow

;;-- Return to caller on error.
On_Error, 2

;;-- All inputs must be given

if N_ELEMENTS(data) eq 0 or N_ELEMENTS(error) eq 0 or $
  N_ELEMENTS(espec) eq 0 or N_ELEMENTS(mean) eq 0 then $
  message, 'USAGE: PCs = VWPCA_NORMGAPPY(data,error,espec,mean,recon=,norm=,conserve_memory=])'

tmp = (SIZE(espec,/dim))        ;number of eigenvectors
if n_elements(tmp) eq 2 then nrecon = tmp[1] else nrecon=1 
nbin = (SIZE(espec,/dim))[0]    ;number of data points
tmp = (SIZE(data,/dim))         ;number of observations to project
if n_elements(tmp) eq 2 then nobj = tmp[1] else nobj=1


;;-- Check for dimensionality mismatches

if (SIZE(data,/dim))[0] ne nbin then message,'VWPCA_NORMGAPPY: "data" must have same dimension as eigenvectors'
if (SIZE(error,/dim))[0] ne nbin then message,'VWPCA_NORMGAPPY: "error" must have same dimension as eigenvectors'
if (SIZE(mean,/dim))[0] ne nbin then message,'VWPCA_NORMGAPPY: "mean" must have same dimension as eigenvectors'

;;------------------------------------------------------------------
;;*** project each galaxy in turn
;;------------------------------------------------------------------

pcs = dblarr(nrecon,nobj)
norm = dblarr(nobj)
if ARG_PRESENT(cov) then cov = dblarr(nrecon,nrecon,nobj)

for j=0L,nobj-1 do begin

;;-- Calculate weighting array from 1-sigma error array
;;   If all bins have error=0 continue to next spectrum

    weight = dblarr(nbin)
    ind = where(error[*,j] ne 0.) 
    if ind[0] ne -1 then weight[ind] = 1./(double(error[ind,j])^2) else begin
        if NOT(KEYWORD_SET(silent)) then $
          print,'VWPCA_NORMGAPPY: error array problem in spectrum (setting pcs=0)',j
        continue
    endelse
    
    ind = where(finite(weight) eq 0)
    if ind[0] ne -1 then begin
        if NOT(KEYWORD_SET(silent)) then $
          print,'VWPCA_NORMGAPPY: error array problem in spectrum (setting pcs=0)',j
        continue
    endif

    data_j = data[*,j] 

;;-- solve partial chi^2/partial N = 0
    Fpr = total(weight*data_j*mean)    ;eq 4 [2]
    Mpr = total(weight*mean*mean)      ;eq 5 [2]
    E = total(rebin(weight*mean,nbin,nrecon)*espec,1) ;eq 6 [2]
    

;;-- Calculate the weighted eigenvectors, multiplied by the eigenvectors
;; (eq. 4-5 [1], 7-8 [2])
    if nrecon gt 1 then begin
       
        if KEYWORD_SET(slow) then begin
            for l=0,nbin-1 do $
              M = M + weight[l]*MATRIX_MULTIPLY(espec[l,*],espec[l,*], /atranspose)
        endif else begin
            M = MATRIX_MULTIPLY(double(espec)*rebin(weight,nbin,nrecon),$
                                double(espec), /atranspose)
        endelse
        
;;-- Calculate the weighted data array, multiplied by the eigenvectors (eq. 4-5 [1])
        F = reform(MATRIX_MULTIPLY(weight*data_j,double(espec),/atranspose))
        
;;-- Calculate new M matrix, this time accounting for unknown normalisation (eq. 11 [2])
        Mnew = Fpr*M - E#F

;;-- Calulate new F matrix, accounting for unknown normalisation
        Fnew = Mpr*F - Fpr*E
        
;;------------------------------------------------------------------
;;*** Solve for Principal Component Amplitudes (eq. 5 [1])        
;;------------------------------------------------------------------
        Minv = invert(Mnew,status,/double)
        
;;-- If status = 1 or 2 there's a definite problem

        if status ne 0 then begin
            if NOT(KEYWORD_SET(silent)) then $
              print, 'VWPCA_NORMGAPPY: problem with matrix inversion (setting pcs=0)', j,status
            continue            ;pcs will be 0.0 for this file
        endif
        
        pcs[*,j] = reform(MATRIX_MULTIPLY(Fnew,Minv))
        norm[j] = Fpr / (Mpr + total(pcs[*,j]*E))

;;-- Calculate covariance matrix, currently just on PCs using ordinary gappy PCA

        if ARG_PRESENT(cov) then begin
            M_gappy = MATRIX_MULTIPLY(double(espec)*rebin(weight*norm[j]^2,nbin,nrecon),$
                                      double(espec), /atranspose)
            cov[*,*,j] = invert(M_gappy,status,/double)
        endif

    endif else begin
;
;-- If we only have one eigenvector

        M = total(weight*espec*espec)
        F = total(weight*data_j*espec)
        Mnew = M*Fpr - E*F
        Fnew = Mpr*F - E*Fpr
        pcs[0,j] = Fnew/Mnew
        norm[j] = Fpr / (Mpr + pcs[0,j]*E)
    endelse
        
;;-- If reconstruction of data array required, fill in regions with weight=0 with PCA reconstruction  

    if KEYWORD_SET(reconstruct) then begin

        bad_pix = WHERE(weight eq 0.,count)
        if count eq 0 then continue
        
        reconstruct = TOTAL(REBIN(TRANSPOSE(pcs[*,j]),N_ELEMENTS(bad_pix),nrecon)*espec[bad_pix,*],2)
        reconstruct += mean[bad_pix]
        data[bad_pix,j] = reconstruct

    endif

endfor

;;-- Calculate new mean if data array has been updated

if KEYWORD_SET(reconstruct) then mean = TOTAL(data,2,/double)/double(nbin)

return,pcs

END
