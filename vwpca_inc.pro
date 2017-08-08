;+
; NAME:
;   VWPCA_INC
;
; AUTHOR:
;   Tamas Budavari <budavari@jhu.edu>
;   Vivienne Wild <wild@iap.fr>
;
; PURPOSE:
;   Incremental and Robust Princial Component Analysis
;
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;   NewSys = VWPCA_INC (EigSys, X, Alpha, Delta)
;
; INPUTS:
;   EigSys = The eigensystem and mean to update; struct with M, W, U, sig2, Vk, Uk, Qk
;            M = mean (nbin) vector
;            W = eigenvalues (nvec) vector
;            U = eigenvectors (nbin,nvec) array
;            Vk = the total weight of all proceeding objects (dbl)
;            Qk = the total weight*r^2 of all proceeding objects (dbl)
;            Uk = the total number of all proceeding objects (dbl)
;            sig2 = the robust estimate of the sample standard
;                   deviation (dbl)
;   X      = The new data (nbin) vector for including into and updating the eigensystem
;            (Must be pre-normalised) 
;   Alpha  = The control parameter of convergence 
;            (e.g. 1.0 - 1.0/double(niter+1), where niter = number of iterations).
;   Delta   = The control parameter of robustness (in units of standard
;            deviation. Set high for fast model updating, but at risk of
;            accepting outliers). 
;
; OUTPUTS:
;   NewSys = The new eigensystem and mean; struct like EigSys.
;
; KEYWORDS:
;   CONSERVE_MEMORY  = perform calculations conserving memory instead
;                      of IDL optimised matrix calculations (slower) 
;
; ADDITIONAL PROGRAMS REQUIRED: 
;   tag_exist.pro IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;
; REFERENCES:
;
;   [1] Budavari, Wild et al. 2008. "Reliable eigenspectra for new generation surveys", ArXiv 0809:0881
;   [3] Y. Li and J. Xu and L. Morphett and R. Jacobs, 2003 
;   "An integrated algorithm of incremental and robust PCA"
;   in the Proceedings of IEEE International Conference of 
;   Image Processing, 2003 
;   http://citeseer.ist.psu.edu/li03integrated.html
;
; MODIFICATION HISTORY:
;
;   2007-02-21  New IDL implementation and works, TB & VW
;   2007-02-23  First bash at the documentation, TB
;   2007-03-13  Documentation finished, IDL opitimised, keywords
;               "conserve_memory" added, VW
;   2008-04-09  V2: adapted algorithm away from [3], TB & VW
;   2008-04-14  Added sigma input, and TrC into Eigsys input, VW
;   2008-10-29  Updated to new algorithm, and tidied for public release, VW
;-
;****************************************************************************************;
; Copyright (C) 2007 Tamas Budavari and Vivienne Wild (MAGPop)                           ;
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

FUNCTION VWPCA_INC, eigsys, x, alpha, delta, conserve_memory=slow

; Return to caller on error.
On_Error, 2

; All inputs must be given
if N_ELEMENTS(eigsys) eq 0 or N_ELEMENTS(x) eq 0 or $
  N_ELEMENTS(alpha) eq 0 or N_ELEMENTS(delta) eq 0 then $
  MESSAGE, 'USAGE: neweigsys = VWPCA_INC(eigsys,x,alpha,delta)'

; Eigsys structure must contain correct tags
if NOT(TAG_EXIST(eigsys,'U')) or NOT(TAG_EXIST(eigsys,'W')) or $
  NOT(TAG_EXIST(eigsys,'m')) or NOT(TAG_EXIST(eigsys,'sig2')) or $
  NOT(TAG_EXIST(eigsys,'Vk')) or NOT(TAG_EXIST(eigsys,'Uk')) or NOT(TAG_EXIST(eigsys,'Qk')) then $
  MESSAGE, 'VWPCA_INC incorrect structure content: eigsys = {U,m,W,sig2,Vk,Uk,Qk}'

U = eigsys.U                    ;eigenvectors as (nbin,nvec) array
W = eigsys.W                    ;eigenvalues (nvec) vector
m = eigsys.m                    ;mean (nbin) vector
Vk_prev = eigsys.Vk
Uk_prev = eigsys.Uk
Qk_prev = eigsys.Qk
sig2 = eigsys.sig2

nbin = N_ELEMENTS(x)
nvec = N_ELEMENTS(W)


;;------
;; check for dimensionality mismatches

if (SIZE(U,/dim))[0] ne nbin or (SIZE(U,/dim))[1] ne nvec then $
  MESSAGE,'VWPCA_INC: eigsys.U must be an (nbin,nvec) matrix'
if (SIZE(W,/dim))[0] ne nvec then $
  MESSAGE,'VWPCA_INC: eigsys.W must be an (nvec) vector'
if (SIZE(m,/dim))[0] ne nbin then $
  MESSAGE,'VWPCA_INC: eigsys.m must be an (nbin) vector'


;;------
;; new observation vector, mean subtracted
y = x - m                       

;;------
;; residual error of vector y, using old eigenbasis
pcs = y ## TRANSPOSE(U)
yrecon = total(rebin(transpose(pcs),nbin,nvec)*U,2)
r = y-yrecon                    ;eq. 10 of [1]

;;------
;; the residual^2
sumr2 = total(r^2)

;;------
;; weight value, to determine weight to give to the new observation
;; vector, depending on it's residual from the original espace
hpi = !PI/2.0
weight1 = 1D / (1.0+ (hpi * sumr2/sig2)^2)
weight2 = (1.0/hpi) * atan(hpi * sumr2/sig2) / (sumr2/sig2)

vk = alpha*vk_prev + weight1          ;update total weight
qk = alpha*qk_prev + weight1*sumr2
uk = alpha*uk_prev + 1                ;update total number

gamma1 = (alpha*vk_prev) / vk   ;eq. 20 [1]
gamma2 = (alpha*qk_prev) / qk   ;eq. 21 [1]
gamma3 = (alpha*uk_prev) / uk   ;eq. 22 [1]

;;------
;; update estimate of sample variance eq. 19 [1]
sig2 = gamma3*sig2 + (1-gamma3)*weight2*sumr2 / delta


;;------
;; update mean  eq. 17 of [1]
Mout = m + (1-gamma1) * y 

;;------
;; derive new eigensystem eq. 18 of [1]

;; first, create new data matrix from old eigenvectors and new
;; observation vector 
if KEYWORD_SET(slow) then begin
    A = dblarr(nbin,nvec+1)
    for i=0L, nbin-1 do begin
        A[i,0:nvec-1] = SQRT(gamma2*W) * U[i,0:nvec-1] 
        A[i,nvec] = SQRT(1-gamma2) * y[i] * SQRT(sig2/sumr2) 
    endfor
endif else begin
;; time gain 0.45
;; vw: this is faster than matrix_multiply
    A = dblarr(nbin,nvec+1)
    A[*,0:nvec-1] = rebin(transpose(SQRT(gamma2*W)),nbin,nvec)*U 
    A[*,nvec] = SQRT(1-gamma2) * y * SQRT(sig2/sumr2)
endelse


;; now, instead of performing SVD on (nbin,nbin) covariance matrix, we
;; perform it on the smaller (p+1,p+1) matrix 
B = MATRIX_MULTIPLY(A,A,/ATRANSPOSE) ;eq. 9 of [3]

;; perform standard Single-Value decomposition (lapack routine)
;; note that the output matrices are back to front from standard convention
LA_SVD, B, Wout, Uout, Vout, /DOUBLE

;; and calculate new eigenvectors of covariance matrix
Uout = MATRIX_MULTIPLY(A,Uout,/BTRANSPOSE);eq 12 of [3]

;;------
;; truncate to desired number of eigenvectors
Uout = Uout[*,0:nvec-1]

;;------
;; here, eigenvalues are singular values  (not squared)
Wout = Wout[0:nvec-1]

;;------
;; (re)normalize eigenvectors
Uout /= (REPLICATE(1.0, nbin) # SQRT(Wout))


output = { M: Mout, W: Wout, U: Uout, sig2: sig2, Vk:Vk, Uk:Uk, Qk:Qk}

return, output

END
