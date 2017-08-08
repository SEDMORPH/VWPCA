;+
; NAME: 
;    VWPCA_SVD
;
; AUTHOR:
;    Vivienne Wild <wild@iap.fr>
;
; PURPOSE: 
;    Carry out a Principal Component Analysis using singular-value-decomposition
; 
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;    evectors = VWPCA_SVD(data,evalues=evalues,pcs=pcs,variances=variances)
; 
; INPUTS:
;    data = The 2-d data array (NBIN x NOBJ) i.e. NOBJ is the total
;           number of objects and NBIN is the total number of attributes
;
; OUTPUTS:
;    evectors  = NBIN x N-EVECTOR matrix containing the
;                eigenvectors. Note that N-EVECTOR = min(NOBJ,NBIN). 
;
; OPTIONAL OUTPUTS:
;    evalues   = NBIN vector containing the corresponding eigenvalues
;    pcs       = N-EVECTOR x NOBJ matrix containing the principal
;                component amplitudes of the input dataset   
;    variances = NBIN vector containing an ESTIMATE of the variance
;                carried by each eigenvector
;   
; EXAMPLE:
;    evectors = vwpca_svd(data)
;
; NOTES: 
;    (1) If only a few principal component amplitudes are required, it is
;        faster to compute these afterwards using:  
;        pcs = data ## TRANSPOSE(evectors[*,0:namplitudes-1]) 
;     
; MODIFICATION HISTORY:
;
;    2006 First implimentation (V. Wild)
;    2008-11-11 Prepared for public release (V. Wild)
;
;-
;****************************************************************************************;
;  Copyright (c) 2008, Vivienne Wild                                                     ;
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

FUNCTION VWPCA_SVD, data,variances=variances,pcs=pcs,evalues=evalues

;;-- normalise the data vector so the SVD is equivalent to e-value
;;   decomposition on the Covariance matrix

nobj = (size(data,/dim))[1]

;;-- perform decomposition

LA_SVD,data/float(sqrt(nobj)), W, U, V
evectors = TRANSPOSE(V)

;-- calculate variances, evalues and pcs if required
;   note that the variances are only an estimate: it assumes we have calculated _all_ eigenvectors

if ARG_PRESENT(variances) then begin
    variances = W^2
    variances = variances/total(variances) 
endif

if ARG_PRESENT(evalues) then evalues = W^2
if ARG_PRESENT(pcs) then pcs = data ## V


return, evectors
END
