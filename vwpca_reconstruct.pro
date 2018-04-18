;+
; NAME: VWPCA_RECONSTRUCT
;
; AUTHORS: 
;    Vivienne Wild <wild@iap.fr>
;
; PURPOSE: 
;   To reconstruct data-vector using the principal component
;   amplitudes of the data vector and eigenbasis
;    
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;   reconstructed_data = VWPCA_RECONSTRUCT(pcs,evects)
;
; INPUTS:
;   pcs      - principal component amplitudes (n-amplitude)
;   evects   - eigenvectors (n-bins x n-evectors matrix)
;
; OUTPUTS:
;   reconstructed data vector
;
; MODIFICATION HISTORY:
;   2006 First implemented in IDL V. Wild
;   2017 Modified header to fix incorrect information on multiple data
;   vectors. It doesn't work for multiple data vectors. 
;-
;****************************************************************************************;
;  Copyright (c) 2006, Vivienne Wild 
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


FUNCTION VWPCA_RECONSTRUCT, pcs, espec

nbin = (size(espec,/dim))[0]
if (size(espec))[0] eq 1 then nrecon=1 else nrecon = (size(espec,/dim))[1]

if nrecon ne 1 then begin
    if n_elements(pcs) ne nrecon then begin
        print,'PCA_RECONSTRUCT: pcs or espec array wrong size'
        return,-1
    endif
    result = total(rebin(Transpose(pcs),nbin,nrecon)*espec,2)
endif else begin
    if n_elements(pcs) ne 1 then begin
        print,'PCA_RECONSTRUCT: pcs or espec array wrong size'
        return,-1
    endif
    result = pcs[0]*espec
endelse



return, result

END
