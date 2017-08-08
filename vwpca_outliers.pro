;+
; NAME: 
;    VWPCA_OUTLIERS
;
; AUTHOR:
;    Vivienne Wild <wild@iap.fr>
;
; PURPOSE: 
;    given a spectral dataset and set of associated parameters
;    remove the outliers from the dataset and return the new data
;    array 
; 
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;    newdata = VWPCA_OUTLIERS(data,pcs,sigma )
; 
; INPUTS:
;    data = The 2-d data array (NBIN x NOBJ) i.e. NOBJ is the total
;           number of objects and NBIN is the total number of attributes
;    pcs  = The 2-d data array (NPARAM x NOBJ) of associated
;           parameters
;    sigma = The required sigma-cut above which outliers are identified
;
; OUTPUTS:
;    newdata   = NBIN x NOBJ_new matrix containing the reduced data array                
;
; OPTIONAL OUTPUTS:
;    good    = indices of data vectors which have been retained
;
; KEYWORDS:
;    silent  = don't tell user anything
; 
; EXAMPLE:
;    data_new = vwpca_outliers(data,pcs,4.5,good=good)
;
;     
; MODIFICATION HISTORY:
;
;    2006 First implementation (V. Wild)
;    2009-01-26 Prepared for public release (V. Wild)
;
;-
;****************************************************************************************;
;  Copyright (c) 2009, Vivienne Wild                                                     ;
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


FUNCTION VWPCA_OUTLIERS,data,pcs,sigma,good=good,silent=silent

if (SIZE(pcs))[0] eq 1 then nrecon=1 else nrecon = (SIZE(pcs,/dim))[0]
ngal = (SIZE(data,/dim))[1]

FOR i=0,nrecon-1 DO BEGIN
    ind = where(ABS(pcs[i,*]) GT sigma*vwpca_djsig(pcs[i,*],sigrej=3.0))
    IF ind[0] ne -1 THEN IF N_ELEMENTS(ind_outlier) EQ 0 THEN ind_outlier=ind ELSE ind_outlier = [ind_outlier,ind]
ENDFOR

IF N_ELEMENTS(ind_outlier) NE 0 THEN ind_outlier = ind_outlier[UNIQ(ind_outlier,SORT(ind_outlier))]

if not(keyword_set(silent)) then PRINT, ' rejected ', N_ELEMENTS(ind_outlier),' objects'
IF N_ELEMENTS(ind_outlier) EQ 0 THEN BEGIN
    good = FINDGEN(ngal)
    RETURN,data
ENDIF

good = VWPCA_SETDIFFERENCE(FINDGEN(ngal),ind_outlier)
newdata = data[*,good]

RETURN,newdata


end
