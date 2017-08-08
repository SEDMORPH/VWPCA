;+
; NAME: VWPCA
;
; AUTHORS: 
;    Vivienne Wild <wild@iap.fr>
;
; PURPOSE: Wrapper for performing a principal component analysis on a data matrix. 
;          Several options exist for the algorithm used to compute the PCA, with all 
;          differences between the algorithms taken account of in this routine. 
;          Principal Component amplitudes may also be calculated if required.
;    
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;    evectors = VWPCA(data,namplitudes,variances,pcs,meanarr,evalues,errarr=,vararr=,norm=norm,
;                     params=,evals_rob=,/EMPCA,/SVD,/PCOMP,/ASTPCA,/ROBUST,/VARWEIGHT,/NOMEANSUB)
;
; INPUTS:
;    data        = [nbin,nobj] matrix on which to perform PCA
;    namplitudes = nb of PCA amplitudes to calculate (=0 if none)
;
; OPTIONAL INPUTS:
;    params      = structure containing parameters for use with EM-PCA or ROBUST PCA
;    errarr      = [nbin,nobj] matrix of errors (either 1-sigma or 0=bad pixel, 1=good pixel)
;
; KEYWORD PARAMETERS:
;    /ASTPCA    = PCA from IDL astro 
;    /PCOMP     = PCA from inbuilt RSI routine 
;    /EMPCA     = PCA using expectation maximisation                        
;    /SVD       = PCA using singular value decomposition
;    /ROBUST    = PCA using robust and iterative routine of Budavari & Wild (2008) 
;    /NOMEANSUB = do not subtract mean (à la Connolly & Szalay). This
;                 may not be compatible with all subroutines.
;    /VARWEIGHT = divide by standard deviation? (currently not
;                 applicable for ROBUST. See also note (2))
;
; OUTPUTS:
;    evectors    = [nbin,nevec] matrix of eigenvectors (returned)
;    variances   = variances accounted for in each eigenvector
;    pcs         = PC amplitudes for each input data vector (if namplitudes !=0)
;    meanarr     = the calculated mean array 
;    evalues     = eigenvalues corresponding to eigenvectors
;
; OPTIONAL OUTPUTS:
;    vararr      = the calculated variance array (if varweight=1)
;    evals_rob   = the eigenvalues at each stage in the robust-PCA
;                  iteration (if ROBUST=1)
;    norm        = the normalisation of the data vectors during
;                  projection (if namplitudes>0, the errarr is input
;                  and varweight=0) 
;
; ADDITIONAL PROGRAMS REQUIRED: 
;   IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;
; REFERENCES: [1] Budavari, Wild et al. (2008), MNRAS
;             [2] Connolly & Szalay (1999, AJ, 117, 2052)
; 
; NOTES: 
;   1) The params structure contains several parameters for personalisation of the ROBUST or EMPCA routines
;           params.niter      = nb of iterations (default=5 for ROBUST; default=20 for EMPCA)
;           params.nevec      = nb of eigenvectors to return/calculate (default=10 both ROBUST and EMPCA)
;           params.nobj_init  = nb of objects used to initialise (ROBUST only; default 10% of dataset)
;           params.delta      = parameter to adjust robustness (read [1])
;           params.memory     = parameter to adjust memory (read [1])
;   2) If the corresponding ERROR array is input then the PC
;      amplitudes are computed using the Normalised-Gappy PCA algorithm
;      of G. Lemson, adapted from Connolly & Szalay. Otherwise a straight matrix 
;      multiplication is performed. Note that in the former case, the
;      PC amplitudes are dependent on the number of eigenvectors used
;      in the calculation (see [2]) The use of the VARWEIGHT keyword
;      is currently incompatible with the use of the Normalised-gappy
;      PCA algorithm. Therefore, if /VARWEIGHT is set, the PC
;      amplitudes are calculated by matrix multiplication.  
;   3) In the case of any truncated method (EMPCA, ROBUST), the variances are
;      approximate only, because they assume that all the eigenvectors have
;      been calculated.
;
; MODIFICATION HISTORY:
;   2008 First implementation in IDL V. Wild
;   2008-12-07 added error check for namp>nobj in SVD
;-
;****************************************************************************************;
;  Copyright (c) 2008, Vivienne Wild 
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

FUNCTION VWPCA, dataarr, namplitudes, variances, pcs, meanarr, evalues, errarr=errarr, params=params, $
                EMPCA=EMPCA, PCOMP=RSI, ASTPCA=ASTPCA, SVD=SVD, ROBUST=ROBUST, $
                vararr=vararr, evals_rob=evalues_it, norm=norm, $
                nomeansub=nomeansub, varweight=varweight

; Return to caller on error.
On_Error, 2

;;******************************************************************
;;*** check inputs
;;******************************************************************

if N_ELEMENTS(SIZE(dataarr,/dim)) ne 2 then begin
    print, 'VWPCA ERROR: input array must be 2D'
    return,-1
endif 

tmp = SIZE(dataarr,/dim)
nbin = double(tmp[0])
nobj = double(tmp[1])

if namplitudes gt nbin then begin
    print, 'VWPCA WARNING: namplitudes must be < number of data points, setting namplitudes=nbin'
    namplitudes=nbin
endif

if keyword_set(SVD) and namplitudes gt nobj then begin
   print, 'VWPCA WARNING: for SVD namplitudes must be < number of objects, setting namplitudes=nobjects'
    namplitudes=nobj
endif

print, 'VWPCA: ',long(nobj),'-object x ',long(nbin),'-data array ',form='(A,I0.0,A,I0.0,A)'
print, 'VWPCA: will return ',namplitudes,' PCA amplitudes',form='(A,I0.0,A)'

if NOT(KEYWORD_SET(EMPCA)) and  NOT(KEYWORD_SET(RSI)) and NOT(KEYWORD_SET(ASTPCA)) $
  and NOT(KEYWORD_SET(SVD)) and NOT(KEYWORD_SET(ROBUST))then begin
    print, 'VWPCA ERROR: please specify which PCA code (/EMPCA, /PCOMP, /ASTPCA, /SVD, /ROBUST)'
    return,-1
endif

if KEYWORD_SET(EMPCA) then EMPCA=1 else EMPCA=0
if KEYWORD_SET(RSI) then RSI=1 else RSI =0
if KEYWORD_SET(ASTPCA) then ASTPCA=1 else ASTPCA=0
if KEYWORD_SET(SVD) then SVD=1 else SVD=0
if KEYWORD_SET(ROBUST) then ROBUST=1 else ROBUST=0

if KEYWORD_SET(nomeansub) then meansub=0 else meansub=1
if KEYWORD_SET(varweight) then varweight=1 else varweight=0

if (ROBUST OR EMPCA OR SVD) AND meansub eq 0 then begin
;- the input matrix MUST be centered (the mean subtracted from each pixel)
    print, 'VWPCA MESSAGE: SVD/ROBUST/EMPCA require a centered matrix, enforcing meansub=1'
    meansub = 1
endif


;;------------------------------------------------------------------
;;*** special parameters for EMPCA
;;------------------------------------------------------------------

if EMPCA then begin
    nevec = 10                  ;default parameters
    niter = 20

    if N_ELEMENTS(params) ne 0 then begin
        if TAG_EXIST(params,'nevec') eq 0 then print, 'VWPCA MESSAGE: computing default nb of eigenvectors ',nevec $
        else nevec = params.nevec
        if TAG_EXIST(params,'niter') eq 0 then print, 'VWPCA MESSAGE: computing default nb of iterations through dataset ',niter $
        else niter = params.niter
     endif else begin
        print, 'VWPCA MESSAGE: using ALL default parameters  '
        print, '                number of iterations: ', niter
        print, '                number of eigenspectra: ', nevec
     endelse 

    if namplitudes gt nevec then begin
        print, 'VWPCA ERROR: cannot have n-amplitudes>n-evectors, returning'
        return,-1
    endif
endif

;;------------------------------------------------------------------
;;*** special parameters for ROBUST PCA
;;------------------------------------------------------------------

if ROBUST then begin

    if varweight then begin
        print, 'VWPCA MESSAGE: for ROBUST, cannot weight data by variance'
        return,-1
    endif

;;--- default parameters
    nevec  = 10
    niter  = 5
; Number of objects to use in initialisation of PCA
    nobj_init = MIN([nobj / 100.,200.])
; Control parameter for robustness (smaller=remove more outliers)
    delta  = 0.5D
; Control parameter for convergence (larger=remember past; 1 = remember one iteration of dataset)
    memory = 1D


    if N_ELEMENTS(params) ne 0 then begin
        if TAG_EXIST(params,'nevec') eq 0 then print, 'VWPCA MESSAGE: computing default nb of eigenvectors ',nevec $
        else nevec = params.nevec
        if TAG_EXIST(params,'niter') eq 0 then print, 'VWPCA MESSAGE: computing default nb of iterations through dataset ',niter $
        else niter = params.niter
        if TAG_EXIST(params,'nobj_init') eq 0 then print, 'VWPCA MESSAGE: using default nb of initial objects ',nobj_init $
        else nobj_init = params.nobj_init
        if TAG_EXIST(params,'delta') eq 0 then print, 'VWPCA MESSAGE: using default delta (robustness parameter) ',delta $
        else delta = params.delta
        if TAG_EXIST(params,'memory') eq 0 then print, 'VWPCA MESSAGE: using default memory (convergence parameter) ',memory $
        else memory = params.memory
     endif else begin
        print, 'VWPCA MESSAGE: using ALL default parameters  '
        print, '                convergence: ', memory
        print, '                robustness: ', delta
        print, '                number of iterations: ', niter
        print, '                number of initial objects: ', nobj_init
        print, '                number of eigenspectra: ', nevec
     endelse 

    if namplitudes gt nevec then begin
        print, 'VWPCA ERROR: cannot have n-amplitudes>n-evectors, returning'
        return,-1
    endif

endif

;;******************************************************************
;;*** data with gaps
;;******************************************************************

if N_ELEMENTS(errarr) ne 0 then GAPPY=1 else GAPPY=0

if GAPPY then begin
    if N_ELEMENTS(SIZE(errarr,/dim)) ne 2 then begin
        print, 'VWPCA MESSAGE: error array must be 2D'
        return,-1
    endif 
    
    tmp = SIZE(errarr,/dim)
    if tmp[0] ne nbin OR tmp[1] ne nobj then begin
        print, 'VWPCA ERROR: error array must be same size as data array'
        return,-1
    endif
endif

;;******************************************************************
;;*** prepare array
;;******************************************************************

;-- create a centralised array (mean subtracted, and variance weighted if requested)
dataarr_c = dataarr              

;;-- subtract mean
if meansub then begin

    meanarr = TOTAL(dataarr,2,/nan) / double(nobj) 
 
    for i=0L,nobj-1 do dataarr_c[*,i] = dataarr[*,i]-meanarr
    if (size(dataarr,/type)) eq 5 then meanarr = double(meanarr) else meanarr=float(meanarr) ;retain same type in output
;-- this is faster but uses more memory
;dataarr_c = TEMPORARY(dataarr_c) - rebin(meanarr,nbin,nobj)

endif else $
  if (size(dataarr,/type)) eq 5 then meanarr = dblarr(nbin) else meanarr = fltarr(nbin)


;;-- divide by variance
if varweight then begin
    
    vararr = SQRT(TOTAL(dataarr_c^2,2,/nan)/(double(nobj)-1))

    for i=0L,nobj-1 do dataarr_c[*,i] = dataarr_c[*,i]/vararr
    if (size(dataarr,/type)) eq 5 then vararr = double(vararr) else vararr=float(vararr) ;retain same type in output
;-- this is faster but uses more memory
;dataarr_c = TEMPORARY(dataarr_c)/rebin(vararr,nbin,nobj)

endif else $
  if (size(dataarr,/type)) eq 5 then vararr = dblarr(nbin)+1D else vararr = fltarr(nbin)+1.


;;******************************************************************
;;*** Perform PCA
;;******************************************************************

;;------------------------------------------------------------------
;;*** SINGULAR VALUE DECOMPOSITION
;;------------------------------------------------------------------

if SVD then begin
    evect= VWPCA_SVD(dataarr_c,variances=variances,evalues=evalues)

;-- we compute namplitude pcs outside of function to save time if
;   namplitudes < n-evectors
    if namplitudes ne 0 then begin
        if gappy and not(varweight) then $
          pcs = vwpca_normgappy(dataarr,errarr, evect[*,0:namplitudes-1], meanarr,norm=norm) $
        else pcs = dataarr_c ## TRANSPOSE(evect[*,0:namplitudes-1]) 
    endif else pcs = 0
endif

;;------------------------------------------------------------------
;;*** RSI - PCOMP (eigenvalue decomposition EIGENQL.PRO)
;;------------------------------------------------------------------


if RSI then begin
;-- In PCOMP All PC amplitudes are calculated regardless of requested
;   number. If namplitudes=0, then ALL pcs are returned. 
    pcs = PCOMP(dataarr_c, coeff=evect, nvariables=namplitudes, eigenvalues=evalues, $
                variances=variances, /cov, /double)

;-- the RSI PCA normalises the eigenvectors by the square of the
;   eigenvalues. Here we undo this for consistency with other
;   routines.
    evect = evect / (REPLICATE(1.0, nbin) # SQRT(evalues))

;-- re-calculate PC amplitudes if we have error information
    if namplitudes ne 0 then begin
        if gappy and not(varweight) then $
          pcs = vwpca_normgappy(dataarr,errarr, evect[*,0:namplitudes-1], meanarr,norm=norm) $
        else pcs = pcs/ REBIN(SQRT(evalues[0:namplitudes-1]), namplitudes,nobj)
    endif else pcs = 0

endif

;;------------------------------------------------------------------
;;*** Expectation Maximisation 
;;------------------------------------------------------------------

if EMPCA then begin
    VWPCA_EM, dataarr_c, nevec, eigenvecs=evect,variances=variances,$
              niter=niter,convergence=conv_out,eigenvalues=evalues

    if namplitudes ne 0 then begin
        if gappy and not(varweight) then $
          pcs = vwpca_normgappy(dataarr,errarr, evect[*,0:namplitudes-1], meanarr,norm=norm) $
        else pcs = dataarr_c ## TRANSPOSE(evect[*,0:namplitudes-1]) 
    endif else pcs = 0
endif

;;------------------------------------------------------------------
;;*** IDL-ASTRO PCA
;;------------------------------------------------------------------

if ASTPCA then begin

;- variances are cumulative and percentage, eigenvalues are scaled by
;  Nobj, pcs are not transposed like their input matrix
    VWPCA_IDLAST, TRANSPOSE(dataarr_c), evalues, evect, cuml_percent, pcs,/cov,/silent

;- re-scale eigenvalues
    evalues = evalues / double(nobj)

;- put evectors back the right way up
    evect = TRANSPOSE(evect)
    
    if namplitudes ne 0 then begin
        if gappy and not(varweight) then $
          pcs = vwpca_normgappy(dataarr,errarr, evect[*,0:namplitudes-1], meanarr,norm=norm) $
        else pcs = pcs[0:namplitudes-1,*]
    endif else pcs = 0   

;- turn cumulative percentages back into variances    
    variances = cuml_percent/100.
    for i = n_elements(variances)-1L,1L,-1L do $
      variances[i] = variances[i]-variances[i-1]
endif
 

;;------------------------------------------------------------------
;;*** ROBUST + ITERATIVE
;;------------------------------------------------------------------

if ROBUST then begin

;- initial starting parameters
    alpha = 1.-1/double(nobj*memory) ; set convergence parameter depending on number of datapoints

;- compute initial eigenbasis using small number of objects
    mean_init = meanarr
    evect_init = VWPCA_SVD(dataarr_c,pcs=pcs,evalues=evals_init)
    evect_init = evect_init[*,0:nevec-1]
    evals_init = evals_init[0:nevec-1]
    pcs = pcs[0:nevec-1,*]

;- compute Euclidean distance between eigenbasis and each initialisation
;  spectrum. Note that this ignores possible gaps in spectra and thus
;  underestimates the true r^2. The correct solution to this is complicated.

    sumr2 = fltarr(nobj_init)
    for i=0L,nobj_init-1 do begin
        if GAPPY then ind = where(errarr[*,i] ne 0) else ind = findgen(nbin)
        sumr2[i] = total((dataarr_c[ind,i] - (vwpca_reconstruct(pcs[*,i],evect_init))[ind])^2)
    endfor

;- estimate starting point for the scale, sig^2
    sig2 = mean(sumr2)          
    hpi = !PI/2.0
    for i=0,nevec-1 do sig2 = sig2*mean((1.0/hpi)*atan(hpi*sumr2/sig2)) / delta

;- calculate starting weights 
    vk = mean(1D / (1.0+ (hpi * sumr2/sig2)^2)) ;Eqn. 20 [1]
    qk = mean(sumr2 / (1.0+ (hpi * sumr2/sig2)^2)) ;Eqn. 21 [1]
    uk = 1.                     ;Eqn. 22 [1]

;- input for robust PCA
    input = { M: mean_init,  W: evals_init[0:nevec-1], U: evect_init[*,0:nevec-1], sig2:sig2, Uk:Uk, Vk:Vk, Qk:Qk }

;- Now perform iterative PCA: store all these parameters for debugging
    evalues_it = fltarr(nevec,nobj*niter-nobj_init+1)
    Vk_it = fltarr(nobj*niter-nobj_init+1)
    sig2_it = fltarr(nobj*niter-nobj_init+1)
    nn = 0L

    for k=0,niter-1 do begin
        print,'VWPCA: iteration number ',k
        if k eq 0 then begin
            evalues_it[*,nn] = evals_init[0:nevec-1]
            Vk_it[nn] = Vk
            sig2_it[nn] = sig2
            nn= nn+1
            start = nobj_init 
        endif else start = 0
        
        for i=start,nobj-1 do begin
            if GAPPY then output = vwpca_incgappy(input,double(dataarr[*,i]),double(errarr[*,i]),$
                                                alpha,delta,/normgappy) $
            else output = vwpca_inc(input,double(dataarr[*,i]),alpha,delta)
            
            ;; store some useful parameters
            evalues_it[*,nn] = output.w
            Vk_it[nn] = output.Vk
            sig2_it[nn] = output.sig2
            nn+=1
            
            input = output
        endfor
    endfor
    
    
;- final output
    evect = output.U
    evalues = output.W
    meanarr = output.M
    variances = output.W/total(output.W) ;approximate!!!!
    
;- calculate pcs if asked for       
    if namplitudes ne 0 then begin
        if GAPPY then begin
            pcs = vwpca_normgappy(dataarr,errarr, evect[*,0:namplitudes-1], meanarr,norm=norm)
        endif else begin
            ;; Because the normalisation of the components of the
            ;; eigensystem is not fixed, it's a good idea to
            ;; renormalise the individual data to the level of the
            ;; meanarr of the eigensystem
            dataarr_c = (dataarr /replicate(mean(meanarr),nbin,nobj)) - rebin(meanarr,nbin,nobj)
            pcs = TRANSPOSE(evect[*,0:namplitudes-1]) # dataarr_c
        endelse
    endif
    
endif

;- that's all folks!
return, evect


END
