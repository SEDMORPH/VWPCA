;+
;  NAME:
;    VWPCA_EM
;
; AUTHORS: 
;    Darren Madgwick
;    Vivienne Wild <wild@iap.fr>
;
;  PURPOSE:
;    Uses the Expectation-Maximization algorithm to solve for the
;    principal components of a given data-set.
;
;  CALLING SEQUENCE:
;    VWPCA_EM, data, k, [eigenvecs=eigenvecs, aver=aver, eigenvalues=eigenvalues,  
;              variances=variances, niter=niter,/silent]
;
;  INPUTS:
;    data   = Matrix of the data-vectors.  This should have 
;             dimensions (d x n), where d is the number of attributes
;             and n is the number of data-vectors. 
;    k      = The number of eigenvectors to be found.  If k << d
;             then this algorithm can be much faster than normal PCA.
;
;  OPTIONAL INPUTS: 
;    niter = number of iterations (default=15)
;
;  KEYWORD PARAMETERS:
;    silent  = print fewer messages
;
;  OPTIONAL OUTPUTS:
;    eigenvecs   = The named variable that will contain the
;                  eigenvectors once they have been calculated.
;    eigenvalues = Corresponding ESTIMATION OF eigenvalues.
;    aver        = Will contain the mean data-vector.  Note that because
;                  we do not construct the covariance matrix the
;                  mean vector must be explicitly subtracted.
;
; NOTES: 
;    (1) note that the variance estimates are approximate as assume all the vectors have been calculated. 
;
; MODIFICATION HISTORY:
;    2004 first implementation in IDL D. Madgwick (em_pcpts.pro)
;    2004-06-18 updated to return variances and include silent keyword V. Wild
;    2008-11-11 included pcpts routine as subroutine V. Wild
;-
;****************************************************************************************;
;  Copyright (c) 2004, Darren Madgwick and Vivienne Wild                                 ;
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

pro pcpts, data, eigenvecs=eigenvecs

  d = (size(data,/dim))[0]	; Dimensions of data-vectors
  n = (size(data,/dim))[1]	; Number of data-vectors


  eigenvecs = correlate(data, /covariance, /double)

  ; Compute the tridiagonal form
  trired, eigenvecs, diagonal, offdiagonal

  ; Compute the eigenvalues and eigenvectors
  triql, diagonal, offdiagonal, eigenvecs
  eigenvalues = diagonal

  diagonal=0
  offdiagonal=0

  eigenindex 	= reverse(sort(eigenvalues))
  eigenvalues 	= eigenvalues[eigenindex]
  eigenvecs 	= eigenvecs[*,eigenindex]

  cnt = where(abs(eigenvalues) lt 1e-6, c)
  if(c gt 0) then begin
    eigenvalues[cnt] = 0.0
    eigenvecs[*,cnt] = 0.0
  endif
  
end


pro VWPCA_EM, data, k, eigenvecs=eigenvecs, aver=aver, eigenvalues=eigenvalues,variances=variances,silent=silent,convergence=conv_out,niter=niter


  d = (size(data,/dim))[0]      ; Dimensions of data-vectors 
  n = (size(data,/dim))[1]      ; Number of data-vectors

;-- Calculate and subtract mean

; Even if this hass been done before, it is necessary for this
; routine, and won't harm to do again
; if you have small datasets this may be faster:
;  aver = total(data,2,/double)/double(n)
;  data = Temporary(data) - rebin(aver,d,n)

  aver = fltarr(d)
  for i=0, d-1 do aver[i]=mean(data[i,*])
  for i=0, n-1 do data[*,i]=data[*,i]-aver 


  X = fltarr(k,n)               ; Matrix of 'hidden' variables
  C = randomn(seed,d,k)         ; Matrix of mixtures

  if NOT KEYWORD_SET(niter) then niter = 15

  conv_out = fltarr(niter)
  for i=0, niter-1 do begin
      X     = float(invert((transpose(C) # C), /double) $ ; # data = memory gobbler, try float
                    # transpose(C)) # data
      C_new = data # float(transpose(X)#invert(X#transpose(X),/double))
      
      frac  = abs((C_new-C) / C)
      
      C     = C_new
      cnt   = where(finite(frac,/NaN) eq 1, ncnt)
      if(ncnt gt 0) then frac[cnt]=0.0
      conv_out[i] = total(frac)/n_elements(frac)
      if NOT KEYWORD_SET(SILENT) then print, 'Fractional convergence: ', i, conv_out[i]
  endfor
  X = 0


;-- Orthogonalize the columns of C
  for i=0, k-1 do begin
    for j=0, i-1 do begin
      dot    = total(C[*,i] * C[*,j],/double)
      C[*,i] = C[*,i] - dot*C[*,j]
    endfor 
        
    ; Normalize
    dot    = total(C[*,i]^2,/double)
    dot    = 1. / sqrt(dot)
    C[*,i] = C[*,i] * dot
  endfor 



;-- Normal PCA now solves for the eigenvectors (bit of a memory eater)
  pcpts, transpose(C)#data, eigenvecs=eigenvecs 
  eigenvecs = C # eigenvecs

  if NOT KEYWORD_SET(silent) then print, 'Estimating eigenvalues...'

; Calculate two rows of the covariance matrix
;;each row should give same result for the eigenvalue
;; (A evec1 = lambda_1 evec) i.e. for evec1 get lambda_1 from
;; covariance matrix of data

;-- use middle of data array, in case of redshifted spectra
    tmp = d/2
    row1 = fltarr(d)
    row2 = row1
    for i=0, d-1 do begin
        row1[i] = correlate(data[tmp,*],data[i,*],/covariance,/double)
        row2[i] = correlate(data[tmp+1,*],data[i,*],/covariance,/double)
    endfor

;-- for each eigenvector get corresponding eigenvalue
    evalues1 = fltarr(k)
    evalues2 = fltarr(k)
    for i=0, k-1 do begin       
        evalues1[i] = (transpose(row1) # eigenvecs[*,i])[0] / eigenvecs[tmp,i] ;these 2 should be the same: lambda_i
        evalues2[i] = (transpose(row2) # eigenvecs[*,i])[0] / eigenvecs[tmp+1,i] 
    endfor

;-- because evecs are estimated, average 2 evals to give better approx. Lower down the evalues are less accurate.
    eigenvalues = 0.5*(evalues1 + evalues2) 

;-- calculate variances
;variances = eigenvalues/TRACE(cov)
;calculate trace of cov matrix
;tr = fltarr(d)
;for i=0,d-1 do tr[i] = total(data[i,*]*data[i,*])/double(n)

tr = total(data*data,2)/double(n)    
variances = eigenvalues/total(tr)
       
end
