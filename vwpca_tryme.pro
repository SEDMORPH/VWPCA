;+
; NAME: 
;    VWPCA_TRYME
;
; AUTHORS: 
;    Vivienne Wild <wild@iap.fr>
;
; PURPOSE:
;    To exemplify the running of PCA using several different
;    decomposition methods on a test dataset of composite stellar
;    population (galaxy) spectra 
;   
; CATEGORY:
;   Statistics
;
; CALLING SEQUENCE:
;   IDL> vwpca_tryme
;
; ADDITIONAL PROGRAMS REQUIRED: 
;   IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;
;   The following are required and included in the distribution:
;   - vwpca_djsig.pro and vwpca_djs_iterstat.pro from SDSS IDL-UTILS
;   - vwpca_setdifference.pro from COYOTE IDL library
;
;
; NOTES:
;   (1) The test data provided in this package contains 
;       ERRARR          FLOAT     = Array[nbin, ngal]
;       NORM            FLOAT     = Array[ngal]
;       SPECARR         FLOAT     = Array[nbin, ngal]
;       WAVE            DOUBLE    = Array[nbin]
 
;   (2) The test data provided here has been normalised i.e. the
;   luminosity information has been removed from the galaxy
;   spectra. It is recommended to find a good way to normalise your
;   own data, e.g the mean/median flux of all/part-of the data
;   vector. How you normalise will affect the results. 
;   
;
; MODIFICATION HISTORY:
;  2008 First implemented in IDL V.Wild
;  2009-01-26 Removed REMOVE_OUTLIERS sub-function to create stand-alone
;             routine. 
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
PRO VWPCA_TRYME

;;-- Collect the data

str = ''
read,str,prompt= 'Would you like to try model (m/M) or real (r/R) spectra?'
if str eq 'm' or str eq 'M' then restore, 'vwpca_modelspec.sav' $
else if str eq 'r' or str eq 'R' then restore, 'vwpca_realspec.sav' $
else begin
    print, 'sorry your response was not understood!'
    return
endelse

nbin = (size(specarr,/dim))[0]
ngal = (size(specarr,/dim))[1]

;;-- number of PC amplitudes to calculate
namplitudes = 4

;;-- fill any gaps (only relavant for the real spectra)

ind = where(errarr eq 0.,count)
if count ne 0 then begin
    for i=0L,nbin-1 do begin
        ind = where(errarr[i,*] eq 0.,compl=compl)
        specarr[i,ind] = median(specarr[i,compl])
    endfor
endif


set_plot,'ps'
device,file='vwpca_tryme.ps',/color,bits_per_pixel=8
device,/portrait,xoffset=2,yoffset=2,ysize=25,xsize=18
!p.charsize=(!x.charsize=(!y.charsize=1.2))
!p.multi=[0,1,5]
loadct,13,nc=30

;;------------------------------------------------------------------
;;***  Now exemplify each variety of PCA using this dataset
;;------------------------------------------------------------------


print, '------------------------------------------------------------------'
print, '*** SINGULAR VALUE DECOMPOSITION ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,/svd)

print, 'Top 10 variances: ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time


plot, wave,meanarr,title='SVD: mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor

for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='SVD: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor

print, '------------------------------------------------------------------'
print, '*** SINGULAR VALUE DECOMPOSITION using error array to calculate PC amplitudes ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,errarr=errarr,norm=norm,/svd)

print, 'Top 10 variances: ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time



plot, wave,meanarr,title='SVD with errors: mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor

for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='SVD with errors: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,(vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr)*norm[i],color=28
endfor


print, '------------------------------------------------------------------'
print, '*** IDL-ASTRO routine ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,/astpca)

print, 'Top 10 variances: ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time



plot, wave,meanarr,title='IDL-AST PCA:  mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor
 
for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='IDL-AST PCA: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor


print, '------------------------------------------------------------------'
print, '*** RSI PCOMP ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,/pcomp)

print, 'Top 10 variances: ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time



plot, wave,meanarr,title='RSI PCOMP: mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor
for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='RSI PCOMP: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor


print, '------------------------------------------------------------------'
print, '*** EXPECTATION MAXIMISATION ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,/empca)

print, 'Top 10 variances (approx): ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time



plot, wave,meanarr,title='EM-PCA: mean + top 4 eigenspectra',/xs
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys
endfor
for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='EM-PCA: e.g. spectrum (black)+PCA reconstruction (red)'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor


print, '------------------------------------------------------------------'
print, '*** SVD trimmed (several iterations to remove outliers) ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

niter = 10
time = systime(1)

specarr_old = specarr
ngood_old = ngal
for j=0,niter-1 do begin
   espec = vwpca(specarr_old,namplitudes,variances,pcs,meanarr,evalues,/svd)
   specarr_old = vwpca_outliers(specarr_old,pcs,4.5,good=good)

   ngood = n_elements(good)
   print,'nb good spectra:',ngood
   if ngood eq ngood_old then break ;thrown none out
   ngood_old = ngood

endfor

print, 'Top 10 variances: ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time

plot, wave,meanarr,title='SVD Trimmed: mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor
for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='SVD trimmed: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor

print, '------------------------------------------------------------------'
print, '*** ROBUST ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,/robust)

print, 'Top 10 variances (approx): ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time



plot, wave,meanarr,title='ROBUST (no error array): mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev=1 else rev=-1
    plot, wave,rev*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor
for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='ROBUST: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor

print, '------------------------------------------------------------------'
print, '*** ROBUST w/ errors ***'
variances = 0 & evalues = 0 & pcs = 0 & espec = 0 & meanarr = 0

params = {niter:5, nevec:10}

;;-- add random gaps to model spectra
if str eq 'm' or str eq 'M' then begin
    rand = randomu(234075609,nbin*ngal)
    ind = where(rand gt 0.9)
    errarr[ind] = 0.0
endif

time = systime(1)
espec = vwpca(specarr,namplitudes,variances,pcs,meanarr,evalues,errarr=errarr,/robust,params=params)

print, 'Top 10 variances (approx): ',variances[0:9],form='(A18,10(F0.4,1X))'
print, 'Top 10 evalues  : ',evalues[0:9],form='(A18, 10(F0.4,1X))'
print, 'Time: ',systime(1)-time



plot, wave,meanarr,title='ROBUST, w/ errors+gaps: mean + top 4 eigenspectra',/xs,xtitle='Wavelength',ytitle='Flux'
rev = lonarr(4)
for i=0,3 do begin
    if max(espec[*,i]) gt abs(min(espec[*,i])) then rev[i]=1 else rev[i]=-1
    plot, wave,rev[i]*espec[*,i],/xs,yr=[-0.1,0.3],/ys,xtitle='Wavelength',ytitle='Flux'
endfor
for i=0,4 do begin
    plot, wave,specarr[*,i],psym=10,/xs,title='ROBUST: e.g. spectrum (black)+PCA reconstruction (red)',xtitle='Wavelength',ytitle='Flux'
    oplot, wave,vwpca_reconstruct(pcs[*,i],espec[*,0:namplitudes-1])+meanarr,color=28
endfor

;;******************************************************************

;;-- close ps file
cleanplot,/silent
device,/portrait
device,/close
set_plot,'x'

END
