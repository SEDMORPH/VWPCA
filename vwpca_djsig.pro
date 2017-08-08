; wrapper for djs_iterstat
; D. Finkbeiner 14 Oct 1999
; DPF 27 Jun 2000 - added keywords
;
; 24 Nov 2008 Changed function name to avoid any incompatibility problems when
; running VWPCA package (V. Wild)

function vwpca_djsig, x, sigrej=sigrej, maxiter=maxiter

  vwpca_djs_iterstat, x, sigma=sigma, sigrej=sigrej, maxiter=maxiter

  return,sigma
end
