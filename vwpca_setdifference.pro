;;*** SETDIFFERENCE.PRO ***

;   a = [2,4,6,8]
;   b = [6,1,3,2]

; SetIntersection(a,b) = [ 2, 6]       ; Common elements
; SetUnion(a,b) = [ 1, 2, 3, 4, 6, 8]  ; Elements in either set
; SetDifference(a,b) = [ 4, 8]         ; Elements in A but not in B
; SetIntersection(a,[3,5,7]) = -1      ; Null Set

; From COYOTE IDL library http://www.dfanning.com/tips/set_operations.html
; 24 Nov 2008 Changed function name to avoid any incompatibility problems when
; running VWPCA package (V. Wild)


FUNCTION VWPCA_SetDifference, a, b  

   ; = a and (not b) = elements in A but not in B

mina = Min(a, Max=maxa)
minb = Min(b, Max=maxb)
IF (minb GT maxa) OR (maxb LT mina) THEN RETURN, a ;No intersection...

if n_elements(b) eq 1 then begin
    ind = where(a eq b[0],compl=compl)
    return,a[compl]
endif

r = Where((Histogram(a, Min=mina, Max=maxa) NE 0) AND $
          (Histogram(b, Min=mina, Max=maxa) EQ 0), count)
IF count eq 0 THEN RETURN, -1 ELSE RETURN, r + mina
END
