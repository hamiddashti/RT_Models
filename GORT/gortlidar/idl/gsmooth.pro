function gsmooth,wave,width,SHIFT=shift

; J. Bryan Blair
; NASA/GSFC Code 924
; Greenbelt, MD 20771
; James.B.Blair@nasa.gov


N=n_elements(wave)
N2 = fix(width*8)
if(keyword_set(SHIFT)) then N2 = fix(width*8+2*abs(shift))
if(N2 gt N/2) then begin
	print,'EXCEED MAX KERNAL SIZE!',N2,N
	N2 = N/2
endif
N2 = fix(N2/2)*2 + 1
shiftval = 0
if(keyword_set(SHIFT)) then shiftval = shift
g1 = gaussf(findgen(N2),1.0,float(fix(N2/2))+shiftval,float(width))
g1 = g1/float(total(g1))
g1 = reverse(g1)
;return,convol(float(wave),g1,/edge_wrap,/center)
return,convol(float(wave),g1,/edge_truncate,/center)
end
