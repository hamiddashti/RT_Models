function gsmooth2,wave,width,mean,sd,SHIFT=shift

mult = 4  

N  = n_elements(wave)
N2 = fix(width*mult)
if(N2 gt N/2) then N2 = N/2

padlen = N2
pad = fltarr(padlen) + mean

N2 = fix(N2/2)*2 + 1
shiftval = 0
if(keyword_set(SHIFT)) then shiftval = shift
g1 = gaussf(findgen(N2),1.0,float(fix(N2/2))+shiftval,float(width))
g1 = g1/float(total(g1))
g1 = reverse(g1)

cwave = fltarr(N+2*padlen)
cwave = convol([pad,float(wave),pad],g1,/edge_truncate,/center)
return,cwave[padlen:N+padlen-1]

end
