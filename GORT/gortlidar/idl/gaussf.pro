function gaussf,x,a0,a1,a2

;       x - array of x values
;       a0 - amplitude of gaussian
;       a1 - center of gaussian
;       a2 - width of gaussian
;       g  - gaussian array

n=n_elements(x)
g=dblarr(n)

for i=0,n-1 do begin
	z=(x(i)-a1)/a2
	;a0=1.0/(a2*sqrt(2.0*!dpi))
	g(i)=a0*exp(-z^2/2.0)
endfor
return,g
end
