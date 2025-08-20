function R_alpha,sigma,delta,gamma,R
x1=-delta/gamma-sqrt((delta/gamma)^2-R/gamma)
x2=-delta/gamma+sqrt((delta/gamma)^2-R/gamma)
return,1-(gaussint(x2/sigma)-gaussint(x1/sigma))
end

function call,s,k,t,r,sigm
if (t le 0)  then return,0>(s-k) else begin
	if (sigm eq 0.0) then return,(s-k*exp(-r*t))>0 else begin
		d1=(alog(s/float(k))+(r+sigm*sigm/2)*t)/(sigm*sqrt(t))
		d2=d1-(sigm*sqrt(t))
		return,s*gaussint(d1)-k*gaussint(d2)*exp(-r*t)
	endelse
endelse
end

function cdelta,s,k,t,r,sigm
if (t le 0.)  then begin
	if (s gt k) then  return,1. else return,0
endif else begin
	if (sigm eq 0.0) then return,1.0
	d1=(alog(s/k)+(r+(sigm*sigm)/2)*t)/(sigm*sqrt(t))
	return,gaussint(d1)
endelse
end

function crho,s,k,t,r,sigm
 d2=(alog(s/k)+(r-(sigm*sigm)/2)*t)/(sigm*sqrt(t))
return,k*t*gaussint(d2)*exp(-r*t)
end

function cgamma,s,k,t,r,sigm
if t le 0 or sigm eq 0. then return,0.
d1=(alog(s/float(k))+(r+sigm*sigm/2)*t)/(sigm*sqrt(t))
return,1/(double(2.50662827463)*sqrt(t)*sigm*exp(d1^2/2)*s)
end

Pro test1
window
for k=92,108 do test,k
end

Pro test,k
x=findgen(100)/100.
s=100
texp=5/365.
r=0.05
vol=0.2
t=1/365.
gamma=-cgamma(s,k,texp,r,vol)
delta=s*gamma+cdelta(s,k,texp,r,vol)
;delta=s*gamma
sigma=vol*sqrt(t)
print,'sigma=',sigma
y=fltarr(100)
for i=0,99 do y(i)=100.*r_alpha(sigma,delta,gamma,x(i))
print,y(50)
oplot,x,y,xtitle='risk ($)',ytitle='loss probability (%)',yrange=[0,20]
end

pro tst
s=1


  
k=100.
t=5/365.
r=0.05
sigm=0.2
print,'call=',call(s,k,t,r,sigm)
print,'delta=', cdelta(s,k,t,r,sigm)
print,'gamma=',cgamma(s,k,t,r,sigm)
end
