
 
function CRYLcall,s,k,t,r,sigma,n,yieldlist,t_divlist
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
u=exp(sigma*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
indexlist=where(t_divlist lt t,count)
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=supp(i)-k>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=supp(i)*pp-k>0
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		for j=0,i do call2(j)=s*u^j*d^(i-j)-k> $
	exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		for j=0,i do call2(j)=s*u^j*d^(i-j)*pp-k> $
	exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end

pro test
s=100.
k=100.
t=2.
r=0.09
sigma=0.20
n=30
t_divlist=[0.5,1.,1.5]          ;     doit avoir 
yieldlist=[0.02,0.02,0.02]      ;     le meme nombre d'elements
print,CRYLcall(s,k,t,r,sigma,n,yieldlist,t_divlist)
end

