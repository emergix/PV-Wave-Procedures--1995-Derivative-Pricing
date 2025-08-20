function defaultbond,q,f,r,n
w=(1-q)/(1+r)
if q ne 0 then $
 s=f*(q/r)*(((1-q)^(n+1)-1)/q-(w^(n+1)-1.)/(1-w))+$
	(1-q)^n*(1/(1+r)^n)*(f*((1+r)^n-1.)/r+1.) $
else s=(f*((1+r)^n-1.)/r+1.)/(1+r)^n
return,s
end

function bond,s,f,r,n
return,((f/(r+s))*((1+r+s)^n-1.)+1.)/(1+r+s)^n
end

pro display_courbe,x,y,titre,xtitre,ytitre
	set_plot,'PS'
	device,/landscape
	plot,x,y,xtitle=xtitre ,ytitle=ytitre,title=titre
	device,/close_file
	spawn,'lpr wave.ps'
	set_plot,'X'
end

pro slides1
n=10
r=0.05
f=0.04
q=findgen(100)/100*0.1
bd=fltarr(100)
for i=0,99 do bd(i)=defaultbond(q(i),f,r,n)
display_courbe,q,bd,'bonds with default','default prob.','bond price'
end

pro slides2
n=10
r=0.05
f=0.04
s=findgen(100)/100*0.1
bd=fltarr(100)
for i=0,99 do bd(i)=bond(s(i),f,r,n)
display_courbe,s,bd,'bonds with spread','spread','bond price'
end





function implicit_spread,q,f1,r,n,x1=x1,x2=x2,xacc=xacc,jmax=jmax
prime=defaultbond(q,f1,r,n)
if keyword_set(x1) eq 0 then x1=0.00000001
if keyword_set(x2) eq 0 then x2=0.05
if keyword_set(xacc) eq 0. then xacc=0.000001
if keyword_set(jmax) eq 0. then jmax=40
fmid=bond(x2,f1,r,n)-prime
f=bond(x1,f1,r,n)-prime
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
fmid=bond(x2,f1,r,n)-prime
	id=id+1
endwhile
if f*fmid ge 0 then begin
	print,'x1 tet x2 n encadre pas la racine : x1=',x1,' x2=',x2
	return,0
	endif
if f lt 0 then begin
	rtbis=x1
	dx=x2-x1
endif else begin
	rtbis=x2
	dx=x1-x2
endelse
while 	abs(dx) ge xacc or fmid ne 0 do begin
	dx=dx*0.5
	xmid=rtbis+dx
fmid=bond(xmid,f1,r,n)-prime
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end


pro slides3
f=0.04
r=0.05
n=10
q=(findgen(100)+1)/100*0.1
bd=fltarr(100)
for i=0,99 do bd(i)=(implicit_spread(q(i),f,r,n)-q(i))/q(i)
display_courbe,q,bd,'bonds with default','default prob.','implicit spread-default prob in %'
end

pro slides4
f=(findgen(100)+1)/100.*0.1
r=0.05
n=10
q= 0.01
bd=fltarr(100)
for i=0,99 do bd(i)=implicit_spread(q,f(i),r,n)
display_courbe,f*100,bd,'bonds with default','fixed rate in %','implicit spread'
end

pro slides5
f=0.06
r=0.05
n=findgen(30)+1
q= 0.01
bd=fltarr(30)
for i=0,29 do bd(i)=implicit_spread(q,f,r,n(i))
display_courbe,n,bd,'bonds with default fixed rate = 6%','maturity','implicit spread'
end
