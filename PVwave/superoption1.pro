
function interpolation,courbe,delai       ;calcule une interpolation / extrapolation
zsize=size(courbe)
n=zsize(1)
if n eq 1 then return,courbe(0,1)
win=where(courbe(*,0) ge delai,ct)
if ct eq 0 then return,courbe(n-1,1)
winmin=min(win)
if winmin eq 0 then return,courbe(0,1)
return,(courbe(winmin,1)-courbe(winmin-1,1))*(delai-courbe(winmin-1,0))/$
	(courbe(winmin,0)-courbe(winmin-1,0))+courbe(winmin-1,1)
end

function intervalue1,courbe,delai1,delai2   	; calcule la valeur estime  ou zero
   						; pour avoir une valeur integree 
                                                ; exacte de courbe
						; (important pour les dividendes)
						; entre les bornes delai1 et delai2

zsize=size(courbe)
n=zsize(1)
t1=courbe(0,0)
t2=courbe(n-1,0)
sum=0.
delai1a=delai1>0
delai2a=delai2<t2
win1=where(courbe(*,0) ge delai1a,ct1)
n1=min(win1)
win2=where(courbe(*,0) ge delai2a,ct2)
n2=min(win2)
sum=sum+((courbe(n1,0)-delai1a)*courbe(n1,1))
sum=sum-((courbe(n2,0)-delai2a)*courbe(n2,1))
if n1 eq n2 then return,sum
for i=n1,n2-1 do sum=sum+((courbe(i+1,0)-courbe(i,0))*courbe(i+1,1))
return,sum/float(delai2-delai1)
end


;calcul d'une option
;americaine a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
;type=1 pour un call , 2 pour un put
	; ctaux et cvol sont les courbes de taux jj forward
	; et de volatilite implicite courte forward tel que calcule
	; par les programmes courbetaux et courbevol
	; cap est un strike qui cappe a la hausse ou a la baisse 
	; la valeur du support pouvant intervenir dans le calcul de l'option
;cvol est bien une volatilite instantanee
function option_binom,otype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist $
     ,ctaux=ctaux,cvol=cvol,cap=cap,courbdiv=courbdiv
	case 1 of
		otype eq 1 : spreadmax = 0
		otype eq 2 : spreadmax = 0
		else:begin
			print,'type d option non connu'
			stop
		end
	endcase
if keyword_set(cap) ne 0 then capflag = 1 else capflag = 0
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
value=fltarr(n+1)
value2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vdeltat=fltarr(n+1)
vp=fltarr(n+1)
vtaux=replicate(r,n)		; on cree vtaux 
vvol=replicate(sigma,n)		; et vvol par default
if (keyword_set(ctaux) ne 0)  then begin
	for i=0,n-1 do vtaux(i)=interpolation(ctaux,deltat*i)/100.
endif
if  (keyword_set(courbdiv) ne 0) then begin
  for i=0,n-1 do vtaux(i)=vtaux(i)-intervalue1(courbdiv,deltat*i,deltat*(i+1))
endif
if keyword_set(cvol) ne 0 then begin
	for i=0,n-1 do begin
		vvol(i)=interpolation(cvol,deltat*i)/100.
	endfor
endif
cste2=t/total(1/vvol(0:n-1)^2)
vdeltat(0:n-1)=cste2/(vvol(0:n-1)^2)
if keyword_set(cvol) ne 0 then begin
tvdeltat=0.
for i=0,n-1 do begin 
	vvol(i)=interpolation(cvol,tvdeltat)/100.
	tvdeltat=tvdeltat+vdeltat(i)
endfor
endif
cste2=t/total(1/vvol(0:n-1)^2)
vdeltat(0:n-1)=cste2/(vvol(0:n-1)^2)
a=exp(vtaux*vdeltat)
u=exp(vvol(1)*sqrt(vdeltat(1)))
d=1/u
vp=(a-d)/(u-d)
if keyword_set(cap) ne 0 then case 1 of
		otype eq 1 : spreadmax = (cap-k)>0.
		otype eq 2 : spreadmax = (k- cap)>0.
		else:print,'type d option non connu'
endcase else spreadmax=s*(u^n)
indexlist=where((t_yieldlist lt t) and (t_yieldlist gt 0),count)
if count gt 0 then ylist=yieldlist(indexlist)
pp=1.
if count gt 0 then for i=0,count-1 do pp=pp*(1-ylist(i))
 for i=0,n do begin
	supp(i)=s*(u^i)*(d^(n-i)) 
	case 1 of
		otype eq 1 :value(i)=(supp(i)*pp-k>0)<spreadmax
		otype eq 2 :value(i)=(k-(supp(i)*pp)>0)<spreadmax
		else:print,'type d option non connu'
	endcase
endfor
for j=0,n do value2(j)=0
tcalcul=t
for i=n-1,0,-1 do begin
 tcalcul=tcalcul-vdeltat(i)
 indexlist=where((t_yieldlist lt t) and (t_yieldlist gt 0),count)
 if count gt 0 then ylist=yieldlist(indexlist)
 pp=1.
 if count gt 0 then for ii=0,count-1 do pp=pp*(1-ylist(ii))
 if tcalcul gt tta then case 1 of
  otype eq 1 : for j=0,i do value2(j)=(s*u^j*d^(i-j)*pp-k> $
   exp(-vtaux(i)*vdeltat(i))*(vp(i)*value(j+1)+(1-vp(i))*value(j)))<spreadmax
  otype eq 2 : for j=0,i do value2(j)=(k-(s*u^j*d^(i-j)*pp)> $
   exp(-vtaux(i)*vdeltat(i))*(vp(i)*value(j+1)+(1-vp(i))*value(j)))<spreadmax
  else:print,'type d option non connu'
 endcase else begin
  for j=0,i do value2(j)=exp(-vtaux(i)*vdeltat(i))*$
	(vp(i)*value(j+1)+(1-vp(i))*value(j))
 endelse
 value(0:i)=value2(0:i)
endfor
return,value2(0)
end

function aj,d
;calcule nb d'annee a partir du jour jjmmaaaa ,
a=d mod 10000
m=((d-a) mod 1000000)/10000
j=(d-a-10000*m)/1000000
return,julday(m,j,a)/365.25
end


function courbetaux
t0=aj(20101993)
v=[	[aj(21101993)-t0,	3.1683],$
	[aj(22101993)-t0,	3.1683],$
	[aj(29101993)-t0,	3.1676],$
	[aj(22111993)-t0,	3.2238],$
	[aj(15121993)-t0,	3.2888],$
	[aj(16031994)-t0,	3.3964],$
	[aj(15061994)-t0,	3.4217],$
	[aj(21091994)-t0,	3.4861],$
	[aj(21121994)-t0,	3.5624],$
	[aj(15031995)-t0,	3.6701],$
	[aj(21061995)-t0,	3.7706],$
	[aj(20091995)-t0,	3.8637],$
	[aj(20121995)-t0,	3.9561],$
	[aj(20031996)-t0,	4.0594],$
	[aj(19061996)-t0,	4.1472],$
	[aj(18091996)-t0,	4.2327],$
	[aj(18121996)-t0,	4.3151],$
	[aj(19031997)-t0,	4.4037],$
	[aj(18061997)-t0,	4.4802],$
	[aj(17091997)-t0,	4.5552],$
	[aj(17121997)-t0,	4.6267],$
	[aj(18031998)-t0,	4.7025],$
	[aj(17061998)-t0,	4.7691],$
	[aj(16091998)-t0,	4.8340],$
	[aj(16121998)-t0,	4.8951],$
	[aj(22041999)-t0,	4.9310],$
	[aj(22101999)-t0,	5.0030],$
	[aj(24042000)-t0,	5.0764],$
	[aj(23102000)-t0,	5.1512],$
	[aj(23042001)-t0,	5.2397],$
	[aj(22102001)-t0,	5.3289],$
	[aj(22042002)-t0,	5.4208],$
	[aj(22102002)-t0,	5.5130],$
	[aj(22042003)-t0,	5.6081],$
	[aj(22102003)-t0,	5.7039]]
return,v
 end



function courbevol
t0=aj(20101993)
v=[	[aj(21101993)-t0,	3.1683],$
	[aj(22101993)-t0,	3.1683],$
	[aj(29101993)-t0,	3.1676],$
	[aj(22111993)-t0,	3.2238],$
	[aj(15121993)-t0,	3.2888],$
	[aj(16031994)-t0,	3.3964],$
	[aj(15061994)-t0,	3.4217],$
	[aj(21091994)-t0,	3.4861],$
	[aj(21121994)-t0,	3.5624],$
	[aj(15031995)-t0,	3.6701],$
	[aj(21061995)-t0,	3.7706],$
	[aj(20091995)-t0,	3.8637],$
	[aj(20121995)-t0,	3.9561],$
	[aj(20031996)-t0,	4.0594],$
	[aj(19061996)-t0,	4.1472],$
	[aj(18091996)-t0,	4.2327],$
	[aj(18121996)-t0,	4.3151],$
	[aj(19031997)-t0,	4.4037],$
	[aj(18061997)-t0,	4.4802],$
	[aj(17091997)-t0,	4.5552],$
	[aj(17121997)-t0,	4.6267],$
	[aj(18031998)-t0,	4.7025],$
	[aj(17061998)-t0,	4.7691],$
	[aj(16091998)-t0,	4.8340],$
	[aj(16121998)-t0,	4.8951],$
	[aj(22041999)-t0,	4.9310],$
	[aj(22101999)-t0,	5.0030],$
	[aj(24042000)-t0,	5.0764],$
	[aj(23102000)-t0,	5.1512],$
	[aj(23042001)-t0,	5.2397],$
	[aj(22102001)-t0,	5.3289],$
	[aj(22042002)-t0,	5.4208],$
	[aj(22102002)-t0,	5.5130],$
	[aj(22042003)-t0,	5.6081],$
	[aj(22102003)-t0,	5.7039]]
return,v
 end


pro test
otype1=1
otype2=2
s=100.
k=100.
t=1.
tta=0.
r=0.07
sigma=0.2
n=50
yieldlist=[0.02,0.02]
t_yieldlist=[1.,1.5]
print,option_binom(otype1,s,k,t,tta,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol()),$
           option_binom(otype2,s,k,t,tta,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol())
print,option_binom(otype1,s,k,t,t,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol()),$
            option_binom(otype2,s,k,t,t,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol())
print,option_binom(otype1,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,cap=k*1.1),$
         option_binom(otype2,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,cap=k*0.9)
print,option_binom(otype1,s,k,t,t,r,sigma,n,yieldlist,t_yieldlist,cap=k*1.1),$
      option_binom(otype2,s,k,t,t,r,sigma,n,yieldlist,t_yieldlist,cap=k*0.9)
end



