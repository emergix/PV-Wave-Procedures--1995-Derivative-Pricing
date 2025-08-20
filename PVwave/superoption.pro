function interpolation,courbe,delai  ;calcule une interpolation / extrapolation
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
   						; pour avoir une valeur integree exacte
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


function moyenne_courbe,courbe,delai
zsize=size(courbe)
n=zsize(1)
t1=courbe(0,0)
f1=courbe(0,1)
if delai lt t1 then return,f1
sigma_s=f1*t1
i=1
while i lt n  do  begin
	if courbe(i,0) ge delai then begin
		sigma_s=sigma_s+(delai-t1)*(2*f1+(courbe(i,1)-f1)*(delai-t1)/(courbe(i,0)-t1))/2.
		return,sigma_s/delai
		endif
	sigma_s=sigma_s+(courbe(i,0)-t1)*(courbe(i,1)+f1)/2.
	t1=courbe(i,0)
	f1=courbe(i,1)
	i=i+1	
endwhile
sigma_s=sigma_s+(delai-t1)*f1
return,sigma_s/delai
end

function courbe_integrale,courbe,t0,x0
;retourne une courbe representant l integrale de courbe valant x0 en t0
zsize=size(courbe)
n=zsize(1)
x=fltarr(n,2)
x(0,0)=courbe(0,0)
x(0,1)=0
for j= 0,n-2 do begin
	x(j+1,0) = courbe(j+1,0)
	x(j+1,1) = courbe(j,1)*(courbe(j+1,0)-courbe(j,0))+x(j,1)
endfor
x0p=interpolation(x,t0)
dx=x0-x0p
for j=0,n-1 do x(j,1)=x(j,1)+dx
return,x
end

function courbe_liss_integrale,courbe,x0
;retourne une courbe representant la 1/t*Somme(u=0,u=t,courbe(u)) et valant x0 en 0
zsize=size(courbe)
n=zsize(1)
x=fltarr(n,2)
x(0,0)=courbe(0,0)
x(0,1)=x0
for j=0,n-2 do begin
	x(j+1,0) = courbe(j+1,0)
	x(j+1,1) = (courbe(j,1)*(courbe(j+1,0)-courbe(j,0))+(x(j,1)*(x(j,0)-x(0,0))))/(x(j+1,0)-x(0,0))
endfor
return,x
end


function courbe_derivee,courbe
;retourne une courbe representant la derivee de courbe
zsize=size(courbe)
n=zsize(1)
y=fltarr(n,2)
for j=0,n-2 do begin
	y(j,0)=courbe(j,0)
	y(j,1)=(courbe(j+1,1)-courbe(j,1))/(courbe(j+1,0)-courbe(j,0))
endfor
y(n-1,0)=courbe(n-1,0)
y(n-1,1)=y(n-2,1)
return,y
end

pro test_d
a=transpose([[0,0],[3,0.2 ]])
print,a
b=courbe_derivee(a)
print,b
end

function courbe_log,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=alog(y(i,1))
return,y
end

function courbe_exp,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=exp(y(i,1))
return,y
end


function courbe_neg,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=-y(i,1)
return,y
end

function courbe_cmul,c,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=c*y(i,1)
return,y
end

function courbe_cplus,c,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=c+y(i,1)
return,y
end


function vect_union,vect1,vect2
n1=n_elements(vect1)
n2=n_elements(vect2)
i1=0
i2=0
vect=fltarr(n1+n2)
vect(0:n1-1)=vect1
vect(n1:n1+n2-1)=vect2
as=sort(vect)
az=vect(as)
ind=1
maxi=az(0)
vect(0)=maxi
for i=1,n1+n2-1 do begin
if  az(i) gt maxi then begin
		vect(ind)=az(i)
		ind=ind+1
		maxi=az(i)
		endif
endfor
a=fltarr(ind)
a(0:ind-1)=vect(0:ind-1)
return,vect(0:ind-1)
end




function courbe_plus,courbe1,courbe2
date1=courbe1(*,0)
date2=courbe2(*,0)
date=vect_union(date1,date2)
n=n_elements(date)
result=fltarr(n,2)
for i=0,n-1 do begin
	result(i,0)=date(i)
	result(i,1)=interpolation(courbe1,date(i))+interpolation(courbe2,date(i))
endfor
return,result
end

function courbe_diff,courbe1,courbe2
date1=courbe1(*,0)
date2=courbe2(*,0)
date=vect_union(date1,date2)
n=n_elements(date)
result=fltarr(n,2)
for i=0,n-1 do begin
	result(i,0)=date(i)
	result(i,1)=interpolation(courbe1,date(i))-interpolation(courbe2,date(i))
endfor
return,result
end

function courbe_mul,courbe1,courbe2
date1=courbe1(*,0)
date2=courbe2(*,0)
date=vect_union(date1,date2)
n=n_elements(date)
result=fltarr(n,2)
for i=0,n-1 do begin
	result(i,0)=date(i)
	result(i,1)=interpolation(courbe1,date(i))*interpolation(courbe2,date(i))
endfor
return,result
end

function courbe_div,courbe1,courbe2
date1=courbe1(*,0)
date2=courbe2(*,0)
date=vect_union(date1,date2)
n=n_elements(date)
result=fltarr(n,2)
for i=0,n-1 do begin
	result(i,0)=date(i)
	result(i,1)=interpolation(courbe1,date(i))/interpolation(courbe2,date(i))
endfor
return,result
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
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol()),option_binom(otype2,s,k,t,tta,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol())
print,option_binom(otype1,s,k,t,t,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol()),option_binom(otype2,s,k,t,t,r,sigma,n,$
	yieldlist,t_yieldlist,ctaux=courbetaux(),cvol=courbevol())
print,option_binom(otype1,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,cap=k*1.1),option_binom(otype2,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,cap=k*0.9)
print,option_binom(otype1,s,k,t,t,r,sigma,n,yieldlist,t_yieldlist,cap=k*1.1),option_binom(otype2,s,k,t,t,r,sigma,n,yieldlist,t_yieldlist,cap=k*0.9)
end

;**********************************************************************************


;courbepd   courbe des  prix des zero coupons domestiques
;courbepf   courbe des  prix des zero coupons  etrangers
;courbevar_s  courbe des variance des taux de change 1 etranger= Sd * domestique
;courbevar_d  courbe des variance des zero coupons domestiques
;courbevar_f  courbe des variance des zero coupons etrangers
;courbecov_sd  courbe des covariances taux de change spot avec zero coupons a terme domestique
;courbecov_sf  courbe des covariances taux de change spot avec zero coupons a terme etranger
;courbecov_df  courbe des covariances zero coupons domestique avec zero coupons etrangers

function currency_option,courbepd,courbepf,courbevar_s,courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
		otype,Sd0,k,T,n
	case 1 of
		otype eq 1 : spreadmax = 0
		otype eq 2 : spreadmax = 0
		else:begin
			print,'type d option non connu'
			stop
		end
	endcase
ctauxd=courbe_neg(courbe_derivee(courbe_log(courbepd)))
ctauxf=courbe_neg(courbe_derivee(courbe_log(courbepf)))
ctaux=courbe_diff(ctauxd,ctauxf)
sig_d= courbevar_d
sig_f= courbevar_f
sig_s= courbevar_s
cov_df= courbecov_df
cov_sf= courbecov_sf
cov_sd= courbecov_sd
integrand=courbe_plus(sig_s,courbe_plus(sig_d,courbe_plus(sig_f,courbe_cmul(2,cov_sf))))
cvol=courbe_diff(integrand,courbe_plus(courbe_cmul(2,cov_sd),courbe_cmul(2,cov_df)))
yield_list=[0]
t_yieldlist=[T]
tta=0.
price=option_binom(otype,Sd0,k,T,tta,0.0,0.,n,[0,0],[T/2.,T],ctaux=ctaux,cvol=cvol)
return,price
end


pro test_currency
courbepd=transpose($
	[	[0,	1.],$
		[1.,	.96],$
		[2.,	.92],$
		[3.,	.88]])

courbepf=transpose($
	[	[0,	1.],$
		[1.,	.95],$
		[2.,	.91],$
		[3.,	.87]])

courbevar_s=transpose($
	[	[0,	.03^2],$
		[1.,	.03^2],$
		[2.,	.03^2],$
		[3.,	.03^2]])

courbevar_f=transpose($
	[	[0,	.03^2],$
		[1.,	.03^2],$
		[2.,	.03^2],$
		[3.,	.03^2]])

courbevar_d=transpose($
	[	[0,	.03^2],$
		[1.,	.03^2],$
		[2.,	.03^2],$
		[3.,	.03^2]])


courbecov_sd=transpose($
	[	[0,	0.],$
		[1.,	.02^2],$
		[2.,	.02^3*2],$
		[3.,	.02^2*3]])


courbecov_sf=transpose($
	[	[0,	0.],$
		[1.,	.01^2],$
		[2.,	.01^3*2],$
		[3.,	.01^2*3]])


courbecov_df=transpose($
	[	[0,	0.],$
		[1.,	.01^2],$
		[2.,	.01^3*2],$
		[3.,	.01^2*3]])
otype=1
Sd0=100.
k=100.
T=2.
n=100
opt=currency_option(courbepd,courbepf,courbevar_s,courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
		otype,Sd0,k,T,n)
print,opt
end



;**********************************************************************************
 
 
function find_correlations,courbepd,courbepf,courbevar_s,courbevar_f,courbevar_d,market_data,n
;all the curves have to have the same x points for the dates this is assumed
; the format for market_data is :
;  [ [otype,Sd0,k,T,price] , ..[] ]
; 
zsize=size(courbepd)
nd=zsize(1)
date1=courbepd(*,0)
maxt=max(date1)
msize=size(market_data)
nm=msize(1)
points=fltarr(nm)
 ; let assume that the correlations are constant with the time
; the unkowns are asd,asf,adf : 3 variables
;  courbecov_sd=asd 
;  courbecov_sf=asf
;  courbecov_df=adf 
ndim=3
psum=fltarr(ndim+1)
ptry=fltarr(ndim+1)
; algorithm that regress on this variables
moy=fltarr(ndim+1)
moys=moyenne_courbe(courbevar_s,maxt)
moyd=moyenne_courbe(courbevar_d,maxt)
moyf=moyenne_courbe(courbevar_f,maxt)
moy(1)=sqrt(moys*moyd)
moy(2)=sqrt(moys*moyf)
moy(3)=sqrt(moyd*moyf)
print,'cov_sd=',moy(1),'cov_sf=',moy(2),'cov_df=',moy(3)
y=fltarr(ndim+2)
p=fltarr(ndim+2,ndim+1)
p(0,1)=1*moy(1)	;initialisation de cov_sd
p(0,2)=0*moy(2)	;initialisation de cov_sf
p(0,3)=0*moy(3)	;initialisation de cov_df
for i=1,ndim do for j=1,ndim do  if i eq j then p(i+1,j)=p(0,j)+moy(i) else p(i+1,j)=0.
for j=1,ndim+1 do begin
	ptry(1:ndim)=p(j,1:ndim)
	;	y(j)=funk(ptry) 					appel a funk
	;#################     definition of funk(ptry)
	asd=ptry(1)
	asf=ptry(2)
	adf=ptry(3)
	courbecov_sd=transpose([[0,asd],[maxt,asd]])
	courbecov_sf=transpose([[0,asf],[maxt,asf]])
	courbecov_df=transpose([[0,adf],[maxt,adf]])
	for i=0,nm-1 do  points(i)=  currency_option(courbepd,courbepf,courbevar_s,$
		courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
		market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
	sum=0
	  for i=0,nm-1 do sum=sum+ (points(i)-market_data(i,4))^2 
	;################### the result is sum
	y(j)=sum
endfor
print,'y initial=',y
ftol=0.0001
;********************* implementation de amoeba p 411 de numerical recipes
nfunk=0
mpts=ndim+1
nmax=5000
 ;get_psum
for j=1,ndim do begin
	sum=0.
	for i=1,mpts do begin
		sum=sum+p(i,j)
		psum(j)=sum
	endfor
endfor
;
endcalcul=0
while (endcalcul eq 0) do begin
	ilo=1
	if y(1) gt y(2) then begin
		inhi=2
		ihi=1
	endif else begin
		inhi=1
		ihi=2
	endelse
	for i=1,mpts do begin
		if y(i) le y(ilo) then ilo=i
	 	if y(i) gt y(ihi) then begin
			inhi=ihi
			ihi=i
		endif else if ((y(i) gt y(inhi)) and (i ne ihi)) then inhi=i
	endfor
	rtol= 2.*abs(y(ihi)-y(ilo)) 
print,points
print,'but=',ftol,'actual=',rtol,'sol=',y(ihi),y(ilo)
	if rtol lt ftol then  begin
		result=fltarr(ndim)
		result(0:ndim-1)=p(ilo,1:ndim)
		endcalcul=1
	endif
if endcalcul eq 0 then begin
	if nfunk gt nmax then begin
		print,'nmax exceeded'
		return,y
	endif
	nfunk=nfunk+2 
	; amotry
	fac=-1.
	fac1=(1.-fac)/float(ndim)
	fac2=fac1-fac
	for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
	;ytry=funk(ptry)								appel a funk
	;#################     definition of funk(ptry)
	asd=ptry(1)
	asf=ptry(2)
	adf=ptry(3)
	courbecov_sd=transpose([[0,asd],[maxt,asd]])
	courbecov_sf=transpose([[0,asf],[maxt,asf]])
	courbecov_df=transpose([[0,adf],[maxt,adf]])
	for i=0,nm-1 do  points(i)=  currency_option(courbepd,courbepf,courbevar_s,$
		courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
		market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
	sum=0
	  for i=0,nm-1 do sum=sum+ (points(i)-market_data(i,4))^2 
	;################### the result is sum
	ytry=sum
	if  ytry lt y(ihi) then begin
		y(ihi)=ytry
		for j=1,ndim do begin
			psum(j)=psum(j)+ptry(j)-p(ihi,j)
			p(ihi,j)=ptry(j)
		endfor
	endif
	if ytry lt y(ilo) then begin
		; amotry
		fac=2.
		ptry=fltarr(ndim+1)
		fac1=(1.-fac)/float(ndim)
		fac2=fac1-fac
		for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
		;	ytry=funk(ptry)							;appel a funk
		;#################     definition of funk(ptry)
		asd=ptry(1)
		asf=ptry(2)
		adf=ptry(3)
		courbecov_sd=transpose([[0,asd],[maxt,asd]])
		courbecov_sf=transpose([[0,asf],[maxt,asf]])
		courbecov_df=transpose([[0,adf],[maxt,adf]])
		for i=0,nm-1 do  points(i)=  currency_option(courbepd,courbepf,courbevar_s,$
			courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
			market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
		sum=0
		  for i=0,nm-1 do sum=sum+ (points(i)-market_data(i,4))^2 
		;################### the result is sum
		ytry=sum
		if  ytry lt y(ihi) then begin
			y(ihi)=ytry
			for j=1,ndim do begin
				psum(j)=psum(j)+ptry(j)-p(ihi,j)
				p(ihi,j)=ptry(j)
			endfor
		endif
		;
	endif else begin
		if ytry ge y(inhi) then begin
			ysave=y(ihi)
			; amotry
			fac=0.5
			fac1=(1.-fac)/float(ndim)
			fac2=fac1-fac
			for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
			;ytry=funk(ptry)						;appel a funk
			;#################     definition of funk(ptry)
			asd=ptry(1)
			asf=ptry(2)
			adf=ptry(3)
			courbecov_sd=transpose([[0,asd],[maxt,asd]])
			courbecov_sf=transpose([[0,asf],[maxt,asf]])
			courbecov_df=transpose([[0,adf],[maxt,adf]])
			for i=0,nm-1 do  points(i)=  currency_option(courbepd,courbepf,courbevar_s,$
				courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
				market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
			sum=0
			  for i=0,nm-1 do sum=sum+ (points(i)-market_data(i,4))^2 
			;################### the result is sum
			ytry=sum
			if  ytry lt y(ihi) then begin
				y(ihi)=ytry
				for j=1,ndim do begin
					psum(j)=psum(j)+ptry(j)-p(ihi,j)
					p(ihi,j)=ptry(j)
				endfor
			endif
			;
			if ytry ge ysave then begin
			 	for i1=1,mpts do if i1 ne ilo then begin
				  for j=1,ndim do begin
					psum(j)=0.5*(p(i1,j)+p(ilo,j)) 
					p(i1,j)=psum(j)
				  endfor
				  ;y(i1)=funk(psum)		;appel a funk
				  ;#################     definition of funk(ptry)
				  asd=ptry(1)
				  asf=ptry(2)
				  adf=ptry(3)
				  courbecov_sd=transpose([[0,asd],[maxt,asd]])
				  courbecov_sf=transpose([[0,asf],[maxt,asf]])
				  courbecov_df=transpose([[0,adf],[maxt,adf]])
				    for i=0,nm-1 do  points(i)=  currency_option(courbepd,courbepf,courbevar_s,$
					courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
					market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
				  sum=0
				  for i=0,nm-1 do sum=sum+ (points(i)-market_data(i,4))^2 
				  ;################### the result is sum
				y(i1)=sum
				endif
				;
				nfunk=nfunk+ndim
				 ;get_psum
				for j=1,ndim do begin
					sum=0.
					for i=1,mpts do begin
						sum=sum+p(i,j)
						psum(j)=sum
					endfor
				endfor
				;
			endif
		endif else nfunk=nfunk-2
	endelse
endif
endwhile
asd=result(0)
asf=result(1)
adf=result(2)
courbecov_sd=transpose([[0,asd],[maxt,asd]])
courbecov_sf=transpose([[0,asf],[maxt,asf]])
courbecov_df=transpose([[0,adf],[maxt,adf]])
correlations=fltarr(3,2,2)
correlations(0,*,*)=courbecov_sd
correlations(1,*,*)=courbecov_sf
correlations(2,*,*)=courbecov_df
return,correlations
end


pro test_covariance
courbepd=transpose($
	[	[0,	1.],$
		[1.,	.96],$
		[2.,	.92],$
		[3.,	.88]])

courbepf=transpose($
	[	[0,	1.],$
		[1.,	.95],$
		[2.,	.91],$
		[3.,	.87]])
 
courbevar_s=transpose($
	[	[0,	.1],$
		[1.,	.2],$
		[2.,	.3],$
		[3.,	.4]])

courbevar_f=transpose($
	[	[0,	.1],$
		[1.,	.1],$
		[2.,	.1],$
		[3.,	.1]])

courbevar_d=transpose($
	[	[0,	.1],$
		[1.,	.1],$
		[2.,	.1],$
		[3.,	.1]])

market_data=transpose([$
		[1,	100.,	100.,	1.,	0.214362 ],$
		[2,	100.,	100.,	1.,	0.220221 ],$
		[1,	100.,	100.,	2.,	0.331588 ],$
		[2,	100.,	100.,	2.,	0.339117 ]])

correlations=find_correlations(courbepd,courbepf,courbevar_s,courbevar_f,courbevar_d,market_data,30)
courbecov_sd=fltarr(2,2)
courbecov_sf=fltarr(2,2)
courbecov_df=fltarr(2,2)
courbecov_sd(*,*)=correlations(0,*,*)
courbecov_sf(*,*)=correlations(1,*,*)
courbecov_df(*,*)=correlations(2,*,*)

;window,1,title="courbecov_sd"
;plot,courbecov_sd(*,0),courbecov_sd(*,1)
;window,2,title="courbecov_sf"
;plot,courbecov_sf(*,0),courbecov_sf(*,1)
;window,3,title="courbecov_df"
;plot,courbecov_df(*,0),courbecov_df(*,1)
end


pro test_marche
courbepd=transpose($
	[	[0,	1.],$
		[1.,	.96],$
		[2.,	.92],$
		[3.,	.88]])

courbepf=transpose($
	[	[0,	1.],$
		[1.,	.95],$
		[2.,	.91],$
		[3.,	.87]])
courbevar_s=transpose($
	[	[0,	.1],$
		[1.,	.2],$
		[2.,	.3],$
		[3.,	.4]])

courbevar_f=transpose($
	[	[0,	.1],$
		[1.,	.1],$
		[2.,	.1],$
		[3.,	.1]])

courbevar_d=transpose($
	[	[0,	.1],$
		[1.,	.1],$
		[2.,	.1],$
		[3.,	.1]])

market_data=transpose([$
		[1,	100.,	100.,	1.,	0.214362 ],$
		[2,	100.,	100.,	1.,	0.220221 ],$
		[1,	100.,	100.,	2.,	0.331588 ],$
		[2,	100.,	100.,	2.,	0.339117 ]])
nm=4
n=30
moy=fltarr(3+1)
points=fltarr(nm)
date1=courbepd(*,0)
maxt=max(date1)
moys=moyenne_courbe(courbevar_s,maxt)
moyd=moyenne_courbe(courbevar_d,maxt)
moyf=moyenne_courbe(courbevar_f,maxt)
moy(1)=sqrt(moys*moyd)
moy(2)=sqrt(moys*moyf)
moy(3)=sqrt(moyd*moyf)
asd=1*moy(1)	;initialisation de cov_sd
asf=1*moy(2)	;initialisation de cov_sf
adf=-1*moy(3)	;initialisation de cov_df
print,'cov_sd=',asd,' / cov_sf=',asf,' / cov_df=',adf
courbecov_sd=transpose([[0,asd],[maxt,asd]])
courbecov_sf=transpose([[0,asf],[maxt,asf]])
courbecov_df=transpose([[0,adf],[maxt,adf]])
for i=0,nm-1 do  points(i)=  currency_option(courbepd,courbepf,courbevar_s,$
		courbevar_d,courbevar_f,courbecov_sd,courbecov_sf,courbecov_df,$
		market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
print,'prix de marche=',points
end
