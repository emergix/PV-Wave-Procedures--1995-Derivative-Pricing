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

function courbe_annualisee,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-2 do begin
	y(i,1)=y(i,1)/(y(i+1,0)-y(i,0))
endfor
y(n-1,1)=y(n-2,1)
return,y
end


function courbe_cmul,c,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=c*y(i,1)
return,y
end

function courbe_floor,floor,courbe
zsize=size(courbe)
n=zsize(1)
y=courbe
for i=0,n-1 do y(i,1)=y(i,1)>floor
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

 ;**********************************************************************************

 


;**********************************************************************************
 
;courbepd   courbe des  prix des zero coupons domestiques
;courbepf   courbe des  prix des zero coupons  etrangers
; courbezeta est la volatilite resultante 1/t * Sum(0,t,Vs^2+Vd^2+Vf^2+2COVsf -2COVsd -2COVdf) annualisee
function currency_option1,courbepd,courbepf,courbe_zeta,$
		otype,Sd0,k,T,n
	case 1 of
		otype eq 1 : spreadmax = 0	; call americain
		otype eq 2 : spreadmax = 0	; put americain
		otype eq 3 : spreadmax = 0	; call europeeen
		otype eq 4 : spreadmax = 0	; put europeen 

		else:begin
			print,'type d option non connu'
			stop
		end
	endcase
ctauxd=courbe_neg(courbe_derivee(courbe_log(courbepd)))
ctauxf=courbe_neg(courbe_derivee(courbe_log(courbepf)))
ctaux=courbe_cmul(100,courbe_annualisee(courbe_diff(ctauxd,ctauxf)))
tta=0.
price=option_binom(otype,Sd0,k,T,tta,0.0,0.,n,[0,0],[T/2.,T],ctaux=ctaux,cvol=courbe_zeta)
return,price
end
  

pro test_marche
courbepd=transpose($
	[	[0,	1.],$
		[1.,	.95],$
		[2.,	.91],$
		[3.,	.864]])

courbepf=transpose($
	[	[0,	1.],$
		[1.,	.97],$
		[2.,	.85],$
		[3.,	.721]])
courbe_zeta=transpose($
	[	[0,	5.],$
		[1.,	4.],$
		[2.,	3.],$
		[3.,	2.]])

 
market_data=transpose([$
		[1,	100.,	100.,	1.,	 1.38297],$
		[1,	100.,	102.,	1.,	 0.618912],$
		[2,	100.,	100.,	1.,	 4.21440],$
		[1,	100.,	100.,	2.,	 1.30657],$
		[1,	100.,	102.,	2.,	 0.552642],$
		[2,	100.,	100.,	2.,	 16.2620],$
		[1,	100.,	100.,	3.,	 1.23729],$
		[1,	100.,	102.,	3.,	 0.490455],$
		[2,	100.,	100.,	3.,	 32.2448],$
		[2,	100.,	90.,	3.,	 32.2448]])

msize=size(market_data)
nm=msize(1)
n=30
points=fltarr(nm)
date1=courbepd(*,0)
maxt=max(date1)
for i=0,nm-1 do  points(i)=  currency_option1(courbepd,courbepf,courbe_zeta,$
 		market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
print,'prix de marche=',points
end

 
function find_zeta ,courbepd,courbepf, zeta0,datelist,market_data,n
;all the curves have to have the same x points for the dates this is assumed
; the format for market_data is :
;  [ [otype,Sd0,k,T,price] , ..[] ]
; 
; zeta 0  : variance implicite initale 
zsize=size(courbepd)
nd=zsize(1)
date1=courbepd(*,0)
maxt=max(date1)
maxd=max(datelist)
ndim=n_elements(datelist)
if maxd gt maxt then begin
	print,'echelle de date dans datelist trop etendue '
	stop
endif
msize=size(market_data)
nm=msize(1)
points=fltarr(nm)
courbe_zeta=dblarr(ndim,2)
courbe_zeta(*,0)=datelist
psum=fltarr(ndim+1)
ptry=fltarr(ndim+1)
; algorithm that regress on this variables
y=fltarr(ndim+2)
p=fltarr(ndim+2,ndim+1)
p(0,*)=alog(zeta0)
for i=1,ndim do for j=1,ndim do  if i eq j then p(i+1,j)=alog(exp(p(0,j))+zeta0/2.1) else p(i+1,j)=alog(p(0,j))
for j=1,ndim+1 do begin
	ptry(1:ndim)=p(j,1:ndim)
	;	y(j)=funk(ptry) 					appel a funk
	;#################     definition of funk(ptry)
	courbe_zeta(*,1)=ptry(1:ndim)
	courbe_zeta=courbe_exp(courbe_zeta)
	for i=0,nm-1 do  points(i)=  currency_option1(courbepd,courbepf,courbe_zeta,$
		market_data(i,0),market_data(i,1),market_data(i,2),market_data(i,3),n)
	sum=0
	  for i=0,nm-1 do sum=sum+ (points(i)-market_data(i,4))^2 
	;################### the result is sum
	y(j)=sum
endfor
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
print,'theoretical prices=',points(1:*)
print,'but=',ftol,' / actual=',rtol,'  / y plus Haut ,y plus Bas =',y(ihi),y(ilo)
print,'best zeta=',exp(transpose(p(ilo,1:*)))
	if rtol lt ftol then  begin
		print,'optimum atteint with ',nfunk,' calls to the evaluation function'
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
	; 										 amotry
	fac=-1.
	fac1=(1.-fac)/float(ndim)
	fac2=fac1-fac
	for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
	;ytry=funk(ptry)								appel a funk
	;#################     definition of funk(ptry)
	courbe_zeta(*,1)=ptry(1:ndim)
	courbe_zeta=courbe_exp(courbe_zeta)
	for i=0,nm-1 do  points(i)=  currency_option1(courbepd,courbepf,courbe_zeta,$
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
		; 									amotry
		fac=2.
		ptry=fltarr(ndim+1)
		fac1=(1.-fac)/float(ndim)
		fac2=fac1-fac
		for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
		;	ytry=funk(ptry)							;appel a funk
		;#################     definition of funk(ptry)
		courbe_zeta(*,1)=ptry(1:ndim)
		courbe_zeta=courbe_exp(courbe_zeta)
		for i=0,nm-1 do  points(i)=   currency_option1(courbepd,courbepf,courbe_zeta,$
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
			; 								amotry
			fac=0.5
			fac1=(1.-fac)/float(ndim)
			fac2=fac1-fac
			for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
			;ytry=funk(ptry)						;appel a funk
			;#################     definition of funk(ptry)
			courbe_zeta(*,1)=ptry(1:ndim)
			courbe_zeta=courbe_exp(courbe_zeta)
			for i=0,nm-1 do  points(i)=  currency_option1(courbepd,courbepf,courbe_zeta,$
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
				  courbe_zeta(*,1)=ptry(1:ndim)
				  courbe_zeta=courbe_exp(courbe_zeta)
				  for i=0,nm-1 do  points(i)=  currency_option1(courbepd,courbepf,courbe_zeta,$
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
return,result
end


pro test_covariance
courbepd=transpose($
	[	[0,	1.],$
		[1.,	.95],$
		[2.,	.91],$
		[3.,	.864]])

courbepf=transpose($
	[	[0,	1.],$
		[1.,	.97],$
		[2.,	.85],$
		[3.,	.721]])
courbe_zeta=transpose($
	[	[0,	5.],$
		[1.,	4.],$
		[2.,	3.],$
		[3.,	2.]])

 
market_data=transpose([$
		[1,	100.,	100.,	1.,	 1.38297],$
		[1,	100.,	102.,	1.,	 0.618912],$
		[2,	100.,	100.,	1.,	 4.21440],$
		[1,	100.,	100.,	2.,	 1.30657],$
		[1,	100.,	102.,	2.,	 0.552642],$
		[2,	100.,	100.,	2.,	 16.2620],$
		[1,	100.,	100.,	3.,	 1.23729],$
		[1,	100.,	102.,	3.,	 0.490455],$
		[2,	100.,	100.,	3.,	 32.2448],$
		[2,	100.,	90.,	3.,	 19.0261]])

zeta0=2.
datelist=[0.,1.,2.,3.]
zeta=find_zeta(courbepd,courbepf,zeta0,datelist,market_data,30)
print,zeta
end

pro show_covariance
n=30
otype=1
t=3.
s=100.
k=100.
nm=20
courbepd=transpose($
	[	[0,	1.],$
		[1.,	.95],$
		[2.,	.91],$
		[3.,	.864]])

courbepf=transpose($
	[	[0,	1.],$
		[1.,	.97],$
		[2.,	.85],$
		[3.,	.721]])

courbe_zeta=transpose($
	[	[0,	4.],$
		[1.,	4.],$
		[2.,	4.],$
		[3.,	4.]])
courbe_zeta1=courbe_zeta
covarray=fltarr(nm)
valarray=fltarr(nm)
f=findgen(nm)/(nm-1.)+0.5
;on varie la premiere periode
for i=0,nm-1 do  begin
		courbe_zeta1(0)=courbe_zeta(0)*f(i)
		valarray(i)=  currency_option1(courbepd,courbepf,courbe_zeta1,$
 		otype,s,k,t,n)
	end

window,1,title='periode 1'
plot,courbe_zeta,valarray,xtitle='zeta volatilite periode 1',ytitle='value of the call'
end
