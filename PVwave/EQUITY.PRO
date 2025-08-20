

;*************************calcul de la coube de volatilite implicite de lindex ***************


;******* calcul des points d'espacement egal en volatilite
function egal_volatilite_pts,vol_curve,n,t
delta_x=fltarr(n)
x=fltarr(n+1)
zsize=size(vol_curve)
nb=zsize(1)
vol2=fltarr(nb,2)
for i=0,nb-1 do begin
 vol2(i,0)=vol_curve(i,0)
 vol2(i,1)=vol_curve(i,1)^2
endfor
c2=integrale_courbec(vol2,0,t)
cn2=c2/float(n)
x(0)=0
for i=1,n do begin
	x(i)=integral_solve_y(vol2,x(i-1),cn2)
	delta_x(i-1)=x(i)-x(i-1)
endfor
structu={,delta_x:delta_x,x:x,v0:cn2}
return,structu
end

;calcul des taux forward a partir de la courbe zero coupon qui se comprend en interpolation continue
;pour trouver une valeur intermediare , donc il est important d avoir comme premier element
;de la courbe zero coupon le taux le plus court (over night)
function taux_forward,couponzeros
s=size(couponzeros)
n=s(1)
res=fltarr(n,2)
res(0,0)=0        ;t01
res(0,1)=couponzeros(0,1)  ;r01
for i=1,n-1 do begin
 res(i,0)=couponzeros(i,0)
 res(i,1)=exp(((alog(1+couponzeros(i,1))*couponzeros(i,0))-(alog(1+couponzeros(i-1,1))*couponzeros(i-1,0)))/$
	(couponzeros(i,0)-couponzeros(i-1,0)))-1.
endfor
return,res
end;


; pour pricer sur futur il faut remplacer s par f * exp(somme de forward rate  - somme des dividend yield entre 0 et T)
; et remplacer la valeur d exercice (s(i,j)-k) par (f(i,j)-k) ou 
;  f(i,j)=s(i,j)*exp(somme de forward rate  - somme des dividend yield entre i et T)
;        =s*u^i*d^(i-j)*exp(somme de forward rate  - somme des dividend yield entre 0 et T)
; ce qui est plus rapide a calculer que pour l option cash car pp ne depends pas de i
;   ***** a faire ****

;calcul d'un call
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
; les taux sont des taux zero coupons pour etre compatible avec HWLN , et les vol des vol forward.
; cashtype = 0 pour une option cash et 1 pour une option sur futur 
;   et s devient un prix futur d 'echeance l echeance de l option : t
function ACRYLoption,optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calule
;par les programmes courbetaux et courbevol
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vdeltat=replicate(t/float(n),n)
tvdeltat=findgen(n+1)*vdeltat(0)
vp=fltarr(n+1)
vtaux=replicate(r*t/float(n),n)
vvol=sigma^2*t/float(n)
if keyword_set(cvol) eq 1 then begin
	struu=egal_volatilite_pts(cvol,n,t)
	vdeltat=struu.delta_x
	tvdeltat=struu.x
	vvol=struu.v0
endif
if keyword_set(ctaux) eq 1 then begin
		ctaux1=taux_forward(ctaux)
		for i=0,n-2 do vtaux(i)=integrale_courbec(ctaux1,tvdeltat(i),tvdeltat(i+1))       ; interpolation de taux
endif
a=exp(vtaux)
u=exp(sqrt(vvol))
d=1/u
vp=(a-d)/(u-d)
indexlist=where(t_yieldlist lt t,count)
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if optype eq 0 then call1(i)=(supp(i)-k)>0 else call1(i)=(k-supp(i))>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if optype eq 0 then call1(i)=supp(i)*pp-k>0 else call1(i)=(k-supp(i)*pp)>0
	endfor
endelse
for j=0,n do call2(j)=0
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			for j=0,i do if optype eq 0 then call2(j)=(s*u^j*d^(i-j)-k)> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=(k-s*u^j*d^(i-j))> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			for j=0,i do if optype eq 0 then call2(j)=(s*u^j*d^(i-j)*pp-k)> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=(k-s*u^j*d^(i-j)*pp)> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end

;cacul de la volatility implicite d'une option
function index_impvol_spot,prime,optype,cashtype,amtype,s,k,t,ctaux,yieldlist,t_yieldlist,n=n,x1=x1,x2=x2,xacc=xacc,jmax=jmax
if prime eq 0 then return,0
if keyword_set(x1) eq 0 then x1=0.005
if keyword_set(x2) eq 0 then x2=0.3
if keyword_set(xacc) eq 0. then xacc=0.001
if keyword_set(jmax) eq 0. then jmax=40
if keyword_set(yieldlist) eq 0. then begin
	yieldlist=[0.] 
	t_yieldlist=[t/2.]
endif
r=interpolation(ctaux,t)
if amtype eq 0 then tta=t else tta=0.
if keyword_set(n) eq 0 then n=100
rmoyen=integrale_courbec(ctaux,0,t)/float(t)
volmin=rmoyen*sqrt(t/float(n))
x1=volmin*1.05
fmid= ACRYLoption(optype,cashtype,s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
f=ACRYLoption(optype,cashtype,s,k,t,tta,r,x1,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
if f ge 0 then print,'prime too small'
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
	fmid= ACRYLoption(optype,cashtype,s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
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
	fmid= ACRYLoption(optype,cashtype,s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end

;addition d une periode a une courbe de volatilite   
function form_vol_curve,impvol,datefin,x,t
nb=n_elements(impvol)
res=fltarr(nb+1,2)
res(0:nb-1,0)=datefin
res(nb,0)=t
res(0:nb-1,1)=impvol
res(nb,1)=x
return,res
end

;calcul d une volatilite forward a partir d un prix et d une courbe imcomplete de volatilite forward
function index_impvol_forward,prime,optype,cashtype,amtype,s,k,t,ctaux,yieldlist,t_yieldlist,impvol,datedebut,datefin,$
		n=n,x1=x1,x2=x2,xacc=xacc,jmax=jmax

if prime eq 0 then return,0
if keyword_set(x1) eq 0 then x1=0.005
if keyword_set(x2) eq 0 then x2=0.3
if keyword_set(xacc) eq 0. then xacc=0.001
if keyword_set(jmax) eq 0. then jmax=40
if keyword_set(yieldlist) eq 0. then begin
	yieldlist=[0.] 
	t_yieldlist=[t/2.]
endif
r=interpolation(ctaux,t)
if amtype eq 0 then tta=t else tta=0.
if keyword_set(n) eq 0 then n=100
rmoyen=integrale_courbec(ctaux,0,t)/float(t)
volmin=rmoyen*sqrt(t/float(n))
x1=volmin*1.05
cvol=form_vol_curve(impvol,datefin,x2,t)	
fmid= ACRYLoption(optype,cashtype,s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,$
	ctaux=ctaux,cvol=form_vol_curve(impvol,datefin,x2,t))-prime
f=ACRYLoption(optype,cashtype,s,k,t,tta,r,x1,n,yieldlist,t_yieldlist,$
	ctaux=ctaux,cvol=form_vol_curve(impvol,datefin,x1,t))-prime	
if f ge 0 then print,'prime too small'
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
	fmid= ACRYLoption(optype,cashtype,s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,$
			ctaux=ctaux,cvol=form_vol_curve(impvol,datefin,x2,t))-prime
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
	fmid= ACRYLoption(optype,cashtype,s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist,ctaux=ctaux,$
			cvol=form_vol_curve(impvol,datefin,xmid,t))-prime
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then begin
;	print,'arg=',optype,cashtype,s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist,ctaux,$
;			form_vol_curve(impvol,datefin,xmid,t),$
;	ACRYLoption(optype,cashtype,s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist,ctaux=ctaux,$
;			cvol=form_vol_curve(impvol,datefin,xmid,t)),prime
				return,rtbis
	end
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end

; calcul de toute la courbe de volatilite forward a partir d une liste de prix d option
function  index_imp_vol_curve,option_list,s,rate_curve,yieldlist,t_yieldlist,n
;option_list a le format suivant :
; 0 for call ,1 for put
; 0 for european , 1 for american
; strike
; maturity
; prix
; 0 pour cash ou 1 pour futur 
; on suppose la liste classee en fonction des maturite par ordre croissant
; avec une seule maturite a chaque fois
sizelist= size(option_list)
nb=sizelist(1)
impvol=fltarr(nb)
datedebut=fltarr(nb)
datefin=fltarr(nb)
datedebut(0)=0.
datefin(0)=option_list(0,3)
for i=0,nb-2 do begin
	datedebut(i+1)=option_list(i,3)
	datefin(i+1)=option_list(i+1,3)
endfor
impvol(0)=index_impvol_spot(option_list(0,4),option_list(0,0),option_list(0,5),$
		option_list(0,1),s,option_list(0,2),$
				option_list(0,3),rate_curve,yieldlist,t_yieldlist,n=n)
for i=1,nb-1 do begin
	impvol(i)=index_impvol_forward(option_list(i,4),option_list(i,0),option_list(0,5),$
		option_list(i,1),s,option_list(i,2),option_list(i,3),rate_curve,$
		yieldlist,t_yieldlist,impvol(0:i-1),datedebut(0:i-1),datefin(0:i-1),n=n)
		end
curve=fltarr(nb,2)
curve(*,0)=datefin
curve(*,1)=impvol
return,curve
end


;*****************************calcul des option barrieres***************
function upoutACRYLoption,optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calule
;par les programmes courbetaux et courbevol
;print,optype,s,k,t,tta,r,sigma,n
;print,yieldlist,t_yieldlist
;print,ctaux
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vp=fltarr(n+1)
vdeltat=replicate(t/float(n),n)
tvdeltat=findgen(n+1)*vdeltat(0)
vtaux=replicate(r*t/float(n),n)
vvol=sigma^2*t/float(n)
if keyword_set(cvol) eq 1 then begin
	struu=egal_volatilite_pts(cvol,n,t)
	vdeltat=struu.delta_x
	tvdeltat=struu.x
	vvol=struu.v0
endif
if keyword_set(ctaux) eq 1 then begin
		ctaux1=taux_forward(ctaux)
		for i=0,n-2 do vtaux(i)=integrale_courbec(ctaux1,tvdeltat(i),tvdeltat(i+1))       ; interpolation de taux
endif
a=exp(vtaux)
u=exp(sqrt(vvol))

d=1/u
vp=(a-d)/(u-d)
indexlist=where(t_yieldlist lt t,count)
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if supp(i) ge barriere then call1(i)=0 else if optype eq 0 then  call1(i)=(supp(i)-k)>0 else call1(i)=(k-supp(i))>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if supp(i)*pp ge barriere then call1(i)=0 else if optype eq 0 then call1(i)=(supp(i)*pp-k)>0 else call1(i)=(k-supp(i)*pp)>0
	endfor
endelse
for j=0,n do call2(j)=0
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			supp1=s*u^j*d^(i-j)
			for j=0,i do if supp1 ge barriere then call2(j)=0 else if optype eq 0 then call2(j)=(supp1-k)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=(k-supp1)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do if (s*u^j*d^(i-j) ge barriere) then call2(j)=0  else $
			 	call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			supp1=s*u^j*d^(i-j)*pp
			for j=0,i do if supp1 ge barriere then call2(j)=0 else if optype eq 0 then call2(j)=(supp1-k)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=(k-supp1)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do if s*u^j*d^(i-j)*pp ge barriere then call2(j)=0 else $
				call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end


function downoutACRYLoption,optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calule
;par les programmes courbetaux et courbevol
;print,optype,s,k,t,tta,r,sigma,n
;print,yieldlist,t_yieldlist
;print,ctaux
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vp=fltarr(n+1)
vdeltat=replicate(t/float(n),n)
tvdeltat=findgen(n+1)*vdeltat(0)
vtaux=replicate(r*t/float(n),n)
vvol=sigma^2*t/float(n)
if keyword_set(cvol) eq 1 then begin
	struu=egal_volatilite_pts(cvol,n,t)
	vdeltat=struu.delta_x
	tvdeltat=struu.x
	vvol=struu.v0
endif
if keyword_set(ctaux) eq 1 then begin
		ctaux1=taux_forward(ctaux)
		for i=0,n-2 do vtaux(i)=integrale_courbec(ctaux1,tvdeltat(i),tvdeltat(i+1))       ; interpolation de taux
endif
a=exp(vtaux)
u=exp(sqrt(vvol))
d=1/u
vp=(a-d)/(u-d)
indexlist=where(t_yieldlist lt t,count)
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if supp(i) le barriere then call1(i)=0 else if optype eq 0 then  call1(i)=(supp(i)-k)>0 else call1(i)=(k-supp(i))>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if supp(i)*pp le barriere then call1(i)=0 else if optype eq 0 then call1(i)=(supp(i)*pp-k)>0 else call1(i)=(k-supp(i)*pp)>0
	endfor
endelse
for j=0,n do call2(j)=0
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			supp1=s*u^j*d^(i-j)
			for j=0,i do if supp1 le barriere then call2(j)=0 else if optype eq 0 then call2(j)=(supp1-k)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=(k-supp1)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do if (s*u^j*d^(i-j) le barriere) then call2(j)=0  else $
			 	call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			supp1=s*u^j*d^(i-j)*pp
			for j=0,i do if supp1 le barriere then call2(j)=0 else if optype eq 0 then call2(j)=(supp1-k)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=(k-supp1)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do if s*u^j*d^(i-j)*pp le barriere then call2(j)=0 else $
				call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end

function downinACRYLoption,optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
return,ACRYLoption(optype,0,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)-$
	downoutACRYLoption(optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
end

function upinACRYLoption,optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
return,ACRYLoption(optype,0,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)-$
	upoutACRYLoption(optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
end

function BACRYLoption,optype,btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
case 1 of
 btype eq 0: return,ACRYLoption(optype,0,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
 btype eq 1: return,upinACRYLoption(optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
 btype eq 2: return,upoutACRYLoption(optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
 btype eq 3: return,downinACRYLoption(optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
 btype eq 4: return,downoutACRYLoption(optype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
else :Print,'BACRYLoption /  unknown type of option :',btype
endcase
end



;*************************calcul du yield equivalent*******************************

function equivalent_index_yield,div_list,s,nstep,zerocurve
;div_list est une liste de dividend de format ; (nb annee , dividende en point d'indice)
; s est la valeur de lindice en points
; nstep est le nombre de steps du vecteur resultat
; zerocurve est une coube zerocoupon
zsize=size(div_list)
ndiv=zsize(1)
lasttime=div_list(ndiv-1,0)
xcurve=fltarr(nstep,2)
rj=fltarr(ndiv)

for j=0,ndiv-1 do rj(j)=interpolation(zerocurve,div_list(j,0))


for i=1,nstep-1 do begin
	xcurve(i,0)=lasttime*i/nstep
	w=0
	ti=xcurve(i,0)
	di=(1+interpolation(zerocurve,ti))^ti
	for j=0,ndiv-1 do w=w+div_list(j,1)/(1+rj(j))^ti
	xcurve(i,1)=(1+(w*di/float(s)))^(1./ti)-1.
end
return,xcurve
end



;*********************** bucket analysis des options ******************************




;on se donne une option , une courbe de taux, une coube de vol 
; on rend un vecteur de derivee pour les taux et 
; un vecteur de derivee pour les vol
;  le grid pour ces vecteur est naturellement le grid de calcul binomial

function bucket_taux_option,optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumptaux=bumptaux
;creation d une courbe de taux equivalente composee de nb buckets
btaux=fltarr(nb,2)
delta_b=t/float(nb)
if keyword_set(ctaux) eq 0 then ctaux=transpose([[0.1,r],[10.,r]])
if keyword_set(bumptaux) eq 0 then bumptaux=0.005
for i=0,nb-1 do begin
	btaux(i,0)=(i+1)*delta_b
	btaux(i,1)=integrale_courbec(ctaux,i*delta_b,(i+1)*delta_b)
end
;calcul des buckets
opt_base=ACRYLoption(optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
bucket_rate=fltarr(nb,2)
for i=0,nb-1 do begin
	bucket_rate(i,0)=(i+1)*delta_b
	btaux(i,1)=btaux(i,1)+bumptaux
	bucket_rate(i,1)=ACRYLoption(optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=btaux,cvol=cvol)-opt_base
	btaux(i,1)=btaux(i,1)-bumptaux
end
return,bucket_rate
end

function bucket_taux_Boption,optype,Btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumptaux=bumptaux
;creation d une courbe de taux equivalente composee de nb buckets
btaux=fltarr(nb,2)
delta_b=t/float(nb)
if keyword_set(ctaux) eq 0 then ctaux=transpose([[0.1,r],[10.,r]])
if keyword_set(bumptaux) eq 0 then bumptaux=0.005
for i=0,nb-1 do begin
	btaux(i,0)=(i+1)*delta_b
	btaux(i,1)=integrale_courbec(ctaux,i*delta_b,(i+1)*delta_b)
end
;calcul des buckets
opt_base=BACRYLoption(optype,btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
bucket_rate=fltarr(nb,2)
for i=0,nb-1 do begin
	bucket_rate(i,0)=(i+1)*delta_b
	btaux(i,1)=btaux(i,1)+bumptaux
	bucket_rate(i,1)=BACRYLoption(optype,btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=btaux,cvol=cvol)-opt_base
	btaux(i,1)=btaux(i,1)-bumptaux
end
return,bucket_rate
end


function bucket_vol_option,optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumpvol=bumpvol
;creation d une courbe de vol equivalente composee de nb buckets
bvol=fltarr(nb,2)
delta_b=t/float(nb)
if keyword_set(cvol) eq 0 then cvol=transpose([[0.1,sigma],[10.,sigma]])
if keyword_set(bumpvol) eq 0 then bumpvol=0.01
for i=0,nb-1 do begin
	bvol(i,0)=(i+1)*delta_b
	bvol(i,1)=integrale_courbec(cvol,i*delta_b,(i+1)*delta_b)
end
;calcul des buckets
opt_base=ACRYLoption(optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
bucket_vol=fltarr(nb,2)
for i=0,nb-1 do begin
	bucket_vol(i,0)=(i+1)*delta_b
	bvol(i,1)=bvol(i,1)+bumpvol
	bucket_vol(i,1)=ACRYLoption(optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=bvol)-opt_base
	bvol(i,1)=bvol(i,1)-bumpvol
end
return,bucket_vol
end


function bucket_vol_Boption,optype,btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumpvol=bumpvol
;creation d une courbe de vol equivalente composee de nb buckets
bvol=fltarr(nb,2)
delta_b=t/float(nb)
if keyword_set(cvol) eq 0 then cvol=transpose([[0.1,sigma],[10.,sigma]])
if keyword_set(bumpvol) eq 0 then bumpvol=0.01
for i=0,nb-1 do begin
	bvol(i,0)=(i+1)*delta_b
	bvol(i,1)=integrale_courbec(cvol,i*delta_b,(i+1)*delta_b)
end
;calcul des buckets
opt_base=BACRYLoption(optype,btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
bucket_vol=fltarr(nb,2)
for i=0,nb-1 do begin
	bucket_vol(i,0)=(i+1)*delta_b
	bvol(i,1)=bvol(i,1)+bumpvol
	bucket_vol(i,1)=BACRYLoption(optype,btype,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=bvol)-opt_base
	bvol(i,1)=bvol(i,1)-bumpvol
end
return,bucket_vol
end





;********************* pricing d une option suivant son type 

function LNOption,r,sigma,n,yieldlist,t_yieldlist,s,portef,nbinstrument,ctaux=ctaux,cvol=cvol

ninst=nbinstrument
instrument_type=portef(1,0:ninst-1)
option_type=portef(2,0:ninst-1)
expiration_date=portef(3,0:ninst-1)
american_starting_date=portef(4,0:ninst-1)
sliding_flag=portef(5,0:ninst-1)
actual_price=portef(6,0:ninst-1)
maturity=portef(7,0:ninst-1)
contract_size=portef(8,0:ninst-1)
nb_of_contracts=portef(9,0:ninst-1)
american_flag=portef(10,0:ninst-1)
strike=portef(11,0:ninst-1)
underlying_type=portef(12,0:ninst-1)
barrier_type=portef(13,0:ninst-1)
barrier=portef(14,0:ninst-1)
notional=portef(15,0:ninst-1)
nb_payments=portef(16,0:ninst-1)
tenor=portef(17,0:ninst-1)
paymntdate=portef(20:39,0:ninst-1)
curveN=portef(40:59,0:ninst-1)
curvek=portef(60:79,0:ninst-1)
sum=0
for i=0,ninst-1 do begin
case instrument_type(i) of
	600: sum=sum+ACRYLoption(0,0,s,strike(i),expiration_date(i),expiration_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	601: sum=sum+ACRYLoption(0,1,s,strike(i),expiration_date(i),expiration_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	602: sum=sum+ACRYLoption(1,0,s,strike(i),expiration_date(i),expiration_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	604: sum=sum+ACRYLoption(1,1,s,strike(i),expiration_date(i),expiration_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	650: sum=sum+ACRYLoption(0,0,s,strike(i),expiration_date(i),0,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	651: sum=sum+ACRYLoption(0,1,s,strike(i),expiration_date(i),0,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	652: sum=sum+ACRYLoption(1,0,s,strike(i),expiration_date(i),0,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	654: sum=sum+ACRYLoption(1,1,s,strike(i),expiration_date(i),0,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	
	700: sum=sum+BACRYLoption(0,0,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	701: sum=sum+BACRYLoption(0,1,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	702: sum=sum+BACRYLoption(0,2,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	703: sum=sum+BACRYLoption(0,3,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	704: sum=sum+BACRYLoption(0,4,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	750: sum=sum+BACRYLoption(0,0,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	751: sum=sum+BACRYLoption(0,1,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	752: sum=sum+BACRYLoption(0,2,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	753: sum=sum+BACRYLoption(0,3,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	754: sum=sum+BACRYLoption(0,4,s,strike(i),barrier(i),expiration_date(i),american_starting_date(i),$
		r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
	else:print,'LNTOption: option type not recognized'
endcase
endfor
return,sum
end
	

function LNTOption,r,sigma,n,yieldlist,t_yieldlist,s,option,ctaux=ctaux,cvol=cvol
if n_elements(option) ne 80 then begin
	print,'error in the dimension of the option:',n_elements(option)
	return,0
	endif
portef=fltarr(80,1)
portef(0:79,0)=option
nbinstrument=1
resulta=LNOption(r,sigma,n,yieldlist,t_yieldlist,s,portef,nbinstrument,ctaux=ctaux,cvol=cvol)
return,resulta
end


; *********************bucket analysis dans le model LN *************


;on se donne une option , une courbe de taux, une coube de vol 
; on rend un vecteur de derivee pour les taux et 
; un vecteur de derivee pour les vol
;  le grid pour ces vecteur est naturellement le grid de calcul binomial

function bucket_taux_LNOption,r,sigma,n,yieldlist,t_yieldlist,s,portef,nbinstrument,nb,ctaux=ctaux,cvol=cvol,bumptaux=bumptaux
;creation d une courbe de taux equivalente composee de nb buckets
btaux=fltarr(nb,2)
t=max(ctaux(*,0))
delta_b=t/float(nb)
if keyword_set(bumptaux) eq 0 then bumptaux=0.005
for i=0,nb-1 do begin
	btaux(i,0)=(i+1)*delta_b
	btaux(i,1)=integrale_courbec(ctaux,i*delta_b,(i+1)*delta_b)
end
;calcul des buckets
opt_base=LNOption(r,sigma,n,yieldlist,t_yieldlist,s,portef,nbinstrument,ctaux=btaux,cvol=cvol)
bucket_rate=fltarr(nb,2)
for i=0,nb-1 do begin
	bucket_rate(i,0)=(i+1)*delta_b
	btaux(i,1)=btaux(i,1)+bumptaux
	bucket_rate(i,1)=LNOption(r,sigma,n,yieldlist,t_yieldlist,s,portef,nbinstrument,ctaux=btaux,cvol=cvol)-opt_base
	btaux(i,1)=btaux(i,1)-bumptaux
end
return,bucket_rate
end


;********************* pricing d un portefeuille d'option sur un arbre***************


; pour pricer sur futur il faut remplacer s par f * exp(somme de forward rate  - somme des dividend yield entre 0 et T)
; et remplacer la valeur d exercice (s(i,j)-k) par (f(i,j)-k) ou 
;  f(i,j)=s(i,j)*exp(somme de forward rate  - somme des dividend yield entre i et T)
;        =s*u^i*d^(i-j)*exp(somme de forward rate  - somme des dividend yield entre 0 et T)
; ce qui est plus rapide a calculer que pour l option cash car pp ne depends pas de i
;   ***** a faire ****

;calcul d'un call
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
; les taux sont des taux forward , et les vol des vol forward.
; cashtype = 0 pour une option cash et 1 pour une option sur futur 
;   et s devient un prix futur d 'echeance l echeance de l option : t
function ACRYLlistoption,optype,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calcule
;par les programmes courbetaux et courbevol
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vdeltat=replicate(t/float(n),n)
tvdeltat=findgen(n+1)*vdeltat(0)
vp=fltarr(n+1)
vtaux=replicate(r*t/float(n),n)
vvol=sigma^2*t/float(n)
if keyword_set(cvol) eq 1 then begin
	struu=egal_volatilite_pts(cvol,n,t)
	vdeltat=struu.delta_x
	tvdeltat=struu.x
	vvol=struu.v0
endif
if keyword_set(ctaux) eq 1 then begin
		for i=0,n-2 do vtaux(i)=integrale_courbec(ctaux,tvdeltat(i),tvdeltat(i+1))       ; interpolation de taux
endif
a=exp(vtaux)
u=exp(sqrt(vvol))
d=1/u
vp=(a-d)/(u-d)
indexlist=where(t_yieldlist lt t,count)
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if optype eq 0 then call1(i)=(supp(i)-k)>0 else call1(i)=(k-supp(i))>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		if optype eq 0 then call1(i)=supp(i)*pp-k>0 else call1(i)=(k-supp(i)*pp)>0
	endfor
endelse
for j=0,n do call2(j)=0
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			for j=0,i do if optype eq 0 then call2(j)=s*u^j*d^(i-j)-k> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=k-s*u^j*d^(i-j)> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			for j=0,i do if optype eq 0 then call2(j)=s*u^j*d^(i-j)*pp-k> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j)) else $
				call2(j)=k-s*u^j*d^(i-j)*pp> $
				exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end


pro testc01
s=2000.
k=2000.
t=12./12.
tta=t
r=0.1
sigma=0.2
n=100
cashtype=0
yieldlist=[0.0,0.0]
t_yieldlist=[.5,1.]
ctaux=transpose([[0.001,0.1],[0.5,0.1],[2.,0.1],[4.,0.057],[8.,0.06],[20.,0.063]])
cvol=transpose([[0.0833333, 0.2],[1.166667, 0.2]])
res=ACRYLoption(0,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist)
print,'prix=',res
;res=ACRYLoption(0,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
;print,'prix=',res
end

pro test07
ctaux=transpose([[0.001,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
yieldlist=[0.00,0.00]
t_yieldlist=[1.5,2.5]
optype=0
cashtype=0
s=100.
k=100.
t=5.
tta=0.
amtype=1
vol=0.19
n=50
cvol=transpose([[0.0001,vol],[10.,vol]])
prime=ACRYLoption(optype,cashtype,s,k,t,tta,0.1,vol,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix=',prime
print,'impvol=',index_impvol_spot(prime,optype,cashtype,amtype,s,k,t,ctaux,yieldlist,t_yieldlist,n=n)

end

pro test01
s=1000.
k=1000.
t=2
tta=2
r=0.1
sigma=0.2
n=50
cashtype=0
yieldlist=[0.01,0.01]
t_yieldlist=[.5,1.]
ctaux=transpose([[0.001,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
cvol=transpose([[0.001,0.19],[2.,0.22],[4.,0.21],[6.,0.2],[8.,0.2]])
res=ACRYLoption(0,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix=',res
end


; test de la volatilite forward
pro testf02
s=2000.
n=30
yieldlist=[0.01,0.01]
t_yieldlist=[.5,1.]
rate_curve=transpose([[0.001,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
option_list=transpose([	[0,	0,	2000.,	1/12.,	48.729,0],$
			[0,	0,	2000.,	2/12.,	69.9536,0],$
			[0,	0,	2000.,	3/12.,	86.4041,0],$
			[0,	0,	2000.,	6/12.,	125.986,0],$
			[0,	0,	2000.,	1.,	184.466,0],$
			[0,	0,	2000.,	2.,	289.06,0],$
			[0,	0,	2000.,	3.,	392.295,0]])

res=index_imp_vol_curve(option_list,s,rate_curve,yieldlist,t_yieldlist,n)
print,'[[',res(0,0),',',res(0,1),'],$'
for i=1,5 do print,'[',res(i,0),',',res(i,1),'],$'
print,'[',res(6,0),',',res(6,1),']]'
end

 ;prix pour les tests de la volatilite forward
pro testf01
s=2000.
n=30
yieldlist=[0.01,0.01]
t_yieldlist=[.5,1.]
vol_curve=transpose([[0.001, 0.19],[0.2, 0.18],[0.5, 0.17],[1., 0.18],[2., 0.19],[3. , 0.2]])
rate_curve=transpose([[0.001,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
option_list=transpose([	[0,	0,	2000.,	1/12.,	45.5739,0],$
			[0,	0,	2000.,	2/12.,	67.0479,0],$
			[0,	0,	2000.,	3/12.,	83.7679,0],$
			[0,	0,	2000.,	6/12.,	125.986,0],$
			[0,	0,	2000.,	1.,	185.324,0],$
			[0,	0,	2000.,	2.,	295.08,0],$
			[0,	0,	2000.,	3.,	402.725,0]])
res=fltarr(7)
for i=0,6 do res(i)=ACRYLoption(option_list(i,0),option_list(i,5),s,option_list(i,2),$
			option_list(i,3),0,0.0565,0.19056,n,yieldlist,t_yieldlist,ctaux=rate_curve,cvol=vol_curve)
print,res
end

pro test03b
s=2000.
n=30
yieldlist=[0.01,0.01]
t_yieldlist=[.5,1.]
vol_curve=transpose([[1/12., 0.19436],[2/12., 0.185526],[3/12., 0.177804],[6./12., 0.171773],[1., 0.177198],[2., 0.183041],[3. , 0.190560]])
rate_curve=transpose([[0.001,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
option_list=transpose([	[0,	0,	2000.,	1/12.,	45.5739,0],$
			[0,	0,	2000.,	2/12.,	67.0479,0],$
			[0,	0,	2000.,	3/12.,	83.7679,0],$
			[0,	0,	2000.,	6/12.,	125.986,0],$
			[0,	0,	2000.,	1.,	185.324,0],$
			[0,	0,	2000.,	2.,	295.08,0],$
			[0,	0,	2000.,	3.,	402.725,0]])
res=fltarr(7)
for i=6,6 do begin
	res(i)=ACRYLoption(option_list(i,0),option_list(i,5),s,option_list(i,2),$
			option_list(i,3),option_list(i,3),0.1,0.2,n,yieldlist,t_yieldlist,ctaux=rate_curve,cvol=vol_curve)
	print,'arg1=',option_list(i,0),option_list(i,5),s,option_list(i,2),$
			option_list(i,3),option_list(i,3),0.1,0.2,n,yieldlist,t_yieldlist,rate_curve,vol_curve,res(i)
end

print,res
end


pro test04
s=1000.
k=1000.
barriere=1250.
t=2
tta=2
r=0.1
sigma=0.2
n=50
yieldlist=[0.01,0.01]
t_yieldlist=[.5,1.]
ctaux=transpose([[0.1,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
cvol=transpose([[0.01,0.19],[2.,0.22],[4.,0.21],[6.,0.2],[8.,0.2]])
print,'prix  call=',BACRYLoption(0,0,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call upin=',BACRYLoption(0,1,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call upout=',BACRYLoption(0,2,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call downin=',BACRYLoption(0,3,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call downout=',BACRYLoption(0,4,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put=',BACRYLoption(1,0,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put upin=',BACRYLoption(1,1,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put upout=',BACRYLoption(1,2,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put downin=',BACRYLoption(1,3,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put downout=',BACRYLoption(1,4,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
barriere=850.
print,'prix  call=',BACRYLoption(0,0,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call upin=',BACRYLoption(0,1,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call upout=',BACRYLoption(0,2,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call downin=',BACRYLoption(0,3,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  call downout=',BACRYLoption(0,4,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put=',BACRYLoption(1,0,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put upin=',BACRYLoption(1,1,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put upout=',BACRYLoption(1,2,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put downin=',BACRYLoption(1,3,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix  put downout=',BACRYLoption(1,4,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)

end


pro test05
zerocurve=transpose([$
				[1/365.,	0.08375],$
				[1/12.,		0.080625],$
				[2/12.,		0.078125],$
				[3/12.,		0.0775],$
				[6/12.,		0.07375],$
				[9/12.,		0.07125],$
				[1.,		0.07125],$
				[2.,		0.0721],$
				[3.,		0.0723],$
				[4.,		0.0730],$
				[5.,		0.0737],$
				[7.,		0.0760],$
				[8.,		0.0768],$
				[9.,		0.0773],$
				[10.,		0.0780]])

s=2000.
nstep=10
div_list=transpose([$	[0.01,	12.],$
			[0.02,	3.],$
			[0.03,	5.],$
			[0.05,	3.],$
			[0.07,	3.],$
			[0.12,	4.],$
			[0.24,	10],$
			[0.35,	13.],$
			[0.55,	12.],$
			[0.66,	3.],$
			[0.88,	5.],$
			[0.97,	8.],$
			[0.98,	9.]$
				])
print,equivalent_index_yield(div_list,s,nstep,zerocurve)
end


pro test08
s=1000.
k=1000.
t=2
tta=2.
r=0.1
sigma=0.2
n=50
cashtype=0
yieldlist=[0.01,0.01,0.01]
t_yieldlist=[.5,1.,1.5]
ctaux=transpose([[0.001,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
cvol=transpose([[0.001,0.19],[2.,0.22],[4.,0.21],[6.,0.2],[8.,0.2]])
res=ACRYLoption(0,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix option=',res
bumptaux=0.005
nb=10
bb=bucket_taux_option(1,0,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumptaux=bumptaux)
vv=bucket_vol_option(1,0,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumpvol=bumpvol)
for i= 0,9 do begin
print,'date=',bb(i,0),'  bucket rate=',bb(i,1),'  bucket vol=',vv(i,1)
endfor
end






pro test08b
s=1000.
k=1000.
barriere=1300.
t=2
tta=1.
r=0.1
sigma=0.2
n=50
cashtype=0
yieldlist=[0.01,0.01,0.1]
t_yieldlist=[.5,1.,1.5]
ctaux=transpose([[0.001,0.04],[0.5,0.051],[2.,0.056],[4.,0.057],[8.,0.06],[20.,0.063]])
cvol=transpose([[0.001,0.19],[2.,0.22],[4.,0.21],[6.,0.2],[8.,0.2]])
res=ACRYLoption(0,cashtype,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
print,'prix option=',res
bumptaux=0.005
nb=10
for typeo=0,4 do begin
print,'type=',typeo
bb=bucket_taux_Boption(1,typeo,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumptaux=bumptaux)
vv=bucket_vol_Boption(1,typeo,s,k,barriere,t,tta,r,sigma,n,yieldlist,t_yieldlist,nb,ctaux=ctaux,cvol=cvol,bumpvol=bumpvol)
for i= 0,9 do begin
print,'date=',bb(i,0),'  bucket rate=',bb(i,1),'  bucket vol=',vv(i,1)
end
endfor
end




pro testb
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
cvol=transpose([[0.25,0.19],[1.,0.20],[3.,0.21],[5.,0.2],[7.,0.19]])
r=0.1
sigma=0.14
n=30
vol=0.2
s=100.
nbinstrument=5
portef=fltarr(80,nbinstrument)
portef(1:17,0)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,1)=[600,0,5.,0,0,0,0,200,-1,0,110.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,2)=[600,0,5.,0,0,0,0,200,-1,0,120.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,3)=[650,1,5.,0,0,0,0,200,5,0,80.,0,0,0,0,0,0]  ;un put 5 an
portef(1:17,4)=[650,1,5.,0,0,0,0,200,5,0,90.,0,0,0,0,0,0]  ;un put 5 an
yieldlist=[0.00,0.00]
t_yieldlist=[0.015,0.025]
nb=20
bucket=bucket_taux_LNOption(r,sigma,n,yieldlist,t_yieldlist,s,portef,nbinstrument,nb,ctaux=rz,cvol=cvol,bumptaux=bumptaux)
print,bucket
end
