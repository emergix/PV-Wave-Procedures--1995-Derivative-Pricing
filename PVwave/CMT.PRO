 
function aj,d
;calcule nb d'annee a partir du jour jjmmaaaa ,
a=d mod 10000
m=((d-a) mod 1000000)/10000
j=(d-a-10000*m)/1000000
return,julday(m,j,a)/365.25
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

function courbetauxf
t0=aj(20101993)
v=[	[aj(21101993)-t0,	2.1683],$
	[aj(22101993)-t0,	2.1683],$
	[aj(29101993)-t0,	2.1676],$
	[aj(22111993)-t0,	2.2238],$
	[aj(15121993)-t0,	2.2888],$
	[aj(16031994)-t0,	2.3964],$
	[aj(15061994)-t0,	2.4217],$
	[aj(21091994)-t0,	2.4861],$
	[aj(21121994)-t0,	2.5624],$
	[aj(15031995)-t0,	2.6701],$
	[aj(21061995)-t0,	2.7706],$
	[aj(20091995)-t0,	2.8637],$
	[aj(20121995)-t0,	2.9561],$
	[aj(20031996)-t0,	3.0594],$
	[aj(19061996)-t0,	3.1472],$
	[aj(18091996)-t0,	3.2327],$
	[aj(18121996)-t0,	3.3151],$
	[aj(19031997)-t0,	3.4037],$
	[aj(18061997)-t0,	3.4802],$
	[aj(17091997)-t0,	3.5552],$
	[aj(17121997)-t0,	3.6267],$
	[aj(18031998)-t0,	3.7025],$
	[aj(17061998)-t0,	3.7691],$
	[aj(16091998)-t0,	3.8340],$
	[aj(16121998)-t0,	3.8951],$
	[aj(22041999)-t0,	3.9310],$
	[aj(22101999)-t0,	4.0030],$
	[aj(24042000)-t0,	4.0764],$
	[aj(23102000)-t0,	4.1512],$
	[aj(23042001)-t0,	4.2397],$
	[aj(22102001)-t0,	4.3289],$
	[aj(22042002)-t0,	4.4208],$
	[aj(22102002)-t0,	4.5130],$
	[aj(22042003)-t0,	4.6081],$
	[aj(22102003)-t0,	4.7039]]
return,v
 end



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


function couponszero,courb_taux
s=size(courb_taux)
n=s(1)
res=fltarr(n,2)
res(0,0)=0        ;t01
res(0,1)=courb_taux(0,1)  ;r01
for i=1,n-1 do begin
 res(i,0)=courb_taux(i-1,0)
  sum=courb_taux(0,1)
 for j=0,i-1 do begin
   if j eq 0 then begin
     somme=((1+courb_taux(i,1))^(courb_taux(j,0))-1)/$
	(1+res(j,1))^(courb_taux(j,0))
   endif else begin
     somme=somme+((1+courb_taux(i,1))^(courb_taux(j,0)$
	-courb_taux(j-1,0))-1)/(1+res(j,1))^(courb_taux(j,0))
   endelse 
 endfor
 res(i,1)=((1+courb_taux(i,1))^(courb_taux(i,0)-courb_taux(i-1,0))$
	/(1-somme))^(1/courb_taux(i,0))-1.
endfor
return,res
end;

function tauxforwardT,zerotaux,t,n
;n nombre de points desire pour la courbe forward
deb=min(zerotaux(*,0))
fin=max(zerotaux(*,0))
len=fin-deb
res=fltarr(n,2)
res(0,0)=0        ;t01
res(0,1)=interpolation(zerotaux,t)  ;r01
for i=1,n-1 do begin
	res(i,0)=i/float(n)*(len-t)
	z2=(1.+interpolation(zerotaux,res(i,0)+t))^(res(i,0)+t)
	z1=(1.+interpolation(zerotaux,res(i,0)))^res(i,0)
	z=z2/z1
	res(i,1)=z^(1/float(t))-1.
endfor
return,res
end;

function tauxswap,zerotaux,x,t,date_list
; x : dans x annees,
;T maturite du coupons zero 
; date_list dates de paiement du taux de swap equivalent
; resultat = taux equivalent verse aux date specifie par date_list
; n est pas annualise
 n=n_elements(date_list)
dt=1.-(((1.+interpolation(zerotaux,x))^x)/((1.+interpolation(zerotaux,(x+t)))^(x+t)))
s=0.
for i=0,n-1 do begin
	 s=s+(((1.+interpolation(zerotaux,x))^x)/((1.+interpolation(zerotaux,(x+date_list(i))))^(x+date_list(i))))
endfor
return,dt/s
end;


function tauxzero_forward,couponzeros
s=size(couponzeros)
n=s(1)
res=fltarr(n,2)
res(0,0)=0        ;t01
res(0,1)=couponzeros(0,1)  ;r01
for i=1,n-1 do begin
 res(i,0)=couponzeros(i,0)
 res(i,1)=exp(((alog(1+couponzeros(i,1))*couponzeros(i,0))-$
		(alog(1+couponzeros(i-1,1))*couponzeros(i-1,0)))/$
	(couponzeros(i,0)-couponzeros(i-1,0)))-1.
endfor
return,res
end;

function forward_tauxzero,forward
s=size(forward)
n=s(1)
res=fltarr(n,2)
res(0,0)=0        ;t01
res(0,1)=forward(0,1)  ;r01
for i=1,n-1 do begin
 res(i,0)=forward(i,0)
 res(i,1)= ((1.+res(i-1,1))^res(i-1,0)*(1.+forward(i,1))^(res(i,0)-res(i-1,0)))$
		^(1./res(i,0))-1.
endfor
return,res
end;

function zero_to_forex,crbzero1,crbzero2,s0
dates=vect_union(crbzero1(*,0),crbzero2(*,0))
nb=n_elements(dates)
change=fltarr(nb,2)
for i=0,nb-1 do begin
change(i,0)=dates(i)
change(i,1)=( (1.+interpolation(crbzero2,dates(i))) / (1.+interpolation(crbzero1,dates(i))) ) ^dates(i) *s0
endfor
return,change
end


; on convertit les taux monetaires en taux actuariels
function transformat_taux,courbe
scourbe=size(courbe)
nb=scourbe(1)
res=fltarr(nb,2)
for i=0,nb-1 do begin
	res(i,0)=courbe(i,0)
if res(i,0) lt 1. then $
	res(i,1)=(1.+courbe(i,1)/100.*courbe(i,0))^(1./courbe(i,0))-1. else $
	res(i,1)=courbe(i,1)/100.
endfor
return,res
end

  
;valeur d un swap CMT-side
 
function CMT_value,T,dt,neche,crbzero 
Date_list=(findgen(neche)+1.)/2.
pmt=fltarr(neche+1,4)
s=0.
 for i=1,neche do begin
	pmt(i,0)=dt*i
	pmt(i,1)=tauxswap(crbzero,pmt(i,0),T,Date_list)
	pmt(i,3)= (1.+interpolation(crbzero,pmt(i,0)))^pmt(i,0)
	s=s+(pmt(i,1)/pmt(i,3))
endfor

return,s
end

;valeur d un swap DIFF-side
 
function DIFF_value,T,dt,neche,crbzero ,crbzerof,ro
; ro n est pas utilise
 ;crbzerof est la courbe zero coupon etrangere
Date_list=(findgen(neche)+1.)/2.
pmt=fltarr(neche+1,4)
s=0.
s0=1.
crbchange=zero_to_forex(crbzero,crbzerof ,s0)
 for i=1,neche do begin
	pmt(i,0)=dt*i
	pmt(i,1)=tauxswap(crbzerof,pmt(i,0),T,Date_list) 
	pmt(i,2)=interpolation(crbchange,pmt(i,0))/s0
	pmt(i,3)= (1.+interpolation(crbzero,pmt(i,0)))^pmt(i,0)
	s=s+(pmt(i,1)*pmt(i,2)/pmt(i,3))
endfor

return,s
end

;valeur d un swap DIFF-side
 
function DIFF2_value,T,dt,neche,crbzero ,crbzerof,ro
;ro est la densite de correlation implicite 
; c est un vecteur de meme dimension que Date_list
;crbzerof est la courbe zero coupon etrangere
Date_list=(findgen(neche)+1.)/2.
pmt=fltarr(neche+1,4)
s=0.
 for i=1,neche do begin
	pmt(i,0)=dt*i
	pmt(i,1)=tauxswap(crbzerof,pmt(i,0),T,Date_list)-ro(i-1)*pmt(i,0)
	pmt(i,3)= (1.+interpolation(crbzero,pmt(i,0)))^pmt(i,0)
	s=s+(pmt(i,1)/pmt(i,3))
endfor

return,s
end


;valeur d un swap float-side
 
function Float_value,dt,neche,crbzero 
Date_list=(findgen(neche)+1.)/2.
pmt=fltarr(neche+1,4)
s=0.
for i=1,neche do begin
	pmt(i,0)=dt*i
	pmt(i,2)=tauxswap(crbzero,pmt(i,0),dt,[dt])
	pmt(i,3)= (1.+interpolation(crbzero,pmt(i,0)))^pmt(i,0)
	s=s +(pmt(i,2)/pmt(i,3))
endfor

return,s
end


 
 
;-----------------------------------------------------------------------------
;                            T E S T S -CMT
 


;bumps localise sur les taux zero coupons:cmt side
;swap pibor 6mois contre taux 5ans pendant 3ans
pro testb1cmt
T=5.
dt=0.5
neche=6
Date_list=[0.5,1.,1.5,2.,2.5,3.]
a=courbetaux()
scourbe=size(a)
nb=scourbe(2)
res=fltarr(nb,2)
b=transpose(a)
crbzero=transformat_taux(b)
 v0=CMT_value(T,dt,neche,crbzero)
s=0.
for i=0,nb-1 do begin
	res(i,0)=a(0,i)
crb=a
crb(1,i)=crb(1,i)+0.01
b=transpose(crb)
crbzero=transformat_taux(b)
 v=CMT_value(T,dt,neche,crbzero)-v0
	res(i,1)=v
	s=s+v
endfor
plot,res(*,0),res(*,1)*1000000.,psym=10,xtitle='years',ytitle='$/M$',title='CMS CMT side'
oplot,res(*,0),replicate(v0,nb),color=108
print,'sensibilite totale pour 1M$ : ',s*1000000.
set_plot,'PS'
device,/landscape
plot,res(*,0),res(*,1)*1000000.,psym=10,xtitle='years',ytitle='$/M$',title='CMS CMT side'
oplot,res(*,0),replicate(v0,nb),color=108
device,/close_file
spawn,'lpr wave.ps'
set_plot,'X'
end

 
 

;bumps localise sur les taux zero coupons:float side
;swap pibor 6mois contre taux 5ans pendant 3ans
pro testb1float
T=5.
dt=0.5
neche=6
Date_list=[0.5,1.,1.5,2.,2.5,3.]
a=courbetaux()
scourbe=size(a)
nb=scourbe(2)
res=fltarr(nb,2)
b=transpose(a)
crbzero=transformat_taux(b)
v0=float_value(dt,neche,crbzero)
s=0.
for i=0,nb-1 do begin
	res(i,0)=a(0,i)
crb=a
crb(1,i)=crb(1,i)+0.01
b=transpose(crb)
crbzero=transformat_taux(b)
v=float_value(dt,neche,crbzero)-v0
	res(i,1)=v
	s=s+v
endfor
plot,res(*,0),res(*,1)*1000000.,psym=10,title='CMS float side'
oplot,res(*,0),replicate(v0,nb),color=108
print,'sensibilite totale pour 1M$ : ',s*1000000.
set_plot,'PS'
device,/landscape
plot,res(*,0),res(*,1)*1000000.,psym=10,title='CMS float side'
oplot,res(*,0),replicate(v0,nb),color=108
device,/close_file
spawn,'lpr wave.ps'
set_plot,'X'


end

;bumps localise sur les taux zero coupons:total
;swap pibor 6mois contre taux 5ans pendant 3ans
pro testb1total
T=5.
dt=0.5
neche=6
Date_list=[0.5,1.,1.5,2.,2.5,3.]
a=courbetaux()
scourbe=size(a)
nb=scourbe(2)
res=fltarr(nb,2)
b=transpose(a)
crbzero=transformat_taux(b)
v0=-cmt_value(T,dt,neche,crbzero)+float_value(dt,neche,crbzero)
s=0.
for i=0,nb-1 do begin
	res(i,0)=a(0,i)
crb=a
crb(1,i)=crb(1,i)+0.01
b=transpose(crb)
crbzero=transformat_taux(b)
v=(-cmt_value(T,dt,neche,crbzero)+float_value(dt,neche,crbzero))-v0
	res(i,1)=v
	s=s+v
endfor
plot,res(*,0),res(*,1)*1000000.,psym=10,xtitle='years',ytitle='$/M$',title='CMS total sensitivity:pays CMT'
oplot,res(*,0),replicate(v0,nb),color=108
print,'sensibilite totale pour 1M$ : ',s*1000000.
set_plot,'PS'
device,/landscape
plot,res(*,0),res(*,1)*1000000.,psym=10,xtitle='years',ytitle='$/M$',title='CMS total sensitivity:pays CMT'
oplot,res(*,0),replicate(v0,nb),color=108
device,/close_file
spawn,'lpr wave.ps'
set_plot,'X'

end



;-------------------------------------------------------------------
;                 Test sur les DIFF


;bumps localise sur les taux zero coupons:diff side
;swap pibor 6mois americain contre taux 5ans pendant 3ans japonais
pro testd1diff
T=5.
dt=0.5
neche=6
ro0=0.001
ro=ro0*replicate(1,neche)
Date_list=[0.5,1.,1.5,2.,2.5,3.]
a=courbetaux()
af=courbetauxf()
scourbe=size(a)
scourbef=size(af)
nb=scourbe(2)
nbf=scourbef(2)
b=transpose(a)
bf=transpose(af)
crbzero=transformat_taux(b)
crbzerof=transformat_taux(bf)
 v0=DIFF_value(T,dt,neche,crbzero,crbzerof,ro)
ro=ro0*replicate(1.,neche)
s=0.
res=fltarr(nb,2)
for i=0,nb-1 do begin
	res(i,0)=a(0,i)
	crb=a
	crb(1,i)=crb(1,i)+0.01
	b=transpose(crb)
	crbzero=transformat_taux(b)
	 v=DIFF_value(T,dt,neche,crbzero,crbzerof,ro)-v0
	res(i,1)=v
	s=s+v
endfor
window,1
plot,res(*,0),res(*,1)*1000000.,psym=10,title='DIFF side :variation des taux US',xtitle='years',ytitle='$/M$'
oplot,res(*,0),replicate(v0,nb),color=108
set_plot,'PS'
device,/landscape
plot,res(*,0),res(*,1)*1000000.,psym=10,title='DIFF side :variation des taux US',xtitle='years',ytitle='$/M$'
oplot,res(*,0),replicate(v0,nb),color=108
device,/close_file
spawn,'lpr wave.ps'
set_plot,'X'
print,'sensibilite totale  varaition taux US pour 1M$ : ',s*1000000.
s=0.
res=fltarr(nbf,2)
for i=0,nb-1 do begin
	res(i,0)=af(0,i)
	crb=af
	crb(1,i)=crb(1,i)+0.01
	bf=transpose(crb)
	crbzerof=transformat_taux(bf)
 	v=DIFF_value(T,dt,neche,crbzero,crbzerof,ro)-v0
	res(i,1)=v
	s=s+v
endfor
window,2
plot,res(*,0),res(*,1)*1000000.,psym=10,title='DIFF side :variation des taux YEN ',xtitle='years',ytitle='$/M$'
oplot,res(*,0),replicate(v0,nb),color=108
set_plot,'PS'
device,/landscape
plot,res(*,0),res(*,1)*1000000.,psym=10,title='DIFF side :variation des taux YEN ',xtitle='years',ytitle='$/M$'
oplot,res(*,0),replicate(v0,nb),color=108
device,/close_file
spawn,'lpr wave.ps'
set_plot,'X'
print,'DIFF side :sensibilite totale  varaition taux YEN pour 1M$ : ',s*1000000.
print,'prix du swap :',v0*1000000.
b=transpose(a)
bf=transpose(af)
crbzero=transformat_taux(b)
crbzerof=transformat_taux(bf)
res=fltarr(neche+1)
roa=ro0*replicate(1.,neche)
for i=0,neche-1 do begin
	ro=roa
	ro(i)=ro(i)*0.9
	 	v=DIFF_value(T,dt,neche,crbzero,crbzerof,ro)-v0
		res(i+1)=v
endfor
ag=findgen(neche+1)/2.
window,3
plot,ag,res*1000000.,title='DIFF side :variation des densite de correlation : chute de 10% de ro',xtitle='years',ytitle='$/M$'
set_plot,'PS'
device,/landscape
plot,ag,res*1000000.,psym=10,title='DIFF side :variation des densite de correlation  ',xtitle='years',ytitle='$/M$'
device,/close_file
spawn,'lpr wave.ps'
set_plot,'X'

end
 



;-------------------------------------------------------------------
;               Test sur les courbes

pro test  
;calcule les taux zero et les taux forward a partir de taux de swap
b=transpose(courbetaux())
a=transformat_taux(b)
plot,a(*,0),a(*,1)
zero=couponszero(a)
forward=tauxforwardT(zero,1.,30)
plot,forward(*,0),forward(*,1),color=108
oplot,zero(*,0),zero(*,1),color=122
oplot,a(*,0),a(*,1)
end

 
pro testz
;bumps localise sur les taux zero coupons:float side
;swap pibor 6mois contre taux 5ans pendant 3ans
T=5.
dt=0.5
neche=6
Date_list=[0.5,1.,1.5,2.,2.5,3.]
a=courbetaux()
scourbe=size(a)
nb=scourbe(2)
res=fltarr(nb,2)
b=transpose(a)
crbzero=transformat_taux(b)
crbforwT=tauxforwardT(crbzero,T,100)
crbforw=tauxforwardT(crbzero,dt,100)
window,1
loadct,5
plot,crbforw(*,0),crbforw(*,1),color=179
oplot,crbzero(*,0),crbzero(*,1),color=91
 end
