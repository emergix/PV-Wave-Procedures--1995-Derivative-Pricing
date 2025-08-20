function int,x
return,x-(x mod double(1.))
end

;calcule le prix d une rente indexee par tg de valeur pc avec un taux sur le terme de taux
; n nombre d annees

;rente discrete versee a la fin d une periode de temps 
function rente_d,pc,n,tg,taux
s=0
for i=1,n do begin
	s=s+pc*((1+tg)/(1+taux))^i
	endfor
return,s
end

;rente continue
function rente_c,pc,n,tg,r
taux=alog((1+tg)/(1+r))
s=((taux+pc)*exp(taux*n)-pc)/taux
return,s
end

function plancheractuariel,coupon,r,n,nominal,prime,t_premiercoupon
;r est un taux discret
;est rembourse a l echeance nominal+prime
alpha=1/(1+r)
return,nominal*(coupon*alpha^t_premiercoupon*(alpha^n-1.)/(alpha-1.)+ $
			alpha^n*(nominal+prime)/nominal)
end

function soustrait_taux,taux,deltataux            ;soustrait un taux fixe a une courbe de taux
ns=size(taux)
n=ns(1)
res=fltarr(n,2)
for i=0,n-1 do begin
	res(i,0)=taux(i,0)
	res(i,1)=taux(i,1)-deltataux
endfor
return,res
end


function mult_div,s0,divlist            ;multiplie  une courbe de taux par un nombre
ns=size(divlist)
n=ns(1)
res=fltarr(n,2)
for i=0,n-1 do begin
	res(i,0)=divlist(i,0)
	res(i,1)=divlist(i,1)*s0
endfor
return,res
end

function div_time_shift,divlist,t,taux   
ns=size(divlist)			        ; calcule une prevision de date/yield pour t annee plus tard ns=size(divlist)
n=ns(1)                                         ;           avec un progression de taux au total
res=fltarr(n,2)
for i=0,n-1 do begin
	
	z=caldat(divlist(i,0))
	a=z-int(z/100.)*100
	j=int(z/10000.)
	mstar=z-a-j*10000
	a=a+1
	if a ge 100 then a = a -100
	ddate= a+mstar+j*10000.
	res(i,0)=jj(ddate)
	res(i,1)=divlist(i,1)*(1+taux)
endfor
return,res
end

function rendement_interne,rate_liste,x1=x1,x2=x2,xacc=xacc,jmax=jmax
;rate_liste est de la forme [[versement1,temps1],...,[versementn,tempsn]]
;retourne le taux qui egalise sigma(versement * 1/(1+taux)^(-temps)) = 0
div=rate_liste(0,*)
tt=rate_liste(1,*)
if keyword_set(x1) eq 0 then x1=0.02
if keyword_set(x2) eq 0 then x2=1.
if keyword_set(xacc) eq 0. then xacc=0.0001
if keyword_set(jmax) eq 0. then jmax=40
fmid=total(div*(1+x2)^(-tt))
f=total(div*(1+x1)^(-tt))
if f*fmid ge 0 then begin
	print,'x1 et x2 n encadre pas la racine : f[',x1,']=',f,' / f[',x2,']=',fmid
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
	fmid=total(div*(1+xmid)^(-tt))
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end

;****************************************************************************
;
;
;
;                  Tables de taux et de volatilites
;
;
;
;
;***************************************************************************


function courbetaux,marche			;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end

function marche_devise,marche			;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end

function devise_correlation,dev1,dev2		;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end

function courbevol,marche			;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end

function taux_emprunt_titre,marche		;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end


function prevision_yield,marche			;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end

function devise_hvolatilite,devise		;fonction qui sera redefinie dans operationnel.pro
print,'fichier operationnel non charge'
end


function dividende_yield,marche,s0
a=prevision_yield(marche,s0)
return,a(*,1)
end		
		
function dividende_date,marche,t0,s0
a=prevision_yield(marche,s0)
return,(a(*,0)-t0)/365.25
end	
		
function aujourdhui
return,jj(010170)+int(systime(1)/86400.)
end

function aa,nb
return,aujourdhui()+365.25*nb
end


;
 ;*******************************************************************************************************
;
;
;
;
;
;
;     Fonctions de calcul analytique classique
;
;            type Black and sholes
;
;
;
;
;*******************************************************************************************************

 function call,s,k,t,r,sigm
if (t le 0)  then return,0>(s-k) else begin
	if (sigm eq 0.0) then return,(s-k*exp(-r*t))>0 else begin
		d1=(alog(s/float(k))+(r+sigm*sigm/2)*t)/(sigm*sqrt(t))
		d2=d1-(sigm*sqrt(t))
		return,s*gaussint(d1)-k*gaussint(d2)*exp(-r*t)
	endelse
endelse
end



function put,s,k,t,r,sigm
if sigm eq 0.0 then return,k*exp(-r*t)-s>0
return,call(s,k,t,r,sigm)+k*exp(-r*t)-s>0
end

; pricing futur
; r normalement est nul
 function callF,F,k,t,r,sigm
if (t le 0)  then return,0>(F-k) else begin
	if (sigm eq 0.0) then return,(F-k)>0 else begin
		d1=(alog(F/float(k))+(sigm*sigm/2)*t)/(sigm*sqrt(t))
		d2=d1-(sigm*sqrt(t))
		return,(F*gaussint(d1)-k*gaussint(d2))*exp(-r*t)
	endelse
endelse
end


function putF,F,k,t,r,sigm
if sigm eq 0.0 then return,((k-F)*exp(-r*t))>0
return,callF(F,k,t,r,sigm)+(k-F)*exp(-r*t)>0
end



function  option,type,s,k,t,r,sigm
case  type of
1:return,call(s,k,t,r,sigm)
-1:return,-call(s,k,t,r,sigm)
2:return,put(s,k,t,r,sigm)
-2:return,-put(s,k,t,r,sigm)
3:return,s
-3:return,-s
else:return,0
endcase
end

;determination du taux de pret implicite
pro tx_pret_titre,c,p,k,t,r,spot,div
;c:prix du call
;p:prix du put
;k:prix d exercice
;t:echeance
;essayer de prendre des options un an
f=(c-p)*exp(r*t)+k
print,'le forward implicite est a (%) :',string(f)
bc=((spot*(1+r)^t)/f)-(div+1.)
print,'le cout du pret de titre implicite est (%):',string(100.*bc)
end



 
function future,s,t,r,div
f=s/((1+r-div)^t)
return,f
end

function divimplicit,f,s,t,r
div=1+r-(s/f)^(1/t)
return,div
end


function yield_dec,t1,t2,yieldlist,datediv_list
;retourne la depreciation passee du capital du aux differentes distributions
; successives de yield entre t1 et t2
indexlist=where((datediv_list le t2) and (datediv_list gt t1),count)
if count eq 0 then return,1.
ylist=yieldlist(indexlist)
p=1.
for i=0,count-1 do p=p*(1-ylist(i))
return,p
end

function cgamma,s,k,t,r,sigm
if t le 0 or sigm eq 0. then return,0.
res=double(0.5*sigm*sqrt(t)+alog(((double(1+r))^t*s)/k)/(sigm*sqrt(t)))
return,1/(double(2.50662827463)*sqrt(t)*sigm*exp(res^2/2)*s)
end

function  option_gamma,type,s,k,t,r,sigm
case  type of
1:return,cgamma(s,k,t,r,sigm)
-1:return,-cgamma(s,k,t,r,sigm)
2:return,cgamma(s,k,t,r,sigm)
-2:return,-cgamma(s,k,t,r,sigm)
3:return,0
-3:return,0
else:return,0
endcase
end


function gammamoyen,s,k,t,r,sigma,n
deltat=double(t)/n
a=exp(r*deltat)
u=exp(sigma*sqrt(deltat))
d=double(1/u)
p=(a-d)/(u-d)
gg=double(0.)
for n1=0,n-1 do begin
  for i=0,n1 do gg=gg+binomtwo(n1,i)*cgamma(s*u^i*d^(n1-i),k,deltat*(n-n1),r,sigma)	
endfor
return,gg/n
end

 
function ctheta,s,k,t,r,sigm
res=0.5*sigm*sqrt(t)+alog(((1+r)^t*s)/k)/(sigm*sqrt(t))
res=2.506628274631*sqrt(t)*exp(res^2/2)
return,-s*sigm/res-r*k*exp(-r*t)*gaussint((alog(s/k)+(r-(sigm*sigm)/2)*t)/$
							(sigm*sqrt(t)))
end
  

function cvega,s,k,t,r,sigm
a1=sigm*sqrt(t)
a2=k/(1+r)^t
a3=alog(s/a2)
 r1=a2*(sqrt(t)/2+a3/(sigm*a1))/exp((-a1/2+a3/a1)^2/2)
r2=s*(sqrt(t)/2-a3/(sigm*a1))/exp((a1/2+a3/a1)^2/2)
 return,(r1+r2)/2.506628274631
end

function crhocall,s,k,t,r,sigm
 d2=(alog(s/k)+(r-(sigm*sigm)/2)*t)/(sigm*sqrt(t))
return,k*t*gaussint(d2)*exp(-r*t)
end

function crhoput,s,k,t,r,sigm
 d2=(alog(s/k)+(r-(sigm*sigm)/2)*t)/(sigm*sqrt(t))
return,-k*t*gaussint(-d2)*exp(-r*t)
end

;calcul de la volatilite historique entre deux indices temporels d un historique normalise
function vol_historique,f,df,debut,fin
nb=fin-debut+1
ret=fltarr(nb)
for i=0,nb-1 do ret(i)=alog(f(4,i+1)/f(4,i))
retbar=avg(ret)
u=0.
for i=0,nb-1 do u=u+(ret(i)-retbar)^2
;vol=sqrt(u/(nb-1))*sqrt(365.*(fin-debut)/(df(fin)-df(debut)))
vol=sqrt(u/(nb-1))*sqrt(365.)
return,vol
end

;*******************************************************************************
;
;
;        declinaision des options analytiques classiques
;
;
;
;********************************************************************************

function Ycall,s1,k,t,r,sigm,yield,t_div
if t_div lt t then s=s1*(1-yield*exp(-r*t_div)) else s=s1
return,call(s,k,t,r,sigm)
end

function Yput,s1,k,t,r,sigm,yield,t_div
if t_div lt t then s=s1*(1-yield*exp(-r*t_div)) else s=s1
return,call(s,k,t,r,sigm)-s+k*exp(-r*t)
end
 

function Dcall,s1,k,t,r,sigm,div,t_div
if t_div lt t then s=s1-div*exp(-r*t_div) else s=s1
return,call(s,k,t,r,sigm)
end

function Dput,s1,k,t,r,sigm,div,t_div
if t_div lt t then s=s1-div*exp(-r*t_div) else s=s1
return,call(s,k,t,r,sigm)-s+k*exp(-r*t)
end

function YLcall,s1,k,t,r,sigm,yieldlist,t_divlist
indexlist=where(t_divlist lt t and t_divlist gt 0,count)
if count eq 0 then return,call(s1,k,t,r,sigm)
ylist=yieldlist(indexlist)
tlist=t_divlist(indexlist)
p=1.
for i=0,count-1 do p=p*(1-ylist(i)*exp(-r*(tlist(i)-tlist(0))))
return,Ycall(s1,k,t,r,sigm,(1.-p),t_divlist(indexlist(0)))
end

function YLcallF,s1,k,t,r,sigm,yieldlist,t_divlist
indexlist=where(t_divlist lt t and t_divlist gt 0,count)
if count eq 0 then return,callF(s1,k,t,r,sigm)
ylist=yieldlist(indexlist)
tlist=t_divlist(indexlist)
p=1.
for i=0,count-1 do p=p*(1-ylist(i)*exp(-r*(tlist(i)-tlist(0))))
return,YcallF(s1,k,t,r,sigm,(1.-p),t_divlist(indexlist(0)))
end


function YLput,s1,k,t,r,sigm,yieldlist,t_divlist
indexlist=where(t_divlist lt t,count)
if count eq 0 then return,put(s1,k,t,r,sigm)
ylist=yieldlist(indexlist)
tlist=t_divlist(indexlist)
p=1.
for i=0,count-1 do p=p*(1-ylist(i)*exp(-r*(tlist(i)-tlist(0))))
return,Yput(s1,k,t,r,sigm,1.-p,t_divlist(indexlist(0)))
end


function YLputF,s1,k,t,r,sigm,yieldlist,t_divlist
indexlist=where(t_divlist lt t,count)
if count eq 0 then return,put(s1,k,t,r,sigm)
ylist=yieldlist(indexlist)
tlist=t_divlist(indexlist)
p=1.
for i=0,count-1 do p=p*(1-ylist(i)*exp(-r*(tlist(i)-tlist(0))))
return,YputF(s1,k,t,r,sigm,1.-p,t_divlist(indexlist(0)))
end


 
function DLcall,s1,k,t,r,sigm,divlist,t_divlist
indexlist=where(t_divlist lt t,count)
if count eq 0 then return,call(s1,k,t,r,sigm)
dlist=divlist(indexlist)
tlist=t_divlist(indexlist)
p=0.
for i=0,count-1 do p=p+dlist(i)*exp(-r*(tlist(i)-tlist(0)))
return,Dcall(s1,k,t,r,sigm,p,t_divlist(indexlist(0)))
end

function DLput,s1,k,t,r,sigm,divlist,t_divlist
indexlist=where(t_divlist lt t,count)
if count eq 0 then return,put(s1,k,t,r,sigm)
dlist=divlist(indexlist)
tlist=t_divlist(indexlist)
p=0.
for i=0,count-1 do p=p+dlist(i)*exp(-r*(tlist(i)-tlist(0)))
return,Dput(s1,k,t,r,sigm,p,t_divlist(indexlist(0)))
end


;******************************************************************************
;
;
;formule de valorisation d'un call rappelable au niveau C
;
;
;******************************************************************************


function RAPcall,s,k,c,t,r,sigma
st=(sigma*sqrt(t))
alpha=r+sigma^2/2
beta=r-sigma^2/2
u1=(s/c)*gaussint((alog(s/c)+alpha*t)/st)+(s/c)^(-2*r^2/sigma^2)*gaussint((alog(s/c)-alpha*t)/st)
v1c=gaussint((alog(c/s)-alpha*t)/st)-(s/c)^(-(1+2*r/sigma^2))*gaussint((alog(s/c)-alpha*t)/st)
v2c=gaussint((alog(c/s)-beta*t)/st)-(s/c)^(1-2*r/sigma^2)*gaussint((alog(s/c)-beta*t)/st)
v1k=gaussint((alog(k/s)-alpha*t)/st)-(s/c)^(-(1+2*r/sigma^2))*gaussint((alog(k*s/c^2)-alpha*t)/st)
v2k=gaussint((alog(k/s)-beta*t)/st)-(s/c)^(1-2*r/sigma^2)*gaussint((alog(k*s/c^2)-beta*t)/st)
u=(c-k)*u1
vc=s*v1c-k*exp(-r*t)*v2c
vk=s*v1k-k*exp(-r*t)*v2k
return,u+vc-vk
end

function RAPcall_delta,s,k,c,t,r,sigma
  s1=s*1.005
  s2=s*0.995
  c1=RAPcall(s1,k,c,t,r,sigma)
  c2=RAPcall(s2,k,c,t,r,sigma)
return,100*(c1-c2)/s
end

function RAPput,s,k,t,r,sigm
if sigm eq 0.0 then return,k*exp(-r*t)-s>0
return,RAPcall(s,k,c,t,r,sigm)+k*exp(-r*t)-s>0
end

;******************************************************************************
;
;
;
;         options lookback
;
;
;
;
;******************************************************************************

;valeur d un call permettant d acheter un stock a une valeur minimum q realise
;qui est le minimum atteint par s depuis le debut de l option
;c est une option dependant du chemin donc il faut reajuster constament q
function CMINcall,s,q,t,r,sigma
b=alog(s/q)
mu=r-sigma^2/2
st=sigma*sqrt(t)
vmin1=q*exp(-r*t)*(gaussint((b+mu*t)/st)-sigma^2/(2*r)* $
	exp(-2*mu*b/sigma^2)*gaussint((-b+mu*t)/st))
vmin2=s*(1+sigma^2/(2*r))*gaussint((-b-(mu+sigma^2)*t)/st)
return,s-(vmin1+vmin2)
end

;valeur d un call permettant de vendre un stock a une valeur maximum m realise
;qui est le maximum atteint par s depuis le debut de l option
;c est une option dependant du chemin donc il faut reajuster constament m
function CMAXput,s,m,t,r,sigma
a=alog(m/s)
mu=r-sigma^2/2
st=sigma*sqrt(t)
vmax1=m*exp(-r*t)*(gaussint((a-mu*t)/st)-sigma^2/(2*r)* $
	exp(2*mu*a/sigma^2)*gaussint((-a-mu*t)/st))
vmax2=s*(1+sigma^2/(2*r))*(1-gaussint((a-(mu+sigma^2)*t)/st))
return,(vmax1+vmax2)-s
end

function CMINcall_delta,s,q,t,r,sigma
  s1=s*1.005
  s2=s*0.995
  c1=CMINcall(s1,q,t,r,sigma)
  c2=CMINcall(s2,q,t,r,sigma)
return,100*(c1-c2)/s
end

function CMAXput_delta,s,m,t,r,sigma
  s1=s*1.005
  s2=s*0.995
  c1=CMAXput(s1,m,t,r,sigma)
  c2=CMAXput(s2,m,t,r,sigma)
return,100*(c1-c2)/s
end

;******************************************************************************
;
;
;
;          formule de calculs de option impliquants des resets europeens
;
;
;
;******************************************************************************


;calcul d un call at the money forwade de td
function OYLcall,td,s,k,t,r,sigma,n,yieldlist,t_divlist
sup=exp((findgen(2*n+1)-n)/float(n)*3*sigma*sqrt(td))*s
vec=fltarr(2*n+1)
vec1=fltarr(2*n+1)
for i=0,2*n do begin
	vec(i)=YLcall(sup(i),sup(i),t,r,sigma,yieldlist,t_divlist)*exp(-((i-n)/float(n))^2)
	vec1(i)=exp(-((i-n)/float(n))^2)
endfor
return,total(vec)/total(vec1)*exp(-r*td)
end

;calcul d un put forwarde de td
function OYLput,td,s,k,t,r,sigma,n,yieldlist,t_divlist
sup=exp((findgen(2*n+1)-n)/float(n)*3*sigma*sqrt(td))*s
vec=fltarr(2*n+1)
vec1=fltarr(2*n+1)
for i=0,2*n do begin
	vec(i)=YLput(sup(i),sup(i),t,r,sigma,yieldlist,t_divlist)*exp(-((i-n)/float(n))^2)
	vec1(i)=exp(-((i-n)/float(n))^2)
endfor
return,total(vec)/total(vec1)*exp(-r*td)
end


function NYLRScall,s01,k01,t01,taux,sig,n,yield_list,datediv_list,datereset_list1
;calcul d un call comportant une liste de resets de type europeens
datereset_list=datereset_list1(sort(datereset_list1))
s=float(s01)
k=float(k01)
t=float(t01)
nr=n_elements(datereset_list)
if (nr eq 0) then return,YLcall(s,k,t,r,sig,yield_list,datediv_list)
sum=0
for ri=0,nr do begin
	if ri eq 0 then begin
		td=0
		tf=datereset_list(ri)
	endif
	if ((ri gt 0) and (ri lt nr-1)) then begin
		td=datereset_list(ri-1)
		tf=datereset_list(ri)
	endif
	if ri eq nr then begin
		td=datereset_list(ri-1)
		tf=t
	endif
sum=sum+OYLcall(td,s,k,tf-td,taux,sig,n,yield_list,datediv_list-td)
endfor
return,sum
end

function NYLRSput,s01,k01,t01,taux,sig,n,yield_list,datediv_list,datereset_list1
;calcul d un put impliquant une liste de reset de type europeens
datereset_list=datereset_list1(sort(datereset_list1))
s=float(s01)
k=float(k01)
t=float(t01)
nr=n_elements(datereset_list)
if (nr eq 0) then return,YLput(s,k,t,r,sig,yield_list,datediv_list)
sum=0
for ri=0,nr do begin
	if ri eq 0 then begin
		td=0
		tf=datereset_list(ri)
	endif
	if ((ri gt 0) and (ri lt nr-1)) then begin
		td=datereset_list(ri-1)
		tf=datereset_list(ri)
	endif
	if ri eq nr then begin
		td=datereset_list(ri-1)
		tf=t
	endif
sum=sum+OYLput(td,s,k,tf-td,taux,sig,n,yield_list,datediv_list-td)
endfor
return,sum
end
;*****************************************************************************************
;
;
;
;                Formule de calcul diffusion-jump (cox p. 370)
;
;
;
;*****************************************************************************************

function jmpcall,s,k,t,r,v,l,d
; l est le parametre de saut 
; d est la variance de la distributuion du parametre de saut
term0=s/10000.
term=term0
som=0
i=0
while abs(term) ge term0 do begin
	term=exp(-l*t)*(l*t)^i*call(s,k,t,r,sqrt(v^2+d^2*i/t))
	som=som+term
end
return,som
end


;*****************************************************************************************
;
;
;
;
;      Fonction de calcul optionnel d 'inspiration Cox Ross
;
;
;
;
;
;
;****************************************************************************************
 
 
function CRcall,s,k,t,r,sigma,n
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
 for i=0,n do begin
	       supp(i)=s*(u^i)*(d^(n-i))
	       call1(i)=supp(i)-k>0
	      endfor
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	for j=0,i do begin
	call2(j)=s*u^j*d^(i-j)-k> $
		exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
	endfor
	for j=0,i do call1(j)=call2(j)
	endfor
return,call2(0)
end


 
function CRput,s,k,t,r,sigma,n
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([k-s,0])
a=exp(r*deltat) 
put1=fltarr(n+1)
put2=fltarr(n+1)
supp=fltarr(n+1)
u=exp(sigma*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
 for i=0,n do begin
	supp(i)=s*(u^i)*(d^(n-i))
	put1(i)=k-supp(i)>0
	endfor
for i=n-1,0,-1 do begin
 	for j=0,n do  put2(j)=0
	for j=0,i do put2(j)=k-s*u^j*d^(i-j)> $
		exp(-r*deltat)*(p*put1(j+1)+(1-p)*put1(j))
	for j=0,i do put1(j)=put2(j)
	endfor
return,put2(0)
end
 
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


function CRYLput,s,k,t,r,sigma,n,yieldlist,t_divlist
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([k-s,0])
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
		call1(i)=k-supp(i)>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=k-supp(i)*pp>0
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		for j=0,i do call2(j)=k-s*u^j*d^(i-j)> $
	exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		for j=0,i do call2(j)=k-s*u^j*d^(i-j)*pp> $
	exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end


;calcul d'un call
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
function ACRYLcall,s,k,t,tta,r,sigma,n,yieldlist,t_divlist
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
if i*deltat gt tta then begin
		for j=0,i do call2(j)=s*u^j*d^(i-j)-k> $
	exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
endif else begin
	for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
if i*deltat gt tta then begin
		for j=0,i do call2(j)=s*u^j*d^(i-j)*pp-k> $
	exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
endif else begin
	for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end


;calcul d'un put
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
function ACRYLput,s,k,t,tta,r,sigma,n,yieldlist,t_divlist
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([k-s,0])
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
		call1(i)=k-supp(i)>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=k-supp(i)*pp>0
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		if i*deltat gt tta then begin
			for j=0,i do call2(j)=k-s*u^j*d^(i-j)> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if i*deltat gt tta then begin
			for j=0,i do call2(j)=k-s*u^j*d^(i-j)*pp> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endelse
	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end

;calcul d'un call
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
function ACRYLcall1,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calule
;par les programmes courbetaux et courbevol
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vdeltat=fltarr(n+1)
vp=fltarr(n+1)
vtaux=replicate(r,n)
vvol=replicate(sigma,n)
if keyword_set(ctaux) ne 0 then begin
for i=0,n-1 do vtaux(i)=interpolation(ctaux,deltat*i)/100.
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
indexlist=where(t_yieldlist lt t,count)
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
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			for j=0,i do call2(j)=s*u^j*d^(i-j)-k> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			for j=0,i do call2(j)=s*u^j*d^(i-j)*pp-k> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end

;calcul d'un call  vrai formule
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
function ACRYLcall2,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calule
;par les programmes courbetaux et courbevol
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vp=fltarr(n+1)
vtaux=replicate(r,n)
vvol=replicate(sigma,n)
vdeltat=replicate(t/float(n),n)
if keyword_set(cvol) ne 0 then begin
	CC=integrale_carre_courbe(cvol,t)/float(n)
	ti=0.
	for i=0,n-1 do begin
		vdeltat(i)=cc/interpolation(cvol,ti)^2
		ti=ti+vdeltat(i)
	endfor
	kk=t/ti
	vdeltat=vdeltat*kk/ti
	for i=0,n-1 do begin
		vtaux(i)=cc/interpolation(ctaux,ti)^2
		ti=ti+vdeltat(i)
	endfor
endif
 a=exp(vtaux*vdeltat)
; ici il faut avoir calcule vdeltat(i) , vtaux(i) et vvol(i) ainsi que u 
a=exp(vtaux*vdeltat)
d=1/u
vp=(a-d)/(u-d)
indexlist=where(t_yieldlist lt t,count)
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
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			for j=0,i do call2(j)=s*u^j*d^(i-j)-k> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			for j=0,i do call2(j)=s*u^j*d^(i-j)*pp-k> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end




;calcul d'un put
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
function ACRYLput1,s,k,t,tta,r,sigma,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol
;ctaux et cvol sont les courbes de taux jj et de volatilite implicite courte forward tel que calule
;par les programmes courbetaux et courbevol
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([k-s,0])
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
vdeltat=fltarr(n+1)
vp=fltarr(n+1)
vtaux=replicate(r,n)
vvol=replicate(sigma,n)
if keyword_set(ctaux) ne 0 then begin
for i=0,n-1 do vtaux(i)=interpolation(ctaux,deltat*i)/100.
endif
if keyword_set(cvol) ne 0 then begin
for i=0,n-1 do vvol(i)=interpolation(cvol,deltat*i)/100.
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
indexlist=where(t_yieldlist lt t,count)
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=k-supp(i)>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=k-supp(i)*pp>0
	endfor
endelse
for j=0,n do call2(j)=0
tcalcul=t
for i=n-1,0,-1 do begin
	tcalcul=tcalcul-vdeltat(i)
	indexlist=where(t_yieldlist lt tcalcul,count)
	if count eq 0 then  begin 
		if tcalcul gt tta then begin
			for j=0,i do call2(j)=k-s*u^j*d^(i-j)> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if tcalcul gt tta then begin
			for j=0,i do call2(j)=k-s*u^j*d^(i-j)*pp> $
				exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-vtaux(i)*vdeltat(i))*(vp(i)*call1(j+1)+(1-vp(i))*call1(j))
		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end



 ;*****************************************************************************************
;
;
;
;
;      calcul d option avec palier (de type explosif)
;
;
;
;
;
;
;****************************************************************************************
 


function CREcall,s,k,t,r,sigma,n,palier
;calcul d un call comportant un palier de type explosif
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
 for i=0,n do begin
	       supp(i)=s*(u^i)*(d^(n-i))
	       call1(i)=(palier-k)<((supp(i)-k)>0)
	      endfor
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	for j=0,i do call2(j)=palier<(s*u^j*d^(i-j)-k> $
		(exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))))
	for j=0,i do call1(j)=call2(j)
	endfor
return,call2(0)
end


function CREput,s,k,t,r,sigma,n,palier
;calcul d un put de type explosif a un palier
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([k-s,0])
a=exp(r*deltat) 
put1=fltarr(n+1)
put2=fltarr(n+1)
supp=fltarr(n+1)
u=exp(sigma*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
 for i=0,n do begin
	supp(i)=s*(u^i)*(d^(n-i))
	put1(i)=(k-palier)<((k-supp(i))>0)
	endfor
for i=n-1,0,-1 do begin
 	for j=0,n do  put2(j)=0
	for j=0,i do put2(j)=palier<((k-s*u^j*d^(i-j))> $
		exp(-r*deltat)*(p*put1(j+1)+(1-p)*put1(j)))
	for j=0,i do put1(j)=put2(j)
	endfor
return,put2(0)
end

function ACRYLEcall,s,k,t,tta,ttl,r,sigma,n,yieldlist,t_divlist,palier
;plafond d explosion du call : au moins superieur a s
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
ttla=ttl < tta
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
		call1(i)=(palier-k)<((supp(i)-k)>0)
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=(palier-k)<((supp(i)*pp-k)>0)
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		if i*deltat gt ttla then begin
			if i*deltat lt tta then zportage=exp(-r*(tta-i*deltat)) else zportage=1.
			for j=0,i do call2(j)=((palier-k)<(s*u^j*d^(i-j)-k))*zportage> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if i*deltat gt ttla then begin
			if i*deltat lt tta then zportage=exp(-r*(tta-i*deltat)) else zportage=1. 
			for j=0,i do call2(j)=((palier-k)<(s*u^j*d^(i-j)*pp-k))*zportage> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endif else begin

			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endelse
	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end
 
function ACRYLEput,s,k,t,tta,ttl,r,sigma,n,yieldlist,t_divlist,palier
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([k-s,0])
ttla=ttl < tta
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
		call1(i)=(k-palier)<(k-supp(i)>0)
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=(k-palier)<(k-supp(i)*pp>0)
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		if i*deltat gt ttla then begin
			if i*deltat lt tta then zportage=exp(-r*(tta-i*deltat)) else zportage=1.
			for j=0,i do call2(j)=((k-palier)<(k-s*u^j*d^(i-j)))*zportage> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if i*deltat gt ttla then begin
			if i*deltat lt tta then zportage=exp(-r*(tta-i*deltat)) else zportage=1.
			for j=0,i do call2(j)=((k-palier)<(k-s*u^j*d^(i-j)*pp))*zportage> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
		endelse
	endelse
	for j=0,i do call1(j)=call2(j)
endfor
return,call2(0)
end


;*********************************************************************************************
;
;
;
;  formule de calculs d options impliquants des paliers  et presentant un caratere americain
;
;
;
;
;*********************************************************************************************

; americaine apartir de tta 
; paliers a partir de ttl
function NACRYLEcall,s01,k01,t01,tta,ttl,taux,sig,n,yield_list,datediv_list,palier_list1
palier_list2=palier_list1(sort(palier_list1))
s=float(s01)
k=float(k01)
t=float(t01)
indexpalier=where(palier_list2 > k,nr)
if (nr eq 0) then return,ACRYLcall(s,k,t,tta,r,sig,n,yield_list,datediv_list)
palier_list=palier_list2(indexpalier)
sum=0
for ri=0,nr-1 do begin
	if ri eq 0 then begin
		kk=k
		palier=palier_list(ri)
	endif else  begin
		kk=palier_list(ri-1)
		palier=palier_list(ri)
	endelse
sum1=ACRYLEcall(s,kk,t,tta,ttl,taux,sig,n,yield_list,datediv_list,palier)
sum=sum+sum1
endfor
sum1=ACRYLcall(s,palier_list(nr-1),t,tta,taux,sig,n,yield_list,datediv_list)
sum=sum+sum1
return,sum
end

function NACRYLEput,s01,k01,t01,tta,ttl,taux,sig,n,yield_list,datediv_list,palier_list1
palier_list2=palier_list1(sort(palier_list1))
indexpalier=where(palier_list2 < k,nr)
s=float(s01)
k=float(k01)
t=float(t01)
indexpalier=where(palier_list2 < k,nr)
if (nr eq 0) then return,ACRYLput(s,k,t,tta,r,sig,n,yield_list,datediv_list)
palier_list=palier_list2(indexpalier)
sum=0
for ri=1,nr do begin
	if ri eq nr then begin
		kk=k
		palier=palier_list(nr-1)
	endif else  begin
		palier=palier_list(ri-1)
		kk=palier_list(ri)
	endelse
sum=sum+ACRYLEput(s,kk,t,tta,ttl,taux,sig,n,yield_list,datediv_list,palier)
endfor
sum=sum+ACRYLput(s,palier_list(0),t,tta,taux,sig,n,yield_list,datediv_list)
return,sum
end

;*********************************************************************************************
;
;
;
;  formule de calculs d options spectrales (levrage controle par une fonction)
;            et presentant un caratere americain
;
;
;
;
;*********************************************************************************************
; americaine apartir de tta 
; paliers a partir de ttl
; les paliers supposes croissants sont de la forme [a1,a2,..., an]
;  les coefficients [b1,b2,...,bn+1]
; concernent la tranche : b1  -> [s01,a1]
;			  b2  -> [a1,a2]
;			 ...
;   			  bn  -> [an-1,an]
;			  bn+1-> [an, infini]
function SPACRYLEcall,s01,k01,t01,tta,ttl,taux,sig,n,yield_list,datediv_list,palier_list2,coeff_list2
n2=n_elements(palier_list2)
n2p=n_elements(coeff_list2)
if n2p ne n2+1 then print,'erreur dans le nombre de coefficients fournit:',n2p,'#',n2,'+1'
s=float(s01)
k=float(k01)
t=float(t01)
indexpalier=where(palier_list2 > k,nr)
if (nr eq 0) then return,ACRYLcall(s,k,t,tta,r,sig,n,yield_list,datediv_list)
palier_list=palier_list2(indexpalier)
coeff_list=coeff_list2(indexpalier)
sum=0
for ri=0,nr-1 do begin
	if ri eq 0 then begin
		kk=k
		palier=palier_list(ri)
	endif else  begin
		kk=palier_list(ri-1)
		palier=palier_list(ri)
	endelse
sum1=ACRYLEcall(s,kk,t,tta,ttl,taux,sig,n,yield_list,datediv_list,palier)
sum=sum+sum1*coeff_list(ri)
endfor
sum1=ACRYLcall(s,palier_list(nr-1),t,tta,taux,sig,n,yield_list,datediv_list)
sum=sum+sum1*coeff_list2(nr)
return,sum
end

function SPACRYLEput,s01,k01,t01,tta,ttl,taux,sig,n,yield_list,datediv_list,palier_list2,coeff_list2
n2=n_elements(palier_list2)
n2p=n_elements(coeff_list2)
if n2p ne n2+1 then print,'erreur dans le nombre de coefficients fournit:',n2p,'#',n2,'+1'
s=float(s01)
k=float(k01)
t=float(t01)
indexpalier=where(palier_list2 < k,nr)
if (nr eq 0) then return,ACRYLput(s,k,t,tta,r,sig,n,yield_list,datediv_list)
palier_list=palier_list2(indexpalier)
coeff_list=coeff_list2(indexpalier+1)
sum=0
for ri=1,nr do begin
	if ri eq nr then begin
		kk=k
		palier=palier_list(nr-1)
	endif else  begin
		kk=palier_list(ri)
		palier=palier_list(ri-1)
	endelse
sum1=ACRYLEput(s,kk,t,tta,ttl,taux,sig,n,yield_list,datediv_list,palier)
sum=sum+sum1*coeff_list(ri-1)
endfor
sum1=ACRYLput(s,palier_list(0),t,tta,taux,sig,n,yield_list,datediv_list)
sum=sum+sum1*coeff_list2(0)
return,sum
end


 

;*********************************************************************************************
;
;
;
;  formule de calculs d options impliquants des resets  et presentant un caratere americain
;
;
;
;
;*********************************************************************************************

;calcul d un call at the money forwade de td
;nc associe au cox
function OACRYLcall,td,s,k,t,ttad,r,sigma,n,nc,yieldlist,t_divlist
sup=exp((findgen(2*n+1)-n)/float(n)*3*sigma*sqrt(td))*s
vec=fltarr(2*n+1)
vec1=fltarr(2*n+1)
for i=0,2*n do begin
	vec(i)=ACRYLcall(sup(i),sup(i),t,ttad,r,sigma,nc,yieldlist,t_divlist)*exp(-((i-n)/float(n))^2)
	vec1(i)=exp(-((i-n)/float(n))^2)
endfor
return,total(vec)/total(vec1)*exp(-r*td)
end

;calcul d un put forwarde de td
function OACRYLput,td,s,k,t,ttad,r,sigma,n,nc,yieldlist,t_divlist
sup=exp((findgen(2*n+1)-n)/float(n)*3*sigma*sqrt(td))*s
vec=fltarr(2*n+1)
vec1=fltarr(2*n+1)
for i=0,2*n do begin
	vec(i)=ACRYLput(sup(i),sup(i),t,ttad,r,sigma,nc,yieldlist,t_divlist)*exp(-((i-n)/float(n))^2)
	vec1(i)=exp(-((i-n)/float(n))^2)
endfor
return,total(vec)/total(vec1)*exp(-r*td)
end

;calcul du call avec reset

function NACRYLRScall,s01,k01,t01,tta,taux,sig,n,nc,yield_list,datediv_list,datereset_list1
datereset_list=datereset_list1(sort(datereset_list1))
s=float(s01)
k=float(k01)
t=float(t01)
nr=n_elements(datereset_list)
if (nr eq 0) then return,ACRYLcall(s,k,t,tta,r,sig,nc,yield_list,datediv_list)
sum=0
for ri=0,nr do begin
	if ri eq 0 then begin
		td=0
		tf=datereset_list(ri)
	endif
	if ((ri gt 0) and (ri lt nr)) then begin
		td=datereset_list(ri-1)
		tf=datereset_list(ri)
	endif
	if ri eq nr then begin
		td=datereset_list(ri-1)
		tf=t
	endif
sum=sum+OACRYLcall(td,s,k,tf-td,tta-td,taux,sig,n,nc,yield_list,datediv_list-td)
endfor
return,sum
end

function NACRYLRSput,s01,k01,t01,tta,taux,sig,n,yield_list,datediv_list,datereset_list1
datereset_list=datereset_list1(sort(datereset_list1))
s=float(s01)
k=float(k01)
t=float(t01)
nr=n_elements(datereset_list)
if (nr eq 0) then return,ACRYLput(s,k,t,tta,r,sig,nc,yield_list,datediv_list)
sum=0
for ri=0,nr do begin
	if ri eq 0 then begin
		td=0
		tf=datereset_list(ri)
	endif
	if ((ri gt 0) and (ri lt nr)) then begin
		td=datereset_list(ri-1)
		tf=datereset_list(ri)
	endif
	if ri eq nr then begin
		td=datereset_list(ri-1)
		tf=t
	endif
sum=sum+OACRYLput(td,s,k,tf-td,tta-td,taux,sig,n,nc,yield_list,datediv_list-td)
endfor
return,sum
end

;*********************************************************************************
;
;
;
;                   options asiatique suivant l approximation de wilkinson
;                      d apres l article de levy (JIMF  1992 ,11,474-491)
;
;
;*********************************************************************************

; calcul des volstar,sstar,kstar equivalents pour le calcul des options asiatiques
function asiatique_call_param,St,At,k,Tech,l,g01,sigma,taux=taux,yieldlist=yieldlist,t_yieldlist=t_yieldlist
if keyword_set(taux) ne 0 then g=interpolation(taux,tech)/100. else g=g01
if l ge tech then begin 
	Emt=St/(g*l)*(exp(g*tech)-1.)
	Emt2=2*St^2/(l^2*(g+sigma^2))*((exp((2*g+sigma^2)*tech)-1.)/(2*g+sigma^2)-(exp(g*tech)-1.)/g)
 endif else begin
	Emt=St/(g*l)*exp(g*(tech-l))*(exp(g*l)-1.)
	Emt2=2*St^2/(l^2*(g+sigma^2))*exp((2*g+sigma^2)*(tech-l))*((exp((2*g+sigma^2)*l)-1.)/(2*g+sigma^2)-(exp(g*l)-1.)/g)
endelse
if l gt tech then begin
	kstar=k-At*(l-tech)/float(l)
	Mt=At-At*(l-tech)/float(l)
endif else begin
	kstar=k
	mt=at
endelse
alpha=2*alog(Emt)-0.5*alog(Emt2)
nu=sqrt(alog(Emt2)-2*alog(Emt))
d1= (0.5*alog(emt2)-alog(kstar))/nu
d2=d1-nu
value=exp(-g*tech)*(emt*gaussint(d1)-(kstar)*gaussint(d2))
mt=emt*exp(-g*tech)
impvol=impvol_call(value,mt,kstar,tech,g) 
return,[mt,kstar,impvol]
end

function asiatique_put_param,St,At,k,Tech,l,g01,sigma,taux=taux,yieldlist=yieldlist,t_yieldlist=t_yieldlist
if keyword_set(taux) ne 0 then g=interpolation(taux,tech)/100. else g=g01
if l ge tech then begin 
	Emt=St/(g*l)*(exp(g*tech)-1.)
	Emt2=2*St^2/(l^2*(g+sigma^2))*((exp((2*g+sigma^2)*tech)-1.)/(2*g+sigma^2)-(exp(g*tech)-1.)/g)
 endif else begin
	Emt=St/(g*l)*exp(g*(tech-l))*(exp(g*l)-1.)
	Emt2=2*St^2/(l^2*(g+sigma^2))*exp((2*g+sigma^2)*(tech-l))*((exp((2*g+sigma^2)*l)-1.)/(2*g+sigma^2)-(exp(g*l)-1.)/g)
endelse
if l gt tech then begin
	kstar=k-At*(l-tech)/float(l)
	Mt=At-At*(l-tech)/float(l)
endif else begin
	kstar=k
	mt=at
endelse
alpha=2*alog(Emt)-0.5*alog(Emt2)
nu=sqrt(alog(Emt2)-2*alog(Emt))
d1= (0.5*alog(emt2)-alog(kstar))/nu
d2=d1-nu
value=exp(-g*tech)*(emt*(gaussint(d1)-1.)-(kstar)*(gaussint(d2)-1.))
mt=emt*exp(-g*tech)
impvol=impvol_put(value,mt,kstar,tech,g) 
return,[mt,kstar,impvol]
end


; l est la longueur de la moyenne
; N est le nombre de point servant  a l etablissement de la moyenne
; tn est l echeance de l option
; t0 est le moment de calcul de la valeur de loption situe a tn - l
;  t est la date courrante de calcul
; St est la valeur du support en t
; At est la valeur de la moyenne en t si t < t0 sinon c est la valeur de la moyenne des cours utiles (de -(N-m) a 0) 
; si t0 est tel que t0 < t alors on peut definir m et zeta tel que
; si h = (tn - t0) / N
; t= t0+ m * h + zeta * h et |zeta| < 1.
; si t < t0 m et zeta ne seront pas defini 
function asiatique_call,St,At,k,Tech,l,N,g,sigma
	h=l/N
	egh=exp(g*h)
	esig2h=exp(sigma^2*h)
if l ge tech then begin 
PRINT,'BRANCHE1'
	t0 = tech-l
	past=l-tech
	m=int(past/h)
	zeta=past-m*h
	Emt=st/(n+1)*exp(g*(1-zeta)*h)*((1.-exp(g*(N-m)*h))/(1.-egh))
	a1=(exp((2*g+sigma^2)*h)-exp((2*g+sigma^2)*(N-m+1)*h))/((1.-egh)*(1.-egh^2*esig2h))
	a2=(exp((g*(N-m+2)+sigma^2)*h)-exp((2*g+sigma^2)*(N-m+1)*h))/((1.-egh)*(1.-egh*esig2h))
	a3=(exp((3*g+sigma^2)*h)-exp((g*(N-m+2)+sigma^2)*h))/((1.-egh)*(1.-egh*esig2h))
	a4=(exp(2*(2*g+sigma^2)*h)-exp((2*g+sigma^2)*(N-m+1)*h))/((1.-egh*esig2h)*(1.-egh^2*esig2h))
	Emt2=St^2/(N+1)^2*exp(-2*zeta*(g+1/(2*sigma^2))*h)*(a1-a2+a3-a4)
endif else begin
PRINT,'BRANCHE2'
	t0 = tech-l
	m=0
	Emt=st/(N+1)*exp(g*t0)*((1.-exp(g*(N+1)*h))/(1.-egh))
	b1=(1.-exp((2*g+sigma^2)*(N+1)*h))/((1.-egh)*(1.-egh^2*esig2h))
	b2=(exp(g*(N+1)*h)-exp((2*g+sigma^2)*(N+1)*h))/((1.-egh)*(1.-egh*esig2h))
	b3=(egh-exp(g*(N+1)*h))/((1.-egh)*(1.-egh*esig2h))
	b4=(exp((2*g+sigma^2)*h)-exp((2*g+sigma^2)*(N+1)*h))/((1.-egh*esig2h)*(1.-egh*egh*esig2h))
	Emt2=St^2/(N+1)^2*exp(2*(g+sigma^2)*t0)*(b1-b2+b3-b4)
endelse
print,'E(Mt)=',emt,'   E(Mt^2)=',emt2
print,' alors q un placement donerai ',At*exp(g*tech)
kstar=k-At*(m+1)/float(N+1)
mt=at-At*(m+1)/float(N+1)
;Mt=emt*exp(-(alpha+nu^2/2.)*tech)		|
;alpha=2*alog(Emt/Mt)-0.5*alog(Emt2/Mt^2)	|  a verifier a des fin de generalisation
;nu=sqrt(alog(Emt2/Mt^2)-2*alog(Emt/Mt))	|	(resoudre par iteration)
;value=call(Mt,Kstar,tech,alpha,nu)		|
alpha=2*alog(Emt)-0.5*alog(Emt2)
nu=sqrt(alog(Emt2)-2*alog(Emt))
d1= (0.5*alog(emt2)-alog(kstar))/nu
d2=d1-nu
value=exp(-g*tech)*(emt*gaussint(d1)-(kstar)*gaussint(d2))
print,'nouveau taux:',alpha,'  nouvelle volat:',nu,'  nouveau k:',kstar,'  nouveau sstar',mt
return,value
end

; version continue (N -> infini)
function asiatique_call_C,St,At,k,Tech,l,g,sigma
vec=asiatique_call_param(St,At,k,Tech,l,g,sigma)
print,'param=',vec
s=vec(0)
k=vec(1)
vol=vec(2)
return,call(s,k,tech,g,vol)
end

function asiatique_put_C,St,At,k,Tech,l,g,sigma
vec=asiatique_put_param(St,At,k,Tech,l,g,sigma)
print,'param=',vec
s=vec(0)
k=vec(1)
vol=vec(2)
return,put(s,k,tech,g,vol)
end


;*****************************************************************************************
;
;
;
;
;      diffusion avec deux jumps (un a la hausse et un a la baisse)
;
;
;
;
;
;
;****************************************************************************************
 
  
function jjACRcall,s,k,t,r,sigma,z1,lambda1,z2,lambda2,n
;call americain
;print,'param=',s,k,t,r,z1,lambda1,z2,lambda2,n,'vol=',sigma
nlimite=r^2*t/sigma^2
if n le nlimite then print,' n=',n,' est trop petit car n doit etre > ',nlimite
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
lambda=1.
a=exp(r*deltat)
call1=dblarr(2*n+1)
call2=dblarr(2*n+1)
supp=dblarr(2*n+1)
u=exp(double(sigma)*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
;print,'u=',u,' d=',d,' p=',p
p=(a-d)/(u-d)
pu=p*exp(-lambda1*t/float(n))
pd=(1-p)*exp(-lambda2*t/float(n))
pzd=(1-p)*(1.-exp(-lambda2*t/float(n)))
pzu=p*(1.-exp(-lambda1*t/float(n)))
kz1=int(alog(z1)/alog(u))
kz2=int(alog(z2)/alog(d))
;print,'kz1=',kz1,'  kz2=',kz2
;print,'p=',p,' pu=',pu,' pzu=',pzu
supp(n)=s
for j=1,n do begin
	supp(n+j)=supp(n+j-1)*u
endfor
for j=-1,-n,-1 do begin
	supp(n+j)=supp(n+j+1)*d
endfor
for j=-n,n do begin
	call2(n+j)=(supp(n+j)-k)>0
endfor
;print,'c=',call2
for j=-n,n do call1(n+j)=0
for i=n-1,0,-1 do begin
	for j=-(n-1),n-1 do begin
		case 1 of 
			j+kz1 gt n-2 and j-kz2 ge -(n-2) :call1(n+j)=supp(n+j)-k > exp(-r*deltat)*($
				pzu*((s*z1*u^(j+1)-k)>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
			j-kz2 lt -(n-2) and j+kz1 le n-2:call1(n+j)=supp(n+j)-k > exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*((s*z2*u^(j-1)-k)>0))
			j+kz1 gt n-2 and j-kz2 lt -(n-2):call1(n+j)=supp(n+j)-k > exp(-r*deltat)*($
				pzu*((s*z1*u^(j+1)-k)>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*((s*z2*u^(n+j-1)-k)>0))
			else:call1(n+j)=supp(n+j)-k > exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
		endcase
	endfor
;print,'c=',call1
	for j=-(n-1),n-1 do begin
		call2(n+j)=call1(n+j)
	endfor
endfor
print,'s=',s,' /k=',k,' /vol=', sigma,' /prime=',call1(n)
return,call1(n)
end
  
function jjECRcall,s,k,t,r,sigma,z1,lambda1,z2,lambda2,n
;call europeen
;print,'param=',s,k,t,r,z1,lambda1,z2,lambda2,n,'vol=',sigma
nlimite=r^2*t/sigma^2
if n le nlimite then print,' n=',n,' est trop petit car n doit etre > ',nlimite
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
lambda=1.
a=exp(r*deltat)
call1=dblarr(2*n+1)
call2=dblarr(2*n+1)
supp=dblarr(2*n+1)
u=exp(double(sigma)*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
;print,'u=',u,' d=',d,' p=',p
p=(a-d)/(u-d)
pu=p*exp(-lambda1*t/float(n))
pd=(1-p)*exp(-lambda2*t/float(n))
pzd=(1-p)*(1.-exp(-lambda2*t/float(n)))
pzu=p*(1.-exp(-lambda1*t/float(n)))
kz1=int(alog(z1)/alog(u))
kz2=int(alog(z2)/alog(d))
;print,'kz1=',kz1,'  kz2=',kz2
;print,'p=',p,' pu=',pu,' pzu=',pzu
supp(n)=s
for j=1,n do begin
	supp(n+j)=supp(n+j-1)*u
endfor
for j=-1,-n,-1 do begin
	supp(n+j)=supp(n+j+1)*d
endfor
for j=-n,n do begin
	call2(n+j)=(supp(n+j)-k)>0
endfor
;print,'c=',call2
for j=-n,n do call1(n+j)=0
for i=n-1,0,-1 do begin
	for j=-(n-1),n-1 do begin
		case 1 of 
			j+kz1 gt n-2 and j-kz2 ge -(n-2) :call1(n+j)= exp(-r*deltat)*($
				pzu*((s*z1*u^(j+1)-k*exp(-r*(t-i*deltat)))>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
			j-kz2 lt -(n-2) and j+kz1 le n-2:call1(n+j)= exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*((s*z2*u^(j-1)-k*exp(-r*(t-i*deltat)))>0))
			j+kz1 gt n-2 and j-kz2 lt -(n-2):call1(n+j)= exp(-r*deltat)*($
				pzu*((s*z1*u^(j+1)-k*exp(-r*(t-i*deltat)))>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+$
				pzd*((s*z2*u^(n+j-1)-k*exp(-r*(t-i*deltat)))>0))
			else:call1(n+j)= exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
		endcase
	endfor
;print,'c=',call1
	for j=-(n-1),n-1 do begin
		call2(n+j)=call1(n+j)
	endfor
endfor
print,'s=',s,' /k=',k,' /vol=', sigma,' /prime=',call1(n)
return,call1(n)
end



function jjACRput,s,k,t,r,sigma,z1,lambda1,z2,lambda2,n
;put americain
nlimite=r^2*t/sigma^2
if n le nlimite then print,' n=',n,' est trop petit car n doit etre > ',nlimite
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
lambda=1.
a=exp(r*deltat)
call1=dblarr(2*n+1)
call2=dblarr(2*n+1)
supp=dblarr(2*n+1)
u=exp(double(sigma)*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
;print,'u=',u,' d=',d,' p=',p
p=(a-d)/(u-d)
pu=p*exp(-lambda1*t/float(n))
pd=(1-p)*exp(-lambda2*t/float(n))
pzd=(1-p)*(1.-exp(-lambda2*t/float(n)))
pzu=p*(1.-exp(-lambda1*t/float(n)))
kz1=int(alog(z1)/alog(u))
kz2=int(alog(z2)/alog(d))
;print,'kz1=',kz1,'  kz2=',kz2
;print,'p=',p,' pu=',pu,' pzu=',pzu
supp(n)=s
for j=1,n do begin
	supp(n+j)=supp(n+j-1)*u
endfor
for j=-1,-n,-1 do begin
	supp(n+j)=supp(n+j+1)*d
endfor
for j=-n,n do begin
	call2(n+j)=(k-supp(n+j))>0
endfor
;print,'c=',call2
for j=-n,n do call1(n+j)=0
for i=n-1,0,-1 do begin
	for j=-(n-1),n-1 do begin
		case 1 of 
			j+kz1 gt n-2 and j-kz2 ge -(n-2) :call1(n+j)=k-supp(n+j) > exp(-r*deltat)*($
				pzu*((k-s*z1*u^(j+1))>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
			j-kz2 lt -(n-2) and j+kz1 le n-2:call1(n+j)=k-supp(n+j) > exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*((k-s*z2*u^(j-1))>0))
			j+kz1 gt n-2 and j-kz2 lt -(n-2):call1(n+j)=supp(n+j)-k > exp(-r*deltat)*($
				pzu*(k-(s*z1*u^(j+1))>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*((k-s*z2*u^(n+j-1))>0))
			else:call1(n+j)=k-supp(n+j) > exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
		endcase
	endfor
;print,'c=',call1
	for j=-(n-1),n-1 do begin
		call2(n+j)=call1(n+j)
	endfor
endfor
print,'s=',s,' /k=',k,' /vol=', sigma,' /prime=',call1(n)
return,call1(n)
end

function jjECRput,s,k,t,r,sigma,z1,lambda1,z2,lambda2,n
;put europeen
nlimite=r^2*t/sigma^2
if n le nlimite then print,' n=',n,' est trop petit car n doit etre > ',nlimite
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
lambda=1.
a=exp(r*deltat)
call1=dblarr(2*n+1)
call2=dblarr(2*n+1)
supp=dblarr(2*n+1)
u=exp(double(sigma)*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
;print,'u=',u,' d=',d,' p=',p
p=(a-d)/(u-d)
pu=p*exp(-lambda1*t/float(n))
pd=(1-p)*exp(-lambda2*t/float(n))
pzd=(1-p)*(1.-exp(-lambda2*t/float(n)))
pzu=p*(1.-exp(-lambda1*t/float(n)))
kz1=int(alog(z1)/alog(u))
kz2=int(alog(z2)/alog(d))
;print,'kz1=',kz1,'  kz2=',kz2
;print,'p=',p,' pu=',pu,' pzu=',pzu
supp(n)=s
for j=1,n do begin
	supp(n+j)=supp(n+j-1)*u
endfor
for j=-1,-n,-1 do begin
	supp(n+j)=supp(n+j+1)*d
endfor
for j=-n,n do begin
	call2(n+j)=(k-supp(n+j))>0
endfor
;print,'c=',call2
for j=-n,n do call1(n+j)=0
for i=n-1,0,-1 do begin
	for j=-(n-1),n-1 do begin
		case 1 of 
			j+kz1 gt n-2 and j-kz2 ge -(n-2) :call1(n+j)=exp(-r*deltat)*($
				pzu*((k*exp(-r*(t-i*deltat))-s*z1*u^(j+1))>0)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
			j-kz2 lt -(n-2) and j+kz1 le n-2:call1(n+j)=exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*((k*exp(-r*(t-i*deltat))-s*z2*u^(j-1))>0))
			j+kz1 gt n-2 and j-kz2 lt -(n-2):call1(n+j)=exp(-r*deltat)*($
				pzu*(k*exp(-r*(t-i*deltat))-(s*z1*u^(j+1))>0)+pu*call2(n+j+1)+$
				pd*call2(n+j-1)+pzd*((k*exp(-r*(t-i*deltat))-s*z2*u^(n+j-1))>0))
			else:call1(n+j)=exp(-r*deltat)*($
				pzu*call2(n+j+kz1+1)+pu*call2(n+j+1)+pd*call2(n+j-1)+pzd*call2(n+j-(kz2+1)))
		endcase
	endfor
;print,'c=',call1
	for j=-(n-1),n-1 do begin
		call2(n+j)=call1(n+j)
	endfor
endfor
print,'s=',s,' /k=',k,' /vol=', sigma,' /prime=',call1(n)
return,call1(n)
end



pro testd3
s=1800
k0=1800
t=0.5
n=200
volglob=0.23
z1=1.1
lambda1=0.2
z2=0.93
lambda2=0.
r=0.1
;p02=crcall(s,k2,t,r,volglob,n)
;voldiff2=impvol_jjcall(p02,s,k0,t,r,z1,lambda1,z2,lambda2,n)
;print,'pour avoir une vol globale de ',volglob,' il faut une vol diffusion de ',voldiff2
npoints=8
voldiff=0.15
kmin=1650.
kmax=1950.
karr=findgen(npoints)*(kmax-kmin)/(float(npoints)-1)+kmin
primearr=fltarr(npoints)
for i=0,npoints-1 do primearr(i)=jjAcrcall(s,karr(i),t,r,voldiff,z1,lambda1,z2,lambda2,n)
volexterne=fltarr(npoints)
for i=0,npoints-1 do volexterne(i)=impvol_call(primearr(i),s,karr(i),t,r)
window,0,title='vol calls'
plot,karr,volexterne
primearrp=fltarr(npoints)
for i=0,npoints-1 do primearrp(i)=jjAcrput(s,karr(i),t,r,voldiff,z1,lambda1,z2,lambda2,n)
volexternep=fltarr(npoints)
for i=0,npoints-1 do volexternep(i)=impvol_crput(primearrp(i),s,karr(i),t,r,n)
window,1,title=' vol puts (americain/americain)'
plot,karr,volexternep
volexternep1=fltarr(npoints)
for i=0,npoints-1 do volexternep1(i)=impvol_put(primearrp(i),s,karr(i),t,r)
window,2,title='puts (europeen/americain) '
plot,karr,volexternep1
end


pro testd31
s=1800
k0=1800
t=0.5
n=7
volglob=0.23
z1=1.1
lambda1=0.
z2=0.93
lambda2=1.
r=0.1
npoints=8
voldiff=0.15
k=1750.
primearrp=jjEcrput(s,k,t,r,voldiff,z1,lambda1,z2,lambda2,n)
pp=ECRput(s,k,t,r,voldiff,n)
print,'prime jj=',primearrp,'   /  prime =',pp
volexternep=impvol_crput(primearrp,s,k,t,r,n)
volexternep1=impvol_put(primearrp,s,k,t,r)
print,'volat. imp. americaine:',volexternep,'   / volat. imp. europeene :',volexternep1
end




;**************************************************************************************
;
;
;
;
;                  call offrant un rendement minimal
;
;
;
;
;**************************************************************************************
;calcul d'un call offrant a l'echeance en cas d exercice  le max{k*(1+g)^tt0,x * k*((sup(t)-sup(0)/sup(0)+1)}
;t est la distance restant a courrir jusqu'a l'echeance
;americain a partir de tta compte depuis le debut de l option 
; tt0 est le delai total entre de debut de l emission et l echeance
;rswap est le taux de swap a t c est a dire le coupons zero a l echeance
function AKCRYLcall,s,k,t,tt0,tta,r,sigma,n,yieldlist,t_divlist,g,x
deltat=max([float(t)/n,0])
if deltat eq 0 then begin 
		value=max([s,k])
		goto,next
		endif
a=exp(r*deltat)
call1=fltarr(n+1)
call2=fltarr(n+1)
supp=fltarr(n+1)
delta=fltarr(n+1)
u=exp(sigma*sqrt(deltat))
d=1/u
p=(a-d)/(u-d)
print,'p=',p
indexlist=where(t_divlist lt t,count)
	tte=0
if count eq 0 then begin
	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=(x*k*(supp(i)-s)/s)>(k*((1+g)^tt0-1))
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i)) 
		call1(i)=(x*k*(supp(i)*pp-s)/s)>k*((1+g)^tt0-1)
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	tte=deltat*(n-1-i)
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
if i*deltat lt (tt0-tta) then begin
		for j=0,i do call2(j)=((x*k*(s*u^j*d^(i-j)-s)/s)>k*((1+g)^(tt0-tte)-1.))> $
					exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
print,'h1a=',exp(r*deltat)*call2(0:i)
			endif else begin
		for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
print,'h1e=',exp(r*deltat)*call2(0:i)
			endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
if i*deltat lt (tt0-tta) then begin
		for j=0,i do call2(j)=((x*k*(s*u^j*d^(i-j)*pp-s)/s)>k*((1+g)^(tt0-tte)-1.))> $
					exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
print,'h2a=',exp(r*deltat)*call2(0:i)
			endif else begin
		for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
print,'h2e=',exp(r*deltat)*call2(0:i)
			endelse
	endelse
	for j=0,i do call1(j)=call2(j)
print,call1(0:i)
endfor
value=call2(0)
next: ;suite du calcul
return,value
end

;****************************************************************************************
;
;
;
;
;                    calcul des obligations convertibles
;
;
;
;
;****************************************************************************************

;calcul d'une obligation convertible offrant  max(s-kt,exp(rt)*opt(t+1))+oblig(taux,echeance)
;lorsque l on exerce ou kt =oblig(taux,echeante-t) 
;americain a partir de tta (distance depuis 0)
;t est la distance restant a courrir jusqu'a l'echeance
;g est le taux garanti  exemple 0.045
;x est le levrage de la convertible exemple 2.
;kn est le nominal de l emprunt
; exemple de prime : 200
; exemple de trigger : 1.3
; suppose que les temps sont exprime en annees
function ACRYLconv,s,kn,t,tta,r,sigma,n,yieldlist,t_divlist,g,t_premiercoupon,x,prime,calltrigger
;si s>1.3*nominal/x alors on exerce a exactement 1.3*nominal/x
;calcul du plancher actuariel
smax=calltrigger*kn/x
obn=plancheractuariel(g,r,t-t_premiercoupon,kn,prime,t_premiercoupon)
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
		k=obn
		call1(i)=x*supp(i)-k>0
		if supp(i) gt smax then call1(i)=x*smax-k>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i))
		k=obn 
		call1(i)=x*supp(i)*pp-k>0
		if supp(i)*pp gt smax then call1(i)=x*smax-k>0
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin

	tte=deltat*(n-1-i)
	if tte lt t-t_premiercoupon then begin
		et_premiercoupon=t_premiercoupon-(t-tte)
		ntte=t-t_premiercoupon
	endif else begin
		ntte=int(tte)
		et_premiercoupon=tte-ntte
	endelse
	k=plancheractuariel(g,r,ntte,kn,prime,et_premiercoupon)

	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		if i*deltat gt tta then begin
			for j=0,i do begin

					call2(j)=x*s*u^j*d^(i-j)-k> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
			if s*u^j*d^(i-j) gt smax then call2(j)=x*smax-k>0
				endfor
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))

		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if i*deltat gt tta then begin
			for j=0,i do begin
					call2(j)=x*s*u^j*d^(i-j)*pp-k> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
			if s*u^j*d^(i-j)*pp gt smax then call2(j)=x*smax-k>0

			endfor
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))

		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
print,'plancher actuariel:',obn
return,call2(0)+obn
end

;calcul d'une obligation convertible offrant  max(s-kt,exp(rt)*opt(t+1))+oblig(taux,echeance)
;lorsque l on exerce ou kt =oblig(taux,echeante-t) 
;americain a partir de tta (distance depuis 0)
; avec deux call trigger separe par une date (delai depuis le debut) datetrigger
;t est la distance restant a courrir jusqu'a l'echeance
;g est le taux garanti 
;x est le levrage de la convertible
;kn est le nominal de l emprunt
function ACRYLconv2,s,kn,t,tta,r,sigma,n,yieldlist,t_divlist,g,t_premiercoupon,x,prime,calltrigger1,calltrigger2,datetrigger,flagtrigger
;si s>1.3*nominal*g alors on exerce a exactement 1.3*g*nominal
;calcul du plancher actuariel
smax1=calltrigger1*kn/g
smax2=calltrigger2*kn/g
;obn=plancheractuariel(x,r,t,kn,prime)
obn=plancheractuariel(g,r,t,kn,prime,t_premiercoupon)
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
		k=obn
		call1(i)=x*supp(i)-k>0
		if supp(i) gt smax2 then call1(i)=x*smax2-k>0
	endfor
endif else begin
	ylist=yieldlist(indexlist)
	pp=1.
	for i=0,count-1 do pp=pp*(1-ylist(i))
 	for i=0,n do begin
		supp(i)=s*(u^i)*(d^(n-i))
		k=obn 
		call1(i)=g*supp(i)*pp-k>0
		if supp(i)*pp gt smax2 then call1(i)=x*smax2-k>0
	endfor
endelse
for j=0,n do call2(j)=0
for i=n-1,0,-1 do begin
	tte=deltat*(n-1-i)
	if flagtrigger eq 0 then if i*deltat gt datetrigger then smax=smax2 else smax=smax1 $
		 else if i*deltat gt datetrigger then $
		 smax=(smax2-smax1)*(i*deltat-datetrigger)/(t-datetrigger)+smax1 $
		else smax=smax1
	if tte lt t-t_premiercoupon then begin
		et_premiercoupon=t_premiercoupon-(t-tte)
		ntte=t-t_premiercoupon
	endif else begin
		ntte=int(tte)
		et_premiercoupon=tte-ntte
	endelse
	k=plancheractuariel(g,r,ntte,kn,prime,et_premiercoupon)	
	indexlist=where(t_divlist lt i*deltat,count)
	if count eq 0 then  begin 
		if i*deltat gt tta then begin
			for j=0,i do begin
					call2(j)=x*s*u^j*d^(i-j)-k> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
			if s*u^j*d^(i-j) gt smax then call2(j)=x*smax-k>0
				endfor
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))

		endelse
	endif else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
		if i*deltat gt tta then begin
			for j=0,i do begin
					call2(j)=x*s*u^j*d^(i-j)*pp-k> $
				exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))
			if s*u^j*d^(i-j)*pp gt smax then call2(j)=x*smax-k>0

			endfor
		endif else begin
			for j=0,i do call2(j)=exp(-r*deltat)*(p*call1(j+1)+(1-p)*call1(j))

		endelse

	endelse
	for j=0,i do call1(j)=call2(j)
endfor
print,'plancher actuariel:',obn
return,call2(0)+obn
end

;****************************************************************************************
;
;
;
;              calcul des obligations indexees
;
;
;
;
;****************************************************************************************


;calcul d'une obligation indexee offrant  
;t est la distance restant a courrir jusqu'a l'echeance
;g est le levrage de l'indexe
;kn est le nominal de l emprunt
;cap est le cap de remboursement en terme de captage du support
;pcgarantie est le pourcentage du capital garanti a l echeance : ex pcgarantie=1. ou 0.95
;sigma est la volatilite du call a la monnaie
;sigma2 est la volatilite du call out
function YLindex,s,kn,t,r,sigma,sigma2,yieldlist,t_divlist,g,cap,pcgarantie
obn=kn/(1+r)^t
print,'plancher actuariel :',obn
k2=s+s/g*(cap-obn*pcgarantie/kn)
w1=YLcall(s,s,t,r,sigma,yieldlist,t_divlist)
w2=YLcall(s,k2,t,r,sigma2,yieldlist,t_divlist)
print,'w1=',100*w1/s,'%,   w2=',100*w2/s,'%'
print,'k2=',100*k2/s,'%'
return,pcgarantie*obn+g*(kn/s)*(w1-w2)
end


;calcul d'une obligation indexee offrant  
;t est la distance restant a courrir jusqu'a l'echeance
;g est le levrage de l'indexe
;kn est le nominal de l emprunt
;cap est le cap de remboursement en terme de captage du support
;pcgarantie est le pourcentage du capital garanti a l echeance : ex pcgarantie=1. ou 0.95
;sigma est la volatilite du call a la monnaie
;sigma2 est la volatilite du call out
function DLindex,s,kn,t,r,sigma,sigma2,divlist,t_divlist,g,cap,pcgarantie
obn=kn/(1+r)^t
print,'plancher actuariel :',obn
k2=s+s/g*(cap-obn*pcgarantie/kn)
w1=DLcall(s,s,t,r,sigma,divlist,t_divlist)
w2=DLcall(s,k2,t,r,sigma2,divlist,t_divlist)
print,'w1=',100*w1/s,'%,   w2=',100*w2/s,'%'
print,'k2=',100*k2/s,'%'
return,pcgarantie*obn+g*(kn/s)*(w1-w2)
end


;***************************************************************************************
;
;
;
;            calcul d options sur deux support
;
;
;
;
;***************************************************************************************

 
function solveu2_u1,u1,M,V1,V2
x1=(double(1.) - 2*M + M^2 + u1^2 - 2*M*u1^2 + M^2*u1^2 + 2*u1*V1 + V2 -2*u1*V2 + u1^2*V2 - sqrt(1. - 2*M + M^2 + V2)*$
    sqrt(1 - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 +M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2) +$
    u1*sqrt(1. - 2*M + M^2 + V2)*sqrt(1 - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 + $
    M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2))/(2*u1*(1. - 2*M + M^2 + V1))

x2=(double(1.) - 2*M + M^2 + u1^2 - 2*M*u1^2 + M^2*u1^2 + 2*u1*V1 + V2 -2*u1*V2 + u1^2*V2 + sqrt(1. - 2*M + M^2 + V2)*$
    sqrt(1. - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 +M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2) -$
    u1*sqrt(1 - 2*M + M^2 + V2)*sqrt(1 - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 +$ 
    M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2))/(2*u1*(1. - 2*M + M^2 + V1))
return,[x1,x2]
end


function CRYLcallmax,s1,s2,k,t,r,sigma1,sigma2,corr,n,yieldlist1,t_divlist1,yieldlist2,t_divlist2
lambda=1.1
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
call1=fltarr(2*n+1,2*n+1)
call2=fltarr(2*n+1,2*n+1)
supp1=fltarr(2*n+1)
supp2=fltarr(2*n+1)
M=(exp(r*double(deltat)))
V1=(exp(double(sigma1)^2*deltat)-1.)*M^2
V2=(exp(double(sigma2)^2*deltat)-1.)*M^2
RR=(exp(corr*double(sigma1)*double(sigma2)*double(deltat)))*M^2
u1=exp(lambda*double(sigma1)*sqrt(double(deltat)))
d1=1./u1
ku2=solveu2_u1(u1,M,V1,V2)
kku2=where(ku2 ge 1.)
u2=min(ku2(kku2))
d2=1./u2
f1=((V1+M^2-M)*u1-(M-1.))/((u1-1.)*(u1^2-1.))
f2=((V2+M^2-M)*u2-(M-1.))/((u2-1.)*(u2^2-1.))
g1=((u1^2)*(V1+M^2-M)-(u1^3)*(M-1.))/((u1-1.)*(u1^2-1.))
g2=((u2^2)*(V2+M^2-M)-(u2^3)*(M-1.))/((u2-1.)*(u2^2-1.))
p1=(u1*u2*(rr-1.)-f1*(u1^2-1.)-f2*(u2^2-1.)+(f2+g2)*(u1*u2-1.))/((u1^2-1)*(u2^2-1))
p2=f1-p1
p2p=(f1*(u1^2-1.)*u2^2+f2*(u2^2-1.)-(f2+g2)*(u1*u2-1.)-u1*u2*(rr-1.))/((u1^2-1.)*(u2^2-1.))
p3=g2-p2
p4=f2-p1
p5=1-(f1+g1)
;print,'prob=',p1,p2,p3,p4,p5
indexlist1=where(t_divlist1 lt t,count1)
indexlist2=where(t_divlist2 lt t,count2)
for i=0,2*n do begin
		supp1(i)=s1*u1^(i-n) 
		supp2(i)=s2*u2^(i-n)
endfor
pp1=1.
pp2=1.
if count1 ne 0 then begin
		ylist1=yieldlist1(indexlist1)
		for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
endif
if count2 ne 0 then begin
		ylist2=yieldlist2(indexlist2)
		for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
endif
epsilon=-1
for j1=-n,n do begin
	nn=(n-((1+epsilon)/2))
	for j2=-nn,nn,2 do begin
		call2(n+j1,n+j2)=max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>0
		endfor
	epsilon=-epsilon
endfor
for i=n-1,0,-1 do begin
	call1(0:2*n,0:2*n)=0
	indexlist1=where(t_divlist1 lt i*deltat,count1)
	indexlist2=where(t_divlist2 lt i*deltat,count2)
	pp1=1.
	pp2=1.
	if count1 ne 0 then begin
			ylist1=yieldlist1(indexlist1)
			for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
	endif
	if count2 ne 0 then begin
			ylist2=yieldlist2(indexlist2)
			for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
	endif
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call1(n+j1,n+j2)=max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>$
			1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		endfor
		epsilon=-epsilon
	endfor
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call2(n+j1,n+j2)=call1(n+j1,n+j2)
		endfor
		epsilon=-epsilon
	endfor
endfor
return,call2(n,n)
end
;***************************************************************************************
;
;
;
;            calcul du call rendement differentiel de deux support S1-S2 > 0
;
;
;
;
;***************************************************************************************

function CRYLcalldiff,s1,s2,t,r,sigma1,sigma2,corr,n,yieldlist1,t_divlist1,yieldlist2,t_divlist2
lambda=1.1
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
call1=fltarr(2*n+1,2*n+1)
call2=fltarr(2*n+1,2*n+1)
supp1=fltarr(2*n+1)
supp2=fltarr(2*n+1)
M=(exp(r*double(deltat)))
V1=(exp(double(sigma1)^2*deltat)-1.)*M^2
V2=(exp(double(sigma2)^2*deltat)-1.)*M^2
RR=(exp(corr*double(sigma1)*double(sigma2)*double(deltat)))*M^2
u1=exp(lambda*double(sigma1)*sqrt(double(deltat)))
d1=1./u1
ku2=solveu2_u1(u1,M,V1,V2)
kku2=where(ku2 ge 1.)
u2=min(ku2(kku2))
d2=1./u2
f1=((V1+M^2-M)*u1-(M-1.))/((u1-1.)*(u1^2-1.))
f2=((V2+M^2-M)*u2-(M-1.))/((u2-1.)*(u2^2-1.))
g1=((u1^2)*(V1+M^2-M)-(u1^3)*(M-1.))/((u1-1.)*(u1^2-1.))
g2=((u2^2)*(V2+M^2-M)-(u2^3)*(M-1.))/((u2-1.)*(u2^2-1.))
p1=(u1*u2*(rr-1.)-f1*(u1^2-1.)-f2*(u2^2-1.)+(f2+g2)*(u1*u2-1.))/((u1^2-1)*(u2^2-1))
p2=f1-p1
p2p=(f1*(u1^2-1.)*u2^2+f2*(u2^2-1.)-(f2+g2)*(u1*u2-1.)-u1*u2*(rr-1.))/((u1^2-1.)*(u2^2-1.))
p3=g2-p2
p4=f2-p1
p5=1-(f1+g1)
;print,'prob=',p1,p2,p3,p4,p5
indexlist1=where(t_divlist1 lt t,count1)
indexlist2=where(t_divlist2 lt t,count2)
for i=0,2*n do begin
		supp1(i)=s1*u1^(i-n) 
		supp2(i)=s2*u2^(i-n)
endfor
pp1=1.
pp2=1.
if count1 ne 0 then begin
		ylist1=yieldlist1(indexlist1)
		for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
endif
if count2 ne 0 then begin
		ylist2=yieldlist2(indexlist2)
		for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
endif
epsilon=-1
for j1=-n,n do begin
	nn=(n-((1+epsilon)/2))
	for j2=-nn,nn,2 do begin
		call2(n+j1,n+j2)=supp1(n+j1)*pp1-supp2(n+j2)*pp2>0
		endfor
	epsilon=-epsilon
endfor
for i=n-1,0,-1 do begin
	call1(0:2*n,0:2*n)=0
	indexlist1=where(t_divlist1 lt i*deltat,count1)
	indexlist2=where(t_divlist2 lt i*deltat,count2)
	pp1=1.
	pp2=1.
	if count1 ne 0 then begin
			ylist1=yieldlist1(indexlist1)
			for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
	endif
	if count2 ne 0 then begin
			ylist2=yieldlist2(indexlist2)
			for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
	endif
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call1(n+j1,n+j2)=supp1(n+j1)*pp1-supp2(n+j2)*pp2>$
			1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		endfor
		epsilon=-epsilon
	endfor
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call2(n+j1,n+j2)=call1(n+j1,n+j2)
		endfor
		epsilon=-epsilon
	endfor
endfor
return,call2(n,n)
end

;***************************************************************************************
;
;
;
;            calcul de la gestion simultanee de deux calls
;
;
;
;
;***************************************************************************************

function CRYLcall2S,s1,s2,k1,k2,t,r,sigma1,sigma2,corr,n,yieldlist1,t_divlist1,yieldlist2,t_divlist2
lambda=1.2
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
call1=fltarr(2*n+1,2*n+1)
call2=fltarr(2*n+1,2*n+1)
supp1=fltarr(2*n+1)
supp2=fltarr(2*n+1)
M=(exp(r*double(deltat)))
V1=(exp(double(sigma1)^2*deltat)-1.)*M^2
V2=(exp(double(sigma2)^2*deltat)-1.)*M^2
RR=(exp(corr*double(sigma1)*double(sigma2)*double(deltat)))*M^2
u1=exp(lambda*double(sigma1)*sqrt(double(deltat)))
d1=1./u1
ku2=solveu2_u1(u1,M,V1,V2)
kku2=where(ku2 ge 1.)
u2=min(ku2(kku2))
d2=1./u2
f1=((V1+M^2-M)*u1-(M-1.))/((u1-1.)*(u1^2-1.))
f2=((V2+M^2-M)*u2-(M-1.))/((u2-1.)*(u2^2-1.))
g1=((u1^2)*(V1+M^2-M)-(u1^3)*(M-1.))/((u1-1.)*(u1^2-1.))
g2=((u2^2)*(V2+M^2-M)-(u2^3)*(M-1.))/((u2-1.)*(u2^2-1.))
p1=(u1*u2*(rr-1.)-f1*(u1^2-1.)-f2*(u2^2-1.)+(f2+g2)*(u1*u2-1.))/((u1^2-1)*(u2^2-1))
p2=f1-p1
p2p=(f1*(u1^2-1.)*u2^2+f2*(u2^2-1.)-(f2+g2)*(u1*u2-1.)-u1*u2*(rr-1.))/((u1^2-1.)*(u2^2-1.))
p3=g2-p2
p4=f2-p1
p5=1-(f1+g1)
print,'prob=',p1,p2,p3,p4,p5
indexlist1=where(t_divlist1 lt t,count1)
indexlist2=where(t_divlist2 lt t,count2)
for i=0,2*n do begin
		supp1(i)=s1*u1^(i-n) 
		supp2(i)=s2*u2^(i-n)
endfor
pp1=1.
pp2=1.
if count1 ne 0 then begin
		ylist1=yieldlist1(indexlist1)
		for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
endif
if count2 ne 0 then begin
		ylist2=yieldlist2(indexlist2)
		for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
endif
epsilon=-1
for j1=-n,n do begin
	nn=(n-((1+epsilon)/2))
	for j2=-nn,nn,2 do begin
		call2(n+j1,n+j2)=((supp1(n+j1)*pp1-k1)>0)+((supp2(n+j2)*pp2-k2)>0)
		endfor
	epsilon=-epsilon
endfor
for i=n-1,0,-1 do begin
	call1(0:2*n,0:2*n)=0
	indexlist1=where(t_divlist1 lt i*deltat,count1)
	indexlist2=where(t_divlist2 lt i*deltat,count2)
	pp1=1.
	pp2=1.
	if count1 ne 0 then begin
			ylist1=yieldlist1(indexlist1)
			for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
	endif
	if count2 ne 0 then begin
			ylist2=yieldlist2(indexlist2)
			for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
	endif
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call1(n+j1,n+j2)=((supp1(n+j1)*pp1-k1)>0)+((supp2(n+j2)*pp2-k2)>0)>$
			1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		endfor
		epsilon=-epsilon
	endfor
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call2(n+j1,n+j2)=call1(n+j1,n+j2)
		endfor
		epsilon=-epsilon
	endfor
endfor
return,call2(n,n)
end

pro testh,rho
a=crylcall(100.,100.,4.,0.1,0.3,100,replicate(0.025,7),0.5+findgen(7))
b=crylcall(100.,100.,4.,0.1,0.3,100,replicate(0.025,7),0.5+findgen(7))
a1=ylcall(100.,100.,4.,0.1,0.3,replicate(0.025,7),0.5+findgen(7))
b1=ylcall(100.,100.,4.,0.1,0.3,replicate(0.025,7),0.5+findgen(7))
c=crylcall2s(100,100,100.,100.,4.,0.1,0.3,0.3,rho,30,replicate(0.025,7),0.5+findgen(7),replicate(0.025,7),0.5+findgen(7))
print,a1+b1,a+b,c
end


;***************************************************************************************
;
;
;
;            calcul du call  differentiel  a plancher de deux support Min{S1-S2>0,S1-k>0} si 
;
;
;
;
;***************************************************************************************

function CRYLcalldifC,s1,s2,k,t,r,sigma1,sigma2,corr,n,yieldlist1,t_divlist1,yieldlist2,t_divlist2
lambda=1.1
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
call1=fltarr(2*n+1,2*n+1)
call2=fltarr(2*n+1,2*n+1)
supp1=fltarr(2*n+1)
supp2=fltarr(2*n+1)
M=(exp(r*double(deltat)))
V1=(exp(double(sigma1)^2*deltat)-1.)*M^2
V2=(exp(double(sigma2)^2*deltat)-1.)*M^2
RR=(exp(corr*double(sigma1)*double(sigma2)*double(deltat)))*M^2
u1=exp(lambda*double(sigma1)*sqrt(double(deltat)))
d1=1./u1
ku2=solveu2_u1(u1,M,V1,V2)
kku2=where(ku2 ge 1.)
u2=min(ku2(kku2))
d2=1./u2
f1=((V1+M^2-M)*u1-(M-1.))/((u1-1.)*(u1^2-1.))
f2=((V2+M^2-M)*u2-(M-1.))/((u2-1.)*(u2^2-1.))
g1=((u1^2)*(V1+M^2-M)-(u1^3)*(M-1.))/((u1-1.)*(u1^2-1.))
g2=((u2^2)*(V2+M^2-M)-(u2^3)*(M-1.))/((u2-1.)*(u2^2-1.))
p1=(u1*u2*(rr-1.)-f1*(u1^2-1.)-f2*(u2^2-1.)+(f2+g2)*(u1*u2-1.))/((u1^2-1)*(u2^2-1))
p2=f1-p1
p2p=(f1*(u1^2-1.)*u2^2+f2*(u2^2-1.)-(f2+g2)*(u1*u2-1.)-u1*u2*(rr-1.))/((u1^2-1.)*(u2^2-1.))
p3=g2-p2
p4=f2-p1
p5=1-(f1+g1)
;print,'prob=',p1,p2,p3,p4,p5
indexlist1=where(t_divlist1 lt t,count1)
indexlist2=where(t_divlist2 lt t,count2)
for i=0,2*n do begin
		supp1(i)=s1*u1^(i-n) 
		supp2(i)=s2*u2^(i-n)
endfor
pp1=1.
pp2=1.
if count1 ne 0 then begin
		ylist1=yieldlist1(indexlist1)
		for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
endif
if count2 ne 0 then begin
		ylist2=yieldlist2(indexlist2)
		for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
endif
epsilon=-1
for j1=-n,n do begin
	nn=(n-((1+epsilon)/2))
	for j2=-nn,nn,2 do begin
		call2(n+j1,n+j2)=Min([supp1(n+j1)*pp1-supp2(n+j2)*pp2>0,supp1(n+j1)*pp1-k>0])
		endfor
	epsilon=-epsilon
endfor
for i=n-1,0,-1 do begin
	call1(0:2*n,0:2*n)=0
	indexlist1=where(t_divlist1 lt i*deltat,count1)
	indexlist2=where(t_divlist2 lt i*deltat,count2)
	pp1=1.
	pp2=1.
	if count1 ne 0 then begin
			ylist1=yieldlist1(indexlist1)
			for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
	endif
	if count2 ne 0 then begin
			ylist2=yieldlist2(indexlist2)
			for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
	endif
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call1(n+j1,n+j2)=Min([supp1(n+j1)*pp1-supp2(n+j2)*pp2>0,supp1(n+j1)*pp1-k>0])>$
			1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		endfor
		epsilon=-epsilon
	endfor
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call2(n+j1,n+j2)=call1(n+j1,n+j2)
		endfor
		epsilon=-epsilon
	endfor
endfor
return,call2(n,n)
end

;***********************************************************************************************
;
;
;                  Fonction  a deux supports
;
;
;***********************************************************************************************

function Norm2d1,x,y,ap,bp,ro
	return,exp(ap*(2*x-ap) + bp*(2*y-bp) + 2*ro*(x-ap)*(y-bp))
end

function Sign,x
	if (x  ge  0) then return,1. else return,-1.
end

function Norm2d2,a,b,ro
; lorsque a le 0 et b le 0 et ro le 0 
	A1=double(0.3253030)
	A2=double(0.4211071)
	A3=double(0.1334425)
	A4=double(0.006374323)
	B1=double(0.1337764)
	B2=double(0.6243247)
	B3=double(1.3425378)
	B4=double(2.2626645)
	pi=double(3.1415926535)
	s=double(0.)
	ap=double(0.)
	bp=double(0.)
	ap=a/sqrt(2.-2.*ro*ro);
	bp=b/sqrt(2.-2.*ro*ro);
	s=(A1*Norm2d1(B1,B1,ap,bp,ro) $
	 + A2*Norm2d1(B1,B2,ap,bp,ro) $
	 + A3*Norm2d1(B1,B3,ap,bp,ro) $
	 + A4*Norm2d1(B1,B4,ap,bp,ro))*A1 $
	+ (A1*Norm2d1(B2,B1,ap,bp,ro) $
	 + A2*Norm2d1(B2,B2,ap,bp,ro) $
	 + A3*Norm2d1(B2,B3,ap,bp,ro) $
	 + A4*Norm2d1(B2,B4,ap,bp,ro))*A2 $
	+ (A1*Norm2d1(B3,B1,ap,bp,ro) $
	 + A2*Norm2d1(B3,B2,ap,bp,ro) $
	 + A3*Norm2d1(B3,B3,ap,bp,ro) $
	 + A4*Norm2d1(B3,B4,ap,bp,ro))*A3 $
	+ (A1*Norm2d1(B4,B1,ap,bp,ro) $
	 + A2*Norm2d1(B4,B2,ap,bp,ro) $
	 + A3*Norm2d1(B4,B3,ap,bp,ro) $
	 + A4*Norm2d1(B4,B4,ap,bp,ro))*A4
	return,s*(sqrt(1-ro*ro))/pi
end

function Norm2d3,a,b,ro
;lorsque b=0  
	if ((a le 0)   and  (ro le 0)) then return,Norm2d2(a,b,ro)
	if ((a le 0)   and  (ro gt 0)) then return,gaussint(a)-Norm2d2(a,b,-ro)
	if ((a gt 0)   and  (ro gt 0)) then return,gaussint(b)-Norm2d2(-a,b,-ro)
	if ((a gt 0)   and  (ro le 0)) then return,gaussint(b)+gaussint(a)+Norm2d2(-a,b,ro)-1.0
end

function Norm2d4,a,b,ro
;lorsque a=0  
	if ((b le 0)   and  (ro le 0)) then return,Norm2d2(a,b,ro)
	if ((b le 0)   and  (ro gt 0)) then return,gaussint(b)-Norm2d2(b,a,-ro)
	if ((b gt 0)   and  (ro gt 0)) then return,gaussint(a)-Norm2d2(-b,a,-ro)
	if ((b gt 0)   and  (ro le 0)) then return,gaussint(a)+gaussint(b)+Norm2d2(-b,a,ro)-1.0
end

function Norm2d,a,b,ro
	l=0.0
	case 1 of
	(b eq 0) :return,Norm2d3(a,b,ro)
	(a eq 0) :return,Norm2d4(a,b,ro)
	((a le 0)   and  (b le 0)  and  (ro le 0)): return,Norm2d2(a,b,ro) 
	((a le 0)   and  (b gt 0)  and  (ro gt 0)): return,gaussint(a)-Norm2d2(a,-b,-ro)
	((a gt 0)   and  (b le 0)  and  (ro gt 0)): return,gaussint(b)-Norm2d2(-a,b,-ro)
	((a gt 0)   and  (b gt 0)  and  (ro le 0)): return,gaussint(a)+gaussint(b)-1.+Norm2d2(-a,-b,ro)
	else: begin
		ro1 = ((ro*a - b) * Sign(a)) / sqrt(a*a - 2.*ro*a*b + b*b)
		ro2 = ((ro*b - a) * Sign(b)) / sqrt(a*a - 2.*ro*a*b + b*b)
		delta1 = (1.-Sign(a)*Sign(b)) / 4.
		return,Norm2d3(a,l,ro1) + Norm2d3(b,l,ro2) - delta1
	end
	endcase
end

function Norm2d2a,a,b,ro
	abcisse=dblarr(8)
	poids=dblarr(8)
	pi=double(3.1415926535)
	s=0.0
	abcisse(0)=	5.29786439318514e-2
	abcisse(1)=	2.67398372167767e-1
	abcisse(2)=	6.16302884182402e-1
	abcisse(3)=	1.06424631211623
	abcisse(4)=	1.58885586227006
	abcisse(5)=	2.18392115309586
	abcisse(6)=	2.86313388370808
	abcisse(7)=	3.68600716272440
	poids(0)=	1.34109188453360e-1
	poids(1)=	2.68330754472640e-1
	poids(2)=	2.75953397988422e-1
	poids(3)=	1.57448282618790e-1
	poids(4)=	4.48141099174625e-2
	poids(5)=	5.36793575602526e-3
	poids(6)=	2.02063649132407e-4
	poids(7)=	1.19259692659532e-6
	ap=a/sqrt(2-2*ro*ro)
	bp=b/sqrt(2-2*ro*ro)
	for i=0,7 do for j=0,7 do s=s+poids(i)*poids(j)*Norm2d1(abcisse(i),abcisse(j),ap,bp,ro)
	return,s*(sqrt(1-ro*ro))/pi
end

function Norm2d3a,a,b,ro
 ; lorsque b=0  
	if (a le 0   and  ro le 0) then return,Norm2d2a(a,b,ro)
	if (a le 0   and  ro gt 0) then return,  gaussint(a)-Norm2d2a(a,b,-ro)
	if (a gt 0   and  ro gt 0) then return,  gaussint(b)-Norm2d2a(-a,b,-ro)
	if (a gt 0   and  ro le 0) then return,  gaussint(b)+  gaussint(a)+Norm2d2a(-a,b,ro)-1.0
end

function Norm2d4a,a,b,ro
;lorsque a=0  
	if ((b le 0)   and  (ro le 0)) then return,Norm2d2a(a,b,ro)
	if ((b le 0)   and  (ro gt 0)) then return,gaussint(b)-Norm2d2a(b,a,-ro)
	if ((b gt 0)   and  (ro gt 0)) then return,gaussint(a)-Norm2d2a(-b,a,-ro)
	if ((b gt 0)   and  (ro le 0)) then return,gaussint(a)+gaussint(b)+Norm2d2a(-b,a,ro)-1.0
end

function Norm2da,a,b,ro
	l=0.0
	case 1 of
	(b eq 0) :return,Norm2d3a(a,b,ro)
	(a eq 0) :return,Norm2d4a(a,b,ro)
	(a le 0   and  b le 0   and  ro le 0): return,Norm2d2a(a,b,ro)
	(a le 0   and  b gt 0   and  ro gt 0): return,  gaussint(a)-Norm2d2a(a,-b,-ro)
	(a gt 0   and  b le 0   and  ro gt 0): return,  gaussint(b)-Norm2d2a(-a,b,-ro)
	(a gt 0   and  b gt 0   and  ro le 0): return,  gaussint(a)+  gaussint(b)-1+Norm2d2a(-a,-b,ro)
	 else : begin
			  ro1 = ((ro*a - b) * Sign(a)) / sqrt (a*a - 2.*ro*a*b + b*b)
			  ro2 = ((ro*b - a) * Sign(b)) / sqrt (a*a - 2.*ro*a*b + b*b)
			  delta = (1.-Sign(a)*Sign(b)) / 4.
			  return,Norm2d3a(a,l,ro1) + Norm2d3a(b,l,ro2) - delta
	end
	endcase
end

function callmin,s1,s2,k,t,r,sig1,sig2,corr
sig=sqrt(sig1^2+sig2^2-2*sig1*sig2*corr)
gam1=alog(s1/k)+(r-sig1^2/2)/(sig1*sqrt(t))
gam2=alog(s2/k)+(r-sig2^2/2)/(sig2*sqrt(t))
alpha1=gam1+sig1*sqrt(t)
alpha2=(alog(s2/s1)-sig^2/2*sqrt(t))/(sig*sqrt(t))
roalpha=(corr*sig2-sig1)/sig
beta1=gam2+sig2*sqrt(t)
beta2=(alog(s1/s2)-sig^2/2*sqrt(t))/(sig*sqrt(t))
robeta=(corr*sig1-sig2)/sig
callmin=s1*Norm2d(alpha1,alpha2,roalpha)+$
	s2*norm2d(beta1,beta2,robeta)-k*exp(-r*t)*norm2d(gam1,gam2,corr)
return,callmin
end

function callmax,s1,s2,k,t,r,sig1,sig2,corr
callmin=callmin(s1,s2,k,t,r,sig1,sig2,corr)
return,call(s1,k,t,r,sig1)+call(s2,k,t,r,sig2)-callmin
end

function putmin,s1,s2,k,t,r,sig1,sig2,corr
callm=callmin(s1,s2,k,t,r,sig1,sig2,corr)
callzero=callmin(s1,s2,0.0,t,r,sig1,sig2,corr)
return,k*exp(-r*t)-callzero+callm
end

function putmax,s1,s2,k,t,r,sig1,sig2,corr
callm=callmax(s1,s2,k,t,r,sig1,sig2,corr)
callzero=callmax(s1,s2,0.0,t,r,sig1,sig2,corr)
return,k*exp(-r*t)-callzero+callm
end

function solveu2_u1,u1,M,V1,V2
x1=(double(1.) - 2*M + M^2 + u1^2 - 2*M*u1^2 + M^2*u1^2 + 2*u1*V1 + V2 -2*u1*V2 + u1^2*V2 - sqrt(1. - 2*M + M^2 + V2)*$
    sqrt(1 - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 +M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2) +$
    u1*sqrt(1. - 2*M + M^2 + V2)*sqrt(1 - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 + $
    M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2))/(2*u1*(1. - 2*M + M^2 + V1))

x2=(double(1.) - 2*M + M^2 + u1^2 - 2*M*u1^2 + M^2*u1^2 + 2*u1*V1 + V2 -2*u1*V2 + u1^2*V2 + sqrt(1. - 2*M + M^2 + V2)*$
    sqrt(1. - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 +M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2) -$
    u1*sqrt(1 - 2*M + M^2 + V2)*sqrt(1 - 2*M + M^2 + 2*u1 - 4*M*u1 + 2*M^2*u1 + u1^2 - 2*M*u1^2 +$ 
    M^2*u1^2 + 4*u1*V1 + V2 - 2*u1*V2 + u1^2*V2))/(2*u1*(1. - 2*M + M^2 + V1))
return,[x1,x2]
end


function ACRYLcallmax,s1,s2,k,t,tta,r,sigma1,sigma2,corr,n,yieldlist1,t_divlist1,yieldlist2,t_divlist2
lambda=1.1
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
call1=fltarr(2*n+1,2*n+1)
call2=fltarr(2*n+1,2*n+1)
supp1=fltarr(2*n+1)
supp2=fltarr(2*n+1)
M=(exp(r*double(deltat)))
V1=(exp(double(sigma1)^2*deltat)-1.)*M^2
V2=(exp(double(sigma2)^2*deltat)-1.)*M^2
RR=(exp(corr*double(sigma1)*double(sigma2)*double(deltat)))*M^2
u1=exp(lambda*double(sigma1)*sqrt(double(deltat)))
d1=1./u1
ku2=solveu2_u1(u1,M,V1,V2)
kku2=where(ku2 ge 1.)
u2=min(ku2(kku2))
d2=1./u2
f1=((V1+M^2-M)*u1-(M-1.))/((u1-1.)*(u1^2-1.))
f2=((V2+M^2-M)*u2-(M-1.))/((u2-1.)*(u2^2-1.))
g1=((u1^2)*(V1+M^2-M)-(u1^3)*(M-1.))/((u1-1.)*(u1^2-1.))
g2=((u2^2)*(V2+M^2-M)-(u2^3)*(M-1.))/((u2-1.)*(u2^2-1.))
p1=(u1*u2*(rr-1.)-f1*(u1^2-1.)-f2*(u2^2-1.)+(f2+g2)*(u1*u2-1.))/((u1^2-1)*(u2^2-1))
p2=f1-p1
p2p=(f1*(u1^2-1.)*u2^2+f2*(u2^2-1.)-(f2+g2)*(u1*u2-1.)-u1*u2*(rr-1.))/((u1^2-1.)*(u2^2-1.))
p3=g2-p2
p4=f2-p1
p5=1-(f1+g1)
;print,'prob=',p1,p2,p3,p4,p5
indexlist1=where(t_divlist1 lt t,count1)
indexlist2=where(t_divlist2 lt t,count2)
for i=0,2*n do begin
		supp1(i)=s1*u1^(i-n) 
		supp2(i)=s2*u2^(i-n)
endfor
pp1=1.
pp2=1.
if count1 ne 0 then begin
		ylist1=yieldlist1(indexlist1)
		for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
endif
if count2 ne 0 then begin
		ylist2=yieldlist2(indexlist2)
		for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
endif
epsilon=-1
for j1=-n,n do begin
	nn=(n-((1+epsilon)/2))
	for j2=-nn,nn,2 do begin
		call2(n+j1,n+j2)=max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>0
		endfor
	epsilon=-epsilon
endfor
for i=n-1,0,-1 do begin
	call1(0:2*n,0:2*n)=0
	indexlist1=where(t_divlist1 lt i*deltat,count1)
	indexlist2=where(t_divlist2 lt i*deltat,count2)
	pp1=1.
	pp2=1.
	if count1 ne 0 then begin
			ylist1=yieldlist1(indexlist1)
			for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
	endif
	if count2 ne 0 then begin
			ylist2=yieldlist2(indexlist2)
			for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
	endif
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
		  if i*deltat gt tta then begin
			call1(n+j1,n+j2)=max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>$
			1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		  endif else begin
			call1(n+j1,n+j2)=1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		  endelse
		endfor
		epsilon=-epsilon
	endfor
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call2(n+j1,n+j2)=call1(n+j1,n+j2)
		endfor
		epsilon=-epsilon
	endfor
endfor
return,call2(n,n)
end

; otype vaut 'call' ou 'put'
; stype vaut 'min' ou 'max' ou 'moyenne' ou spread (positif ou pas suivant k)
; r est un taux qui servira a calculer la probabilite de hausse p
; ctaux sera utilisee pour actualiser d un pas de temps sur l autre et doit etre obtenu 
function ACRYLbioption,otype,stype,s1,s2,k,t,tta,r,sigma1,sigma2,corr,n,yieldlist1,t_divlist1,yieldlist2,t_divlist2,ctaux=ctaux
lambda=1.1
if otype ne 'call' and otype ne 'put' then print,' erreur de otype:',otype
if stype ne 'max' and stype ne 'min' and stype ne 'moyenne' and stype ne 'spread' then print,' erreur de stype:',stype
deltat=max([float(t)/n,0])
if deltat eq 0 then return,max([s-k,0])
call1=fltarr(2*n+1,2*n+1)
call2=fltarr(2*n+1,2*n+1)
supp1=fltarr(2*n+1)
supp2=fltarr(2*n+1)
M=(exp(r*double(deltat)))
V1=(exp(double(sigma1)^2*deltat)-1.)*M^2
V2=(exp(double(sigma2)^2*deltat)-1.)*M^2
RR=(exp(corr*double(sigma1)*double(sigma2)*double(deltat)))*M^2
u1=exp(lambda*double(sigma1)*sqrt(double(deltat)))
d1=1./u1
ku2=solveu2_u1(u1,M,V1,V2)
kku2=where(ku2 ge 1.)
u2=min(ku2(kku2))
d2=1./u2
vtaux=replicate(r,n)
if keyword_set(ctaux) ne 0 then begin
for i=0,n-1 do vtaux(i)=interpolation(ctaux,deltat*i)/100.
endif
f1=((V1+M^2-M)*u1-(M-1.))/((u1-1.)*(u1^2-1.))
f2=((V2+M^2-M)*u2-(M-1.))/((u2-1.)*(u2^2-1.))
g1=((u1^2)*(V1+M^2-M)-(u1^3)*(M-1.))/((u1-1.)*(u1^2-1.))
g2=((u2^2)*(V2+M^2-M)-(u2^3)*(M-1.))/((u2-1.)*(u2^2-1.))
p1=(u1*u2*(rr-1.)-f1*(u1^2-1.)-f2*(u2^2-1.)+(f2+g2)*(u1*u2-1.))/((u1^2-1)*(u2^2-1))
p2=f1-p1
p2p=(f1*(u1^2-1.)*u2^2+f2*(u2^2-1.)-(f2+g2)*(u1*u2-1.)-u1*u2*(rr-1.))/((u1^2-1.)*(u2^2-1.))
p3=g2-p2
p4=f2-p1
p5=1-(f1+g1)
;print,'prob=',p1,p2,p3,p4,p5
indexlist1=where(t_divlist1 lt t,count1)
indexlist2=where(t_divlist2 lt t,count2)
for i=0,2*n do begin
		supp1(i)=s1*u1^(i-n) 
		supp2(i)=s2*u2^(i-n)
endfor
pp1=1.
pp2=1.
if count1 ne 0 then begin
		ylist1=yieldlist1(indexlist1)
		for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
endif
if count2 ne 0 then begin
		ylist2=yieldlist2(indexlist2)
		for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
endif
epsilon=-1
for j1=-n,n do begin
	nn=(n-((1+epsilon)/2))
	for j2=-nn,nn,2 do begin
case 1 of
	otype eq 'call' and stype eq 'max':call2(n+j1,n+j2)=max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>0
	otype eq 'call' and stype eq 'min':call2(n+j1,n+j2)=min([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>0
	otype eq 'call' and stype eq 'moyenne':call2(n+j1,n+j2)=(supp1(n+j1)*pp1+supp2(n+j2)*pp2)/2.-k>0
	otype eq 'call' and stype eq 'spread':call2(n+j1,n+j2)=(supp1(n+j1)*pp1-supp2(n+j2)*pp2)-k>0
	otype eq 'put' and stype eq 'max':call2(n+j1,n+j2)=k-max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])>0
	otype eq 'put' and stype eq 'min':call2(n+j1,n+j2)=k-min([supp1(n+j1)*pp1,supp2(n+j2)*pp2])>0
	otype eq 'put' and stype eq 'moyenne':call2(n+j1,n+j2)=k-(supp1(n+j1)*pp1+supp2(n+j2)*pp2)/2.>0
	otype eq 'put' and stype eq 'spread':call2(n+j1,n+j2)=k-(supp1(n+j1)*pp1-supp2(n+j2)*pp2)>0
endcase
		endfor
	epsilon=-epsilon
endfor
for i=n-1,0,-1 do begin
	call1(0:2*n,0:2*n)=0
	indexlist1=where(t_divlist1 lt i*deltat,count1)
	indexlist2=where(t_divlist2 lt i*deltat,count2)
	pp1=1.
	pp2=1.
	if count1 ne 0 then begin
			ylist1=yieldlist1(indexlist1)
			for j=0,count1-1 do pp1=pp1*(1-ylist1(j))
	endif
	if count2 ne 0 then begin
			ylist2=yieldlist2(indexlist2)
			for j=0,count2-1 do pp2=pp2*(1-ylist2(j))
	endif
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
		  if i*deltat gt tta then begin
			kappa= exp(-vtaux(i)*deltat)*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
			case 1 of 
				otype eq 'call' and stype eq 'max':call1(n+j1,n+j2)=$
					max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>kappa
				otype eq 'call' and stype eq 'min':call1(n+j1,n+j2)=$
					min([supp1(n+j1)*pp1,supp2(n+j2)*pp2])-k>kappa
				otype eq 'call' and stype eq 'moyenne':call1(n+j1,n+j2)=$
					(supp1(n+j1)*pp1+supp2(n+j2)*pp2)/2.-k>kap    pa
				otype eq 'call' and stype eq 'spread':call1(n+j1,n+j2)=$
					(supp1(n+j1)*pp1-supp2(n+j2)*pp2)-k>kappa
				otype eq 'put' and stype eq 'max':call1(n+j1,n+j2)=$
					k-max([supp1(n+j1)*pp1,supp2(n+j2)*pp2])>kappa
				otype eq 'put' and stype eq 'min':call1(n+j1,n+j2)=$
					k-min([supp1(n+j1)*pp1,supp2(n+j2)*pp2])>kappa
				otype eq 'put' and stype eq 'moyenne':call1(n+j1,n+j2)=$
					k-(supp1(n+j1)*pp1+supp2(n+j2)*pp2)/2.>kappa
				otype eq 'put' and stype eq 'spread':call1(n+j1,n+j2)=$
					k-(supp1(n+j1)*pp1-supp2(n+j2)*pp2)>kappa

			endcase
		  endif else begin
			call1(n+j1,n+j2)=1./M*(p1*call2(n+j1+1,n+j2+1)+p2*call2(n+j1+1,n+j2-1)+$
			p3*call2(n+j1-1,n+j2-1)+p4*call2(n+j1-1,n+j2+1)+p5*call2(n+j1,n+j2))
		  endelse
		endfor
		epsilon=-epsilon
	endfor
	epsilon=-1
	for j1=-i,+i do begin
	nn=(i-((1+epsilon)/2))
		for j2=-nn,nn,2 do begin
			call2(n+j1,n+j2)=call1(n+j1,n+j2)
		endfor
		epsilon=-epsilon
	endfor
endfor
return,call2(n,n)
end
  



;***************************************************************************************************
;
;
;              calcul de parametre implicite (volatilite ,taux,...)
;
;
;**************************************************************************************************


function impvol_call,prime,s,k,t,r,x1=x1,x2=x2,xacc=xacc,jmax=jmax,yieldlist=yieldlist,t_yieldlist=t_yieldlist,ctaux=ctaux,tta=tta,n=n
if prime eq 0 then return,0
if keyword_set(x1) eq 0 then x1=0.001
if keyword_set(x2) eq 0 then x2=0.3
if keyword_set(xacc) eq 0. then xacc=0.001
if keyword_set(jmax) eq 0. then jmax=40
case 1 of
	keyword_set(yieldlist) eq 0.  and keyword_set(ctaux ) eq 0. : yieldflag=0  
	keyword_set(yieldlist) eq 1.  and keyword_set(ctaux ) eq 0. : yieldflag=1 
	keyword_set(ctaux ) eq 1. : yieldflag=2  
endcase  
if keyword_set(yieldlist) eq 0. then begin
	yieldlist=[0.] 
	t_yieldlist=[t/2.]
endif
if keyword_set(tta) eq 0. then begin
		ttaflag=0 
		tta=0.
		endif else ttaflag=1
if keyword_set(n) eq 0 then n=100
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :begin
		fmid=call(s,k,t,r,x2)-prime
		f=call(s,k,t,r,x1)-prime
	end
	yieldflag eq 1 and ttaflag eq 0 :begin
		fmid=YLcall(s,k,t,r,x2,yieldlist,t_yieldlist)-prime
		f=YLcall(s,k,t,r,x1,yieldlist,t_yieldlist)-prime
	end
	yieldflag eq 2  :begin
		fmid= ACRYLcall1(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
		f=ACRYLcall1(s,k,t,tta,r,x1,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	end
	else:begin
		fmid=ACRYLcall(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist)-prime
		f=ACRYLcall(s,k,t,tta,r,x1,n,yieldlist,t_yieldlist)-prime
	end
endcase		
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :fmid=call(s,k,t,r,x2)-prime
	yieldflag eq 1 and ttaflag eq 0 :fmid=YLcall(s,k,t,r,x2,yieldlist,t_yieldlist)-prime
	yieldflag eq 2  :fmid= ACRYLcall1(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	else:fmid=ACRYLcall(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist)-prime
endcase
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
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :fmid=call(s,k,t,r,xmid)-prime
	yieldflag eq 1 and ttaflag eq 0 :fmid=YLcall(s,k,t,r,xmid,yieldlist,t_yieldlist)-prime
	yieldflag eq 2  :fmid= ACRYLcall1(s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	else:fmid=ACRYLcall(s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist)-prime
endcase
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end


function imptaux_call,prime,s,k,t,v,x1=x1,x2=x2,xacc=xacc,jmax=jmax,yieldlist=yieldlist,t_yieldlist=t_yieldlist,ctaux=ctaux,tta=tta,n=n
if prime eq 0 then return,0
if keyword_set(x1) eq 0 then x1=0.
if keyword_set(x2) eq 0 then x2=0.2
if keyword_set(xacc) eq 0. then xacc=0.0001
if keyword_set(jmax) eq 0. then jmax=40
if keyword_set(yieldlist) eq 0. then yieldflag=0 else yieldflag=1
if keyword_set(tta) eq 0. then begin
		ttaflag=0 
		tta=0
		endif else ttaflag=1
if keyword_set(n) eq 0 then n=100
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :begin
		fmid=call(s,k,t,x2,v)-prime
		f=call(s,k,t,x1,v)-prime
	end
	yieldflag eq 1 and ttaflag eq 0 :begin
		fmid=YLcall(s,k,t,x2,v,yieldlist,t_yieldlist)-prime
		f=YLcall(s,k,t,x1,v,yieldlist,t_yieldlist)-prime
	end
	else:begin
		fmid=ACRYLcall(s,k,t,tta,x2,v,n,yieldlist,t_yieldlist)-prime
		f=ACRYLcall(s,k,t,tta,x1,v,n,yieldlist,t_yieldlist)-prime
	end
endcase
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :fmid=call(s,k,t,x2,v)-prime
	yieldflag eq 1 and ttaflag eq 0 :fmid=YLcall(s,k,t,x2,v,yieldlist,t_yieldlist)-prime
	else:fmid=ACRYLcall(s,k,t,tta,x2,v,n,yieldlist,t_yieldlist)-prime
endcase
	id=id+1
endwhile
;print,'fmid=',fmid,' /f=',f
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
	case 1 of
		yieldflag eq 0 and ttaflag eq 0 :fmid=call(s,k,t,xmid,v)-prime
		yieldflag eq 1 and ttaflag eq 0 :fmid=YLcall(s,k,t,xmid,v,yieldlist,t_yieldlist)-prime
		else:fmid=ACRYLcall(s,k,t,tta,xmid,v,n,yieldlist,t_yieldlist)-prime
	endcase
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end


function impvol_put,prime,s,k,t,r,x1=x1,x2=x2,xacc=xacc,jmax=jmax,yieldlist=yieldlist,t_yieldlist=t_yieldlist,ctaux=ctaux,tta=tta,n=n
if prime eq 0 then return,0
if keyword_set(x1) eq 0 then x1=0.001
if keyword_set(x2) eq 0 then x2=0.3
if keyword_set(xacc) eq 0. then xacc=0.001
if keyword_set(jmax) eq 0. then jmax=40
if keyword_set(yieldlist) eq 0. then yieldlist=[0.] & t_yieldlist=[t/2.]
case 1 of
	keyword_set(yieldlist) eq 0.  and keyword_set(ctaux ) eq 0. : yieldflag=0  
	keyword_set(yieldlist) eq 1.  and keyword_set(ctaux ) eq 0. : yieldflag=1 
	keyword_set(ctaux ) eq 1. : yieldflag=2  
endcase  

if keyword_set(tta) eq 0. then begin
		ttaflag=0 
		tta=0
		endif else ttaflag=1
if keyword_set(n) eq 0 then n=100
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :begin
		fmid=put(s,k,t,r,x2)-prime
		f=put(s,k,t,r,x1)-prime
	end
	yieldflag eq 1 and ttaflag eq 0 :begin
		fmid=YLput(s,k,t,r,x2,yieldlist,t_yieldlist)-prime
		f=YLput(s,k,t,r,x1,yieldlist,t_yieldlist)-prime
	end
	yieldflag eq 2  :begin
		fmid= ACRYLput1(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
		f=ACRYLput1(s,k,t,tta,r,x1,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	end
	else:begin
		fmid=ACRYLput(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist)-prime
		f=ACRYLput(s,k,t,tta,r,x1,n,yieldlist,t_yieldlist)-prime
	end
endcase
		
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :fmid=put(s,k,t,r,x2)-prime
	yieldflag eq 1 and ttaflag eq 0 :fmid=YLput(s,k,t,r,x2,yieldlist,t_yieldlist)-prime
	yieldflag eq 2  :fmid= ACRYLput1(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	else:fmid=ACRYLput(s,k,t,tta,r,x2,n,yieldlist,t_yieldlist)-prime
endcase
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
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :fmid=put(s,k,t,r,xmid)-prime
	yieldflag eq 1 and ttaflag eq 0 :fmid=YLput(s,k,t,r,xmid,yieldlist,t_yieldlist)-prime
	yieldflag eq 2  :fmid= ACRYLput1(s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist,ctaux=ctaux)-prime
	else:fmid=ACRYLput(s,k,t,tta,r,xmid,n,yieldlist,t_yieldlist)-prime
endcase
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end


function imptaux_put,prime,s,k,t,v,x1=x1,x2=x2,xacc=xacc,jmax=jmax,yieldlist=yieldlist,t_yieldlist=t_yieldlist,tta=tta,n=n
if prime eq 0 then return,0
if keyword_set(x1) eq 0 then x1=0.
if keyword_set(x2) eq 0 then x2=0.2
if keyword_set(xacc) eq 0. then xacc=0.0001
if keyword_set(jmax) eq 0. then jmax=40
if keyword_set(yieldlist) eq 0. then yieldflag=0 else yieldflag=1
if keyword_set(tta) eq 0. then begin
		ttaflag=0 
		tta=0
		endif else ttaflag=1
if keyword_set(n) eq 0 then n=100
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :begin
		fmid=put(s,k,t,x2,v)-prime
		f=put(s,k,t,x1,v)-prime
	end
	yieldflag eq 1 and ttaflag eq 0 :begin
		fmid=YLput(s,k,t,x2,v,yieldlist,t_yieldlist)-prime
		f=YLput(s,k,t,x1,v,yieldlist,t_yieldlist)-prime
	end
	else:begin
		fmid=ACRYLput(s,k,t,tta,x2,v,n,yieldlist,t_yieldlist)-prime
		f=ACRYLput(s,k,t,tta,x1,v,n,yieldlist,t_yieldlist)-prime
	end
endcase
id=0
while f*fmid gt 0. and id lt 7 do begin
	x2=2*x2
case 1 of
	yieldflag eq 0 and ttaflag eq 0 :fmid=put(s,k,t,x2,v)-prime
	yieldflag eq 1 and ttaflag eq 0 :fmid=YLput(s,k,t,x2,v,yieldlist,t_yieldlist)-prime
	else:fmid=ACRYLput(s,k,t,tta,x2,v,n,yieldlist,t_yieldlist)-prime
endcase
	id=id+1
endwhile
;print,'fmid=',fmid,' /f=',f
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
	case 1 of
		yieldflag eq 0 and ttaflag eq 0 :fmid=put(s,k,t,xmid,v)-prime
		yieldflag eq 1 and ttaflag eq 0 :fmid=YLput(s,k,t,xmid,v,yieldlist,t_yieldlist)-prime
		else:fmid=ACRYLput(s,k,t,tta,xmid,v,n,yieldlist,t_yieldlist)-prime
	endcase
	if fmid le 0 then rtbis=xmid
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis
endwhile
	if abs(dx) lt xacc or fmid eq 0 then return,rtbis else print,'nb max d iteration atteint:',jmax
return,0
end


  
function impvoltauxcall,s,k,t,tta,n,yieldlist,t_yieldlist,ctaux1,cvol1
;calcule la vol et le taux convenant a une prise ne compte des courbe de vol
;et des courbe de taux ainsi que des yieldidendes et du caractere americain 
; de l option
callbut=ACRYLcall1(s,k,t,tta,0.,0.,n,yieldlist,t_yieldlist,ctaux=ctaux1,cvol=cvol1)
vol_operationnelle=interpolation(cvol1,t)/100.
taux_operationnel=imptaux_call(callbut,s,k,t,vol_operationnelle)
;print,'call en courbe vol/taux=',callbut,'  call simplifie=',call(s,k,t,taux_operationnel,vol_operationnelle)
return,[vol_operationnelle,taux_operationnel]
end

function impvoltauxput,s,k,t,tta,n,yieldlist,t_yieldlist,ctaux1,cvol1
;calcule la vol et le taux convenant a une prise ne compte des courbe de vol
;et des courbe de taux ainsi que des yieldidendes et du caractere americain 
; de l option
putbut=ACRYLput1(s,k,t,tta,0.,0.,n,yieldlist,t_yieldlist,ctaux=ctaux1,cvol=cvol1)
vol_operationnelle=interpolation(cvol1,t)/100.
taux_operationnel=imptaux_put(putbut,s,k,t,vol_operationnelle)
return,[vol_operationnelle,taux_operationnel]
end

function ACRYLimpvoltauxcall,s,k,t,tta,n,yieldlist,t_yieldlist,ctaux1,cvol1
;calcule la vol et le taux convenant a une prise ne compte des courbe de vol
;et des courbe de taux 
callbut=ACRYLcall1(s,k,t,tta,0.,0.,n,yieldlist,t_yieldlist,ctaux=ctaux1,cvol=cvol1)
vol_operationnelle=interpolation(cvol1,t)/100.
taux_operationnel=imptaux_call(callbut,s,k,t,vol_operationnelle,yieldlist=yieldlist,t_yieldlist=t_yieldlist,n=n,tta=tta)
return,[vol_operationnelle,taux_operationnel]
end

function ACRYLimpvoltauxput,s,k,t,tta,n,yieldlist,t_yieldlist,ctaux1,cvol1
;calcule la vol et le taux convenant a une prise ne compte des courbe de vol
;et des courbe de taux 
putbut=ACRYLput1(s,k,t,tta,0.,0.,n,yieldlist,t_yieldlist,ctaux=ctaux1,cvol=cvol1)
vol_operationnelle=interpolation(cvol1,t)/100.
taux_operationnel=imptaux_put(putbut,s,k,t,vol_operationnelle,yieldlist=yieldlist,t_yieldlist=t_yieldlist,n=n,tta=tta)
return,[vol_operationnelle,taux_operationnel]
end



;****************************************************************************************
;
;
;                  calcul des deltas
;
;
;
;****************************************************************************************
function call_delta,s,k,t,r,sigm
if (t le 0.)  then begin
	if (s gt k) then  return,1. else return,0
endif else begin
	if (sigm eq 0.0) then return,1.0
	d1=(alog(s/k)+(r+(sigm*sigm)/2)*t)/(sigm*sqrt(t))
	return,gaussint(d1)
endelse
end

function put_delta,s,k,t,r,sigm
return,call_delta(s,k,t,r,sigm)-1
end

function delta,type,s,k,t,r,sigm
case  type of
1:return,call_delta(s,k,t,r,sigm)
2:return,put_delta(s,k,t,r,sigm)
else:return,0
endcase
end

function YLcall_delta,s,k,t,r,sigma,yieldlist,t_divlist
  s1=s*1.005
  s2=s*0.995
  c1=YLcall(s1,k,t,r,sigma,yieldlist,t_divlist)
  c2=YLcall(s2,k,t,r,sigma,yieldlist,t_divlist)
return,100*(c1-c2)/s
end

function YLput_delta,s,k,t,r,sigma,yieldlist,t_divlist
  s1=s*1.005
  s2=s*0.995
  c1=YLput(s1,k,t,r,sigma,yieldlist,t_divlist)
  c2=YLput(s2,k,t,r,sigma,yieldlist,t_divlist)
return,100*(c1-c2)/s
end


function CRcall_delta,s,k,t,r,sigma,n
  s1=s*1.005
  s2=s*0.995
  c1=CRcall(s1,k,t,r,sigma,n)
  c2=CRcall(s2,k,t,r,sigma,n)
return,100*(c1-c2)/s
end

function CRput_delta,s,k,t,r,sigma,n
  s1=s*1.005
  s2=s*0.995
  c1=CRput(s1,k,t,r,sigma,n)
  c2=CRput(s2,k,t,r,sigma,n)
return,100*(c1-c2)/s
end

function CREcall_delta,s,k,t,r,sigma,n,palier
  s1=s*1.005
  s2=s*0.995
  c1=CREcall(s1,k,t,r,sigma,n,palier)
  c2=CREcall(s2,k,t,r,sigma,n,palier)
return,100*(c1-c2)/s
end

function CREput_delta,s,k,t,r,sigma,n,palier
  s1=s*1.005
  s2=s*0.995
  c1=CREput(s1,k,t,r,sigma,n,palier)
  c2=CREput(s2,k,t,r,sigma,n,palier)
return,100*(c1-c2)/s
end

function CRYLcall_delta,s,k,t,r,sigma,n,yieldlist,t_divlist
  s1=s*1.005
  s2=s*0.995
  c1=CRYLcall(s1,k,t,r,sigma,n,yieldlist,t_divlist)
  c2=CRYLcall(s2,k,t,r,sigma,n,yieldlist,t_divlist)
return,100*(c1-c2)/s
end

function CRYLput_delta,s,k,t,r,sigma,n,yieldlist,t_divlist
  s1=s*1.005
  s2=s*0.995
  c1=CRYLput(s1,k,t,r,sigma,n,yieldlist,t_divlist)
  c2=CRYLput(s2,k,t,r,sigma,n,yieldlist,t_divlist)
return,100*(c1-c2)/s
end

 
;***************************************************************************
;
;
;                      fonctions de synthese
;
;
;
;
;***************************************************************************

function calcul_option,opt_type,opt_caracteristique,support,s0,k,t,t0,displayflag,$
americain=americain,debut_americain=debut_americain,asiatique=asiatique, resets1= resets,plafond_explosif=plafond_explosif,$
plancher_explosif=plancher_explosif,niveau_rappel=niveau_rappel, paliers= paliers,debut_paliers=debut_paliers,$
lookback=lookback,n_cox=n_cox,moyenne=moyenne,garantie_change=garantie_change,spectre=spectre
taux=courbetaux(marche_devise(support))
taux_emprunt=taux_emprunt_titre(support)
taux=soustrait_taux(taux,taux_emprunt)
vol=courbevol(support)
div_yield=dividende_yield(support,s0)
if keyword_set(garantie_change) ne 0 then begin
	devise1=marche_devise(support)
	devise2=garantie_change
	if devise1 ne devise2 then corr=devise_correlation(support,devise2) else corr=0.
	volmoy1=interpolation(vol,t)/100.
	volmoy2=devise_volatilite(devise2)/100.
	cov=corr*volmoy1*volmoy2*100.
	taux2=courbetaux(devise2)
	delta_taux=interpolation(taux2,t)-interpolation(taux,t)+taux_emprunt
	print,'ecart de taux:',delta_taux
print,'covariance devise/marche qui sera retire des taux de financements option :',cov
	taux=soustrait_taux(taux,(cov-delta_taux))
endif
div_date=dividende_date(support,t0,s0)
if keyword_set(n_cox) eq 0 then n_cox=150 
if keyword_set(moyenne) eq 0 then moyenne=s0 
if keyword_set(asiatique) eq 1 then begin
	volmoy=interpolation(vol,t)/100.
if opt_type eq 'call'then vec=asiatique_call_param(s0,moyenne,k,t,asiatique,0.,volmoy,taux=taux,yieldlist=div_yield,t_yieldlist=div_date) else $
		vec=asiatique_put_param(s0,moyenne,k,t,asiatique,0.,volmoy,taux=taux,yieldlist=div_yield,t_yieldlist=div_date)
	s01=vec(0)
	k01=vec(1)
	vol1=(vec(2)/volmoy)*vol
endif else vol1=vol & s01=s0 & k01=k
if  keyword_set(americain) eq 0 and keyword_set(debut_americain) eq 0 then tta=t
if  keyword_set(americain) eq 1 and keyword_set(debut_americain) eq 0 then tta=0
if  keyword_set(debut_americain) eq 1 then tta=(debut_americain-t0)/365.25
if keyword_set(debut_paliers) eq 0 then ttl=0 else ttl=(debut-paliers-t0)/365.25
if opt_caracteristique eq 0 and opt_type eq 'call' then $
	val=ACRYLcall1(s01,k01,t,tta,0.,0.,n_cox,div_yield,div_date,ctaux=taux,cvol=vol1)
if opt_caracteristique eq 0 and opt_type eq 'put' then $
	val=ACRYLput1(s01,k01,t,tta,0.,0.,n_cox,div_yield,div_date,ctaux=taux,cvol=vol1)
if opt_caracteristique eq 1 and opt_type eq 'call' and tta eq t then begin
	 resets=( resets-t0)/365.25
	vt=ACRYLimpvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=NYLRScall(s01,k01,t,vt(1),vt(0),n_cox,div_yield,div_date, resets)
	endif
if opt_caracteristique eq 1 and opt_type eq 'put' and tta eq t then begin
	 resets=( resets-t0)/365.25
	vt=ACRYLimpvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=NYLRSput(s01,k01,t,vt(1),vt(0),n_cox,div_yield,div_date, resets)
	endif
if opt_caracteristique eq 1 and opt_type eq 'call' and tta ne t then begin
	 resets=( resets-t0)/365.25
	vt=ACRYLimpvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=NACRYLRScall(s01,k01,t,tta,vt(1),vt(0),int(n_cox/2.),int(n_cox/2.),div_yield,div_date, resets)
	endif
if opt_caracteristique eq 1 and opt_type eq 'put' and tta ne t then begin
	 resets=( resets-t0)/365.25
	vt=ACRYLimpvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=NACRYLRSput(s01,k01,t,tta,vt(1),vt(0),int(n_cox/2.),int(n_cox/2.),div_yield,div_date, resets)
	endif
if opt_caracteristique eq 2 and opt_type eq 'call' then begin
	vt=ACRYLimpvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=ACRYLEcall(s01,k01,t,tta,ttl,vt(1),vt(0),n_cox,div_yield,div_date,plafond_explosif+s01-s0)
	endif
if opt_caracteristique eq 2 and opt_type eq 'put' then begin
	vt=ACRYLimpvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=ACRYLEput(s01,k01,t,tta,ttl,vt(1),vt(0),n_cox,div_yield,div_date,plancher_explosif+s01-s0)
	endif
if opt_caracteristique eq 3 and opt_type eq 'call' then begin
	vt=ACRYLimpvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=NACRYLEcall(s01,k01,t,tta,ttl,vt(1),vt(0),n_cox,div_yield,div_date, paliers+s01-s0)
	endif
if opt_caracteristique eq 3 and opt_type eq 'put' then begin
	vt=ACRYLimpvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=NACRYLEput(s01,k01,t,tta,ttl,vt(1),vt(0),n_cox,div_yield,div_date, paliers+s01-s0)
	endif
if opt_caracteristique eq 4 and opt_type eq 'call' then begin
	vt=impvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=RAPcall(s01,k01,niveau_rappel+s01-s0,t,vt(1),vt(0))
	endif
if opt_caracteristique eq 4 and opt_type eq 'put' then begin
	vt=impvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=RAPput(s01,k01,niveau_rappel+s01-s0,t,vt(1),vt(0))
	endif
if opt_caracteristique eq 5 and opt_type eq 'call' then begin
	vt=impvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=CMINcall(s01,s01,t,vt(1),vt(0))
	endif
if opt_caracteristique eq 5 and opt_type eq 'put' then begin
	vt=impvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=CMAXput(s01,s01,t,vt(1),vt(0))
	endif
if opt_caracteristique eq 6 and opt_type eq 'call' then begin
	vt=ACRYLimpvoltauxcall(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=SPACRYLEcall(s01,k01,t,tta,ttl,vt(1),vt(0),n_cox,div_yield,div_date, paliers+s01-s0,spectre)
	endif
if opt_caracteristique eq 6 and opt_type eq 'put' then begin
	vt=ACRYLimpvoltauxput(s01,k01,t,tta,n_cox,div_yield,div_date,taux,vol1)
	val=SPACRYLEput(s01,k01,t,tta,ttl,vt(1),vt(0),n_cox,div_yield,div_date, paliers+s01-s0,spectre)
	endif

valeur=val
return,valeur
end

pro calcul_call,support,s0=s0,k=k,top=top,t01=t0 ,resultat=resultat,$
americain=americain,debut_americain=debut_americain,asiatique=asiatique,resets=resets,plafond_explosif=plafond_explosif,$
plancher_explosif=plancher_explosif,niveau_rappel=niveau_rappel,paliers=paliers,debut_paliers=debut_paliers,$
lookback=lookback,n_cox=n_cox,moyenne=moyenne,garantie_change=garantie_change,spectre=spectre
displayflag=1
if keyword_set(t0) eq 0 then t0=aujourdhui() else t0=t01
if keyword_set(n_cox) eq 0 then n_cox=150 
t=(top-t0)/365.25
opt_type='call'
opt_determine=0
opt_caracteristique=0
if top lt 10193. then begin
 print,'erreur dans l echeance:',top
 opt_determine=2
endif
if keyword_set( resets) ne 0 then  begin
	opt_caracteristique=1 
	opt_determine=opt_determine+1 
endif
if keyword_set(plafond_explosif) ne 0 then  begin
	opt_caracteristique=2 
	opt_determine=opt_determine+1 
if keyword_set(plancher_explosif) ne 0 then opt_determine=2
endif
if keyword_set(  paliers) ne 0 then  begin
	opt_caracteristique=3
	opt_determine=opt_determine+1
endif
if keyword_set(  spectre) ne 0 then  begin
	opt_caracteristique=6
	if keyword_set(paliers) eq 0 then begin
		opt_determine=2
		print,'Il manque les niveaux du spectre : le mot-cle Paliers=[...]'
	endif
endif
if keyword_set(niveau_rappel) ne 0 then  begin
	opt_caracteristique=4
	opt_determine=opt_determine+1
endif
if keyword_set(lookback) ne 0 then begin
	opt_caracteristique=5
	opt_determine=opt_determine+1
	if lookback ne 'min' then opt_determine=2
endif
if opt_determine gt 1 then print, 'calcul impossible' else begin
val=calcul_option(opt_type,opt_caracteristique,support,s0,k,t,t0,displayflag,$
americain=americain,debut_americain=debut_americain,asiatique=asiatique,resets=resets,plafond_explosif=plafond_explosif,$
plancher_explosif=plancher_explosif,niveau_rappel=niveau_rappel,paliers=paliers,lookback=lookback,n_cox=n_cox,moyenne=moyenne,$
garantie_change=garantie_change,spectre=spectre)
pour=val/s0*100.
print,'valeur du call:',val,'    soit en pourcentage:',pour,'%'
if keyword_set(resultat) ne 0 then begin
	resultat=val
	print,'resultat dans resultat'
end
endelse
return
end

pro calcul_put,support,s0=s0,k=k,top=top,t01=t0 ,$
americain=americain,debut_americain=debut_americain,asiatique=asiatique, resets= resets,plafond_explosif=plafond_explosif,$
plancher_explosif=plancher_explosif,niveau_rappel=niveau_rappel,paliers=paliers,debut_paliers=debut_paliers,$
lookback=lookback,n_cox=n_cox,moyenne=moyenne,garantie_change=garantie_change,spectre=spectre,resultat=resultat
displayflag=1
if keyword_set(t0) eq 0 then t0=aujourdhui() else t0=t01
if keyword_set(n_cox) eq 0 then n_cox=150 
t=(top-t0)/365.25
opt_type='put'
opt_determine=0 &opt_caracteristique=0
if top lt 10193. then begin
 print,'erreur dans l echeance:',top
 opt_determine=2
endif
if keyword_set( resets) ne 0 then  begin
	opt_caracteristique=1 
	opt_determine=opt_determine+1
endif
if keyword_set(plancher_explosif) ne 0 then  begin
	opt_caracteristique=2 
	opt_determine=opt_determine+1
if keyword_set(plafond_explosif) ne 0 then opt_determine=2
endif
if keyword_set( paliers) ne 0 then  begin
	opt_caracteristique=3
	opt_determine=opt_determine+1
endif
if keyword_set(  spectre) ne 0 then  begin
	opt_caracteristique=6
	if keyword_set(paliers) eq 0 then begin
		opt_determine=2
		print,'Il manque les niveaux du spectre : le mot-cle Paliers=[...]'
	endif
endif
if keyword_set(niveau_rappel) ne 0 then  begin
	opt_caracteristique=4
	opt_determine=opt_determine+1
endif
if keyword_set(lookback) ne 0 then begin
	 opt_caracteristique=5
	opt_determine=opt_determine+1
	if lookback ne 'max' then opt_determine=2
endif
if opt_determine gt 1 then print, 'calcul impossible' else begin
val=calcul_option(opt_type,opt_caracteristique,support,s0,k,t,t0,displayflag,$
americain=americain,debut_americain=debut_americain,asiatique=asiatique, resets= resets,plafond_explosif=plafond_explosif,$
plancher_explosif=plancher_explosif,niveau_rappel=niveau_rappel, paliers= paliers,lookback=lookback,n_cox=n_cox,moyenne=moyenne,$
garantie_change=garantie_change,spectre=spectre)
pour=val/s0*100.
print,'valeur du put:',val,'    soit en pourcentage:',pour,'%'
if keyword_set(resultat) ne 0 then begin
	resultat=val
	print,'resultat dans resultat'
end
endelse
return
end




pro calcul_bioption,otype,stype,support1,support2,s01=s01,s02=s02,k=k,top=top,t01=t0 ,$
americain=americain,debut_americain=debut_americain,correlation=correlation,n_cox=n_cox
; ce qui suppose que les indices sur lesquelles on calcule l option sont compte dans cette monnaie
; de meme que la matrice de covariance et les dividendes 
; otype vaut 'call' ou 'put'
; stype vaut 'min' ou 'max' ou 'moyenne'
; on suppose le financement dans la monnaie associee au support1
displayflag=1
if keyword_set(t0) eq 0 then t0=aujourdhui() else t0=t01
if keyword_set(correlation) eq 0 then corr=0. else corr=correlation
if keyword_set(n_cox) eq 0 then n_cox=150 
if keyword_set(financement) eq 0 then financement='frf'
t=(top-t0)/365.25
vol1=courbevol(support1)
div_yield1=dividende_yield(support1,s01)
div_date1=dividende_date(support1,t0,s01)
vol2=courbevol(support2)
div_yield2=dividende_yield(support2,s02)
div_date2=dividende_date(support2,t0,s02)
ctaux=courbetaux(marche_devise(support1))
opt_determine=0 &opt_caracteristique=0
if top lt 10193. then begin
 print,'erreur dans l echeance:',top
 opt_determine=2
endif
if  keyword_set(americain) eq 0 and keyword_set(debut_americain) eq 0 then tta=t
if  keyword_set(americain) eq 1 and keyword_set(debut_americain) eq 0 then tta=0
if  keyword_set(debut_americain) eq 1 then tta=(debut_americain-t0)/365.25
if opt_determine gt 1 then print, 'calcul impossible' else begin
	if otype eq 'call' then vt1=ACRYLimpvoltauxcall(s01,s01*k,t,tta,n_cox,div_yield1,div_date1,ctaux,vol1) $
  	else vt1=ACRYLimpvoltauxput(s01,s01*k,t,tta,n_cox,div_yield1,div_date1,ctaux,vol1)
	r=vt1(1)
	sigma1=vt1(0)
if otype eq 'call' then begin
	callbut=ACRYLcall1(s02,s02*k,t,tta,0.,0.,n_cox,div_yield2,div_date2,ctaux=ctaux,cvol=vol2)
	sigma2=impvol_call(callbut,s02,s02*k,t,r,n=n_cox,yieldlist=div_yield2,t_yieldlist=div_date2,tta=tta)
endif else begin
	putbut=ACRYLput1(s02,s02*k,t,tta,0.,0.,n_cox,div_yield2,div_date2,ctaux=ctaux,cvol=vol2)
	sigma2=impvol_put(putbut,s02,s02*k,t,r,n=n_cox,yieldlist=div_yield2,t_yieldlist=div_date2,tta=tta)
endelse
;	print,'r1=',r,'   sigma1=',sigma1,'    sigma2=',sigma2
 	if keyword_set(asiatique) eq 1 then begin
	sigma11=sigma1
	s011=s01
	sigma22=sigma2
	s022=s02
 	endif else begin
	sigma11=sigma1
	s011=s01
	sigma22=sigma2
	s022=s02
	endelse
	if  keyword_set(americain) eq 0 and keyword_set(debut_americain) eq 0 then tta=t
	if  keyword_set(americain) eq 1 and keyword_set(debut_americain) eq 0 then tta=0
	if  keyword_set(debut_americain) eq 1 then tta=(debut_americain-t0)/365.25
	if keyword_set( resets) ne 0 then  begin
	opt_caracteristique=1 
	opt_determine=opt_determine+1
	endif
endelse
if opt_determine gt 1 then print, 'calcul impossible' else begin
n_cox2=int(3./2.*float(n_cox)^(2./3.))
s011=1.0
s022=1.0
val=ACRYLbioption(otype,stype,s011,s022,k,t,tta,r,sigma11,sigma22,corr,n_cox2,div_yield1,div_date1,div_yield2,div_date2,ctaux=ctaux)
pour=val*100.
print,'valeur du ',otype,' ',stype,': en pourcentage de supports:',pour,'%'
endelse
end



;***************************************************************************
;
;
;                      Tests des fonctions
;
;
;
;
;***************************************************************************

  
;representation des dividendes, pondere par les capitalisations exprimees en points d index
; pevu pour 92 en date du 12/2/92 par Capel S.A.
pro rep_divid
divlist=[0.32,0.16,0.1,0.45,0.79,0.37,0.42,0.59,0.41,0.40,1.55,1.89,0.69,0.75,2.1, $
1.45,0.12,0.32,1.78,0.34,1.90,0.8,0.78,1.16,0.89,0.93,2.95,2.05,0.1,6.06,1.98,0.17, $
1.32,1.45,0.37,0.11,1.37,1.79,0.49,0.42,0.30,0.36,0.46]
datelist=[0116,130,203,326,425,510,521,603,603,604,606,606,607,610,610,613,614,618,618,$
618,624,628,628,701,701,701,701,701,701,702,702,703,703,704,704,709,710,715,716,725,726,$
731]
n=n_elements(datelist)
newjj=fltarr(n)
for i=0,n-1 do begin
	m=int(datelist(i)/100.)
	d=datelist(i)-m*100
	newjj(i)=(jj(d*10000.+m*100.+92)-jj(010192))/365.
	endfor
window,0,title='dividendes cac40 (en points d indice) / au cours de l annee'
plot,newjj,divlist,psym=4,xrange=[0,1.]
intdiv=fltarr(n)
intdiv(0)=divlist(0)
for i=1,n-1 do begin
intdiv(i)=intdiv(i-1)+divlist(i)
endfor
window,1,title='somme des versements de dividendes'
plot,newjj,intdiv,xrange=[0,1.]
s=0
for i=0,n-1 do s=s+newjj(i)*divlist(i)
s=s/intdiv(n-1)
print,'distance moyenne du dividende annuel :'+string(s)+' annees
d=caldat(jj(010192)+int(s*365.))
print,'soit en 92 le '+string(d)
end
 
