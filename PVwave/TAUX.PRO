;debut

 
function volforward,courb_vol
s=size(courb_vol)
n=s(1)
res=fltarr(n,3)
res(0,0)=0
res(0,1)=courb_vol(0,0)
res(0,2)=courb_vol(0,1)
for i=1,n-1 do begin
res(i,0)=courb_vol(i-1,0)
res(i,1)=courb_vol(i,0)-courb_vol(i-1,0)
vol=(courb_vol(i,1)^2*courb_vol(i,0)-courb_vol(i-1,1)^2*courb_vol(i-1,0))/(courb_vol(i,0)-courb_vol(i-1,0))
if vol lt 0 then begin
vol=0
print,'il y a erreur dans la vol',i+1
endif
res(i,2)=sqrt(vol)
endfor
return,res
end;


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
     somme=((1+courb_taux(i,1))^(courb_taux(j,0))-1)/(1+res(j,1))^(courb_taux(j,0))
   endif else begin
     somme=somme+((1+courb_taux(i,1))^(courb_taux(j,0)-courb_taux(j-1,0))-1)/(1+res(j,1))^(courb_taux(j,0))
   endelse 
 endfor
 res(i,1)=((1+courb_taux(i,1))^(courb_taux(i,0)-courb_taux(i-1,0))/(1-somme))^(1/courb_taux(i,0))-1.
endfor
return,res
end;

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


pro printcourbetaux, courb_taux,pr
; pr = 1 --> print + ecran
; pr = 2 --> ecran
; pr > 2 --> rien
ss=size(courb_taux)
n=ss(1)
if pr eq 1  then  begin
	set_plot,'PS'
	plot,courb_taux(*,0),courb_taux(*,1),title='zero coupons'
;	device,/landscape
	device,/close_file
	spawn,'lpr wave.ps'
endif
set_plot,'X'
if pr le 2 then  begin
		window,0,title='zero coupons'
	plot,courb_taux(*,0),courb_taux(*,1),title='zero coupons'
endif
end

function transformat_taux,courbe
scourbe=size(courbe)
nb=scourbe(2)
res=fltarr(nb,2)
for i=0,nb-1 do begin
	res(i,0)=courbe(0,i)
	res(i,1)=courbe(1,i)/100.
endfor
return,res
end

function ech,str_mois,otype
case 1 of
otype eq 11 or otype eq 12 :case 1 of
	str_mois eq 'mar3' : a=jj(310393)
	str_mois eq 'jun3' : a=jj(300693)
	str_mois eq 'sep3' : a=jj(310993)
	str_mois eq 'dec3' : a=jj(311293)
	else : print,'erreur de specification d echeance:',str_mois
	endcase
else : print,'erreur de specification de type pour l echeance:',otype
endcase
return,a
end

function transformat_prix,courbe,t0
scourbe=size(courbe)
if scourbe(0) eq 1 then nb=1 else nb=scourbe(2)
res=fltarr(nb,4)
for i=0,nb-1 do begin
	res(i,0)=courbe(0,i)			; le type d option
	res(i,1)=(courbe(1,i)-t0)/365.			; temps restant a courir jusque a l echeance
	res(i,2)=1./(courbe(2,i)/100.)-1.	; taux annuel sous entendu
	res(i,3)=(1+double(courbe(3,i)/100.))^((courbe(1,i)-t0)/365.)-1. ; calcul de la prime
endfor
return,res
end


function tauxjjforward,courb_taux
s=size(courb_taux)
n=s(1)
res=fltarr(n,3)
coupons=couponszero(courb_taux)
res(0,0)=0        ;t01
res(0,1)=courb_taux(0,0)
res(0,2)=coupons(0,1)
for i=1,n-1 do begin
 res(i,0)=courb_taux(i-1,0)
 res(i,1)=courb_taux(i,0)-courb_taux(i-1,0)
 a=(1.+coupons(i,1))^courb_taux(i,0)
 b=(1.+coupons(i-1,1))^courb_taux(i-1,0)
 res(i,2)=(a/b)^(1/(courb_taux(i,0)-courb_taux(i-1,0)))-1
endfor
return,res
end;


;*******************************************************************************
;
;
;
;        construction de la courbe de taux selon Augros
;
;
;
;*******************************************************************************

; le record reuter des option pibor3mois mars93 est  Q0#PIBH3+
; le record reuter des option pibor3mois juin93 est  Q0#PIBM3+
; le record reuter des option pibor3mois sept93 est  Q0#PIBV3+
; le record reuter des option pibor3mois dec93  est  Q0#PIBZ3+


pro visubin,bin,p
ss=size(bin)
nn=ss(1)
for i=0,(nn-1)<p do begin
	ff=fltarr(ss(2),ss(3))
	ff(0:ss(2)-1,0:ss(3)-1)=bin(i,0:ss(2)-1,0:ss(3)-1)
	print,'i=',transpose(ff)
endfor
end

pro displaybin,bin,p=p,w=w
ss=size(bin)
nd=ss(0)
if keyword_set(w) then w1=w else w1=0
if keyword_set(p) then p1=p else p1=0
if nd eq 3 then begin
	ff=fltarr(ss(2),ss(3))
	ff(0:ss(2)-1,0:ss(3)-1)=bin(p1,0:ss(2)-1,0:ss(3)-1)
	ndd=ss(2)
endif
if nd eq 2 then begin
	ff=fltarr(ss(1),ss(2))
	ff(0:ss(1)-1,0:ss(2)-1)=bin(0:ss(1)-1,0:ss(2)-1)
	ndd=ss(1)
endif
super=fltarr(ndd)
infer=fltarr(ndd)
ulper=fltarr(ndd)
for i=0,ss(2)-1 do begin
	infer(i)=ff(i,0)
	super(i)=ff(i,i)
endfor
ymax=max(super)
ymin=min(infer)
epsilon=ymax-ymin
window,w1,title='arbre binomial '+string(w1)
plot,super,yrange=[ymin-epsilon/2.,ymax+epsilon/10.]
oplot,infer
end


function format_coupons,duree,coupon,delai_prem_coupon,method_flag
; le but de cette fonction est de fabriquer un vecteur de dividendes de la forme :
;  [ ... [delai , coupon] ...]
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle

if duree le delai_prem_coupon then begin
		vect_divid=[[delai_prem_coupon,coupon]] 
endif else begin
		vect_divid_d=fltarr(500)
		vect_divid_t=fltarr(500)
		ndiv=0
		delta_t_rest=duree-delai_prem_coupon
		while delta_t_rest ge 0. do begin
			vect_divid_t(ndiv)=delai_prem_coupon+ndiv
			vect_divid_d(ndiv)=coupon
			delta_t_rest=delta_t_rest-1.
			ndiv=ndiv+1.
		endwhile
		delta_t_rest=delta_t_rest+1.
		if delta_t_rest gt 0.00001*duree then begin
			if method_flag eq 0 then cp=(1.+coupon)^delta_t_rest-1. $
				else cp=coupon*delta_t_rest
			vect_divid_t(ndiv)=duree
			vect_divid_d(ndiv)=cp
			ndiv=ndiv+1
		endif
			vect_divid=fltarr(ndiv,2)
		vect_divid(*,0)=vect_divid_t(0:ndiv-1)
		vect_divid(*,1)=vect_divid_d(0:ndiv-1)
endelse
return,vect_divid
end

function calcule_p_delta,t,n,vol,q,r
tstep=t/double(n)
p=((1.+r)^tstep-exp(-vol*sqrt(tstep)))/(exp(vol*sqrt(tstep))-exp(-vol*sqrt(tstep)))
delta=exp(-vol*sqrt(tstep)*q*(1-q)/double(n))
;print,'p=',p,' / 1000X(1- delta)=',1000.*(1.-delta)
return,[p,delta]
end



function binom_taux1,G,p,delta,t
;G est une liste dont chaque element est un yield 
; comprennant n+1 element de la sorte dont le premier n est pas utilise
; t est un numero compris entre 0 et n-1 qui reference un temps
;p et delta sont les parametre du modele
n=n_elements(g)
gni=dblarr(n,n)
gni(0,0)=G(t-1)
if t lt n then begin
	h=1/(p+(1.-p)*delta^dindgen(n))
	z=dblarr(n-t+1)
	z(1)=h(t)
	for i=2,n-t do z(i)=z(i-1)*h(t+i-1)/h(i-1)
	for ni=1,n-t do begin
		gni(ni,0)=G(t+ni-1)/G(ni-1)*z(ni)*delta^(t*ni)
		for i=1,ni do gni(ni,i)=gni(ni,0)/(delta^(t*i))
	endfor
endif
return,gni
end


function binom_taux2,tauxzero,p,delta,duration,t1,tstep
;fonction qui retourne un reseau des valeurs d un coupon zero d echeance duration ,
;  le reseau allant de 0 a t1 avec n etapes
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duree pour laquelle on calcule les differentes valeurs du coupon zero dans les differents
;	etats du monde
; t1 est l echeance maximale du treillis
; tstep est la maille temporelle du treillis
; p et delta sont les parametres du modele
;retourne deux arbres : un pour la duration demandee et l autre pour la duration de 1 pas de reseau (tstep)
ss=size(tauxzero)
nn=ss(1)
dt=tauxzero(*,0)
cp1=tauxzero(*,1)
cp=dblarr(nn)
for i=0,nn-1 do cp(i)=1./(1.+double(cp1(i)))^dt(i)
n=ceilling((double(t1+duration)/tstep))+2
n1=ceilling(double(t1)/tstep)+1
tmps=(1.+dindgen(n))*tstep
cps=dblarr(n)
for i=0,n-1 do begin
	if tmps(i) le dt(0) then cps(i)=(1.+cp1(0))^(-tmps(i))
	;if tmps(i) ge dt(nn-1) then cps(i)=cp(nn-1)
	if tmps(i) ge dt(nn-1) then cps(i)=(1.+cp1(nn-1))^(-tmps(i))
	if (tmps(i) gt dt(0)) and (tmps(i) lt dt(nn-1)) then begin
		indless=where(dt lt tmps(i))
		indt1=max(indless)
		indsup=where(dt ge tmps(i))
		indt2=min(indsup)
		;cps(i)=cp(indt1)+((cp(indt2)-cp(indt1))/(dt(indt2)-dt(indt1)))*(tmps(i)-dt(indt1))
		tt1=double(dt(indt1))
		tt2=double(dt(indt2))
				; on interpolle les taux
		cps(i)=(1.+(cp1(indt1)*(1.-(tmps(i)-tt1)/(tt2-tt1))+cp1(indt2)*((tmps(i)-tt1)/(tt2-tt1))))^(-tmps(i))  

	endif
endfor
;window,0,title='taux coupon zero'
;plot,dt,cp1
;window,10,title='actualisation coupons zero'
;plot,dt,cp
;window,11,title='actualisation interpolee'
;plot,tmps,cps
nt1=int(duration/tstep)>1
	bin1=binom_taux1(cps,p,delta,nt1)
	if (nt1 lt n) and (duration/tstep gt 1.) then begin 
		bin2=binom_taux1(cps,p,delta,nt1+1)
		bin=bin1+((duration/tstep)-nt1)*(bin2-bin1)
	endif else bin=bin1
bin1=binom_taux1(cps,p,delta,1)
result=dblarr(2,n1+1,n1+1)
result(0,*,*)=bin(0:n1,0:n1)
result(1,*,*)=bin1(0:n1,0:n1)
return,result
end


function binom_tauxn,tauxzero,p,delta,vect_divid,t1,tstep
;fonction qui retourne un ensemble de reseau des valeurs d un coupon zero d echeance variable ,
;-----------------------------------------------------------------------I
;  le reseau allant de 0 a t1 avec n1 etapes				I
; a chaque coupons [ti,ci] corresponds un arbre de duration associe ti  I
;-----------------------------------------------------------------------I
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duree pour laquelle on calcule les differentes valeurs du coupon zero dans les differents
;	etats du monde
; t1 est l echeance maximale du treillis
; n1 est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
;retourne autant d arbres qu il y a de coupons plus un arbre pour le taux court
ss=size(tauxzero)
nn=ss(1)
dt=tauxzero(*,0)
cp1=tauxzero(*,1)
cp=dblarr(nn)
for i=0,nn-1 do cp(i)=1./(1.+cp1(i))^dt(i)
vect_duration=vect_divid(*,0)
duration=max(vect_duration)
n_vect=n_elements(vect_duration)
t=double(t1+duration)
n=ceilling((double(t1+duration)/tstep))+2
n1=ceilling(double(t1)/tstep)+1
tmps=(1.+dindgen(n))*tstep
cps=dblarr(n)
for i=0,n-1 do begin
	if tmps(i) le dt(0) then cps(i)=(1.+cp1(0))^(-tmps(i))
	if tmps(i) ge dt(nn-1) then cps(i)=(1.+cp1(nn-1))^(-tmps(i))
	if (tmps(i) gt dt(0)) and (tmps(i) lt dt(nn-1)) then begin
		indless=where(dt lt tmps(i))
		indt1=max(indless)
		indsup=where(dt ge tmps(i))
		indt2=min(indsup)
		tt1=double(dt(indt1))
		tt2=double(dt(indt2))
			 ; on interpolle les taux
		cps(i)=(1.+(cp1(indt1)*(1.-(tmps(i)-tt1)/(tt2-tt1))+cp1(indt2)*((tmps(i)-tt1)/(tt2-tt1))))^(-tmps(i)) 
	endif
endfor
result=dblarr(n_vect+1,n1+1,n1+1)
for i=0,n_vect-1 do begin
	nt1=int(vect_duration(i)/tstep)>1
	bin1=binom_taux1(cps,p,delta,nt1)
	if (nt1 lt n) and (vect_duration(i)/tstep gt 1.) then begin 
		bin2=binom_taux1(cps,p,delta,nt1+1)
		bin=bin1+((vect_duration(i)/tstep)-nt1)*(bin2-bin1)
	endif else bin=bin1
	result(i,*,*)=bin(0:n1,0:n1)
endfor
bin1=binom_taux1(cps,p,delta,1)
result(n_vect,*,*)=bin1(0:n1,0:n1)
return,result
end




function binom_tauxnk,tauxzero,p,delta,durationmax,nmax,t1,tstep
;fonction qui retourne un reseau des valeurs des coupon zero d echeance variable entre 0 et tmax par pas de tmax/(n1/nmax) ,
;-------------------------------------------------------------------------------I
;  le reseau allant de 0 a t1 avec n1 etapes 				        I
; les arbres ont une duration qui va de 0 a max{t1,tmax} par pas de nmax*t1/n1  I
;-------------------------------------------------------------------------------I
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; t1 est l echeance maximale du treillis
; n1 est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; tmax est la duration max envisagee
; nmax est le nombre d etape elementaire constituant le pas de temps separant deux arbres successifs de duration
ss=size(tauxzero)
nn=ss(1)
dt=tauxzero(*,0)
cp1=tauxzero(*,1)
cp=dblarr(nn)
for i=0,nn-1 do cp(i)=1./(1.+cp1(i))^dt(i)
n1=ceilling(double(t1)/tstep)+1
n2=ceilling(double(durationmax)/tstep)+1
n_arbre=ceilling(n2/double(nmax))+1
n=n1+(n_arbre*nmax)+1
tmps=(1.+dindgen(n))*tstep
cps=dblarr(n)
for i=0,n-1 do begin
	if tmps(i) le dt(0) then cps(i)=(1.+cp1(0))^(-tmps(i))
	if tmps(i) ge dt(nn-1) then cps(i)=(1.+cp1(nn-1))^(-tmps(i))
	if (tmps(i) gt dt(0)) and (tmps(i) lt dt(nn-1)) then begin
		indless=where(dt lt tmps(i))
		indt1=max(indless)
		indsup=where(dt ge tmps(i))
		indt2=min(indsup)
		tt1=double(dt(indt1))
		tt2=double(dt(indt2))
 			; on interpolle les taux
		cps(i)=(1.+(cp1(indt1)*(1.-(tmps(i)-tt1)/(tt2-tt1))+cp1(indt2)*((tmps(i)-tt1)/(tt2-tt1))))^(-tmps(i)) 
	endif
endfor
result=dblarr(n_arbre+1,n1+1,n1+1)
for i=0,n_arbre do begin
	nt=(i*nmax)>1
	bin=binom_taux1(cps,p,delta,nt)
	result(i,*,*)=bin(0:n1,0:n1)
endfor
return,result
end



function taux_call,k,n,tstep,duration,tauxzero,p,delta,method_flag
;option sur FRA garantissant un taux k et payant en precompte
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le taux-prix d exercice; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis jusqu a l echeance
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=tstep*n
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
;displaybin,bin,p=0,w=1
;displaybin,bin,p=1,w=2
if method_flag eq 0 then gmax=(1+k)^duration else gmax=1+k*duration
tree_call=dblarr(n+1,n+1)
; fixation de la valeur a l echeance de l option
for i=0,n do begin
	if method_flag eq 0 then tree_call(n,i)=(1./bin(0,n,i)-gmax)*bin(0,n,i) > 0. $         ;methode actuarielle
		else begin
			tm=(1./bin(0,n,i)-1.)/t       			 ; methode interbancaire
			tree_call(n,i)=((tm-k)*t/(1.+tm*t))>0. 	 	 ; methode interbancaire
		endelse
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	for i=0,ni do begin
		tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
	endfor
endfor
;displaybin,tree_call
return,tree_call(0,0)
end


function taux_put,k,n,tstep,duration,tauxzero,p,delta,method_flag
;option sur FRA garantissant un taux k et payant en precompte
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le taux-prix d exercice; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=n*tstep
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
if method_flag eq 0 then gmax=(1+k)^duration else gmax=1+k*duration
tree_call=dblarr(n+1,n+1)
; fixation de la valeur a l echeance de l option
for i=0,n do begin
	if method_flag eq 0 then tree_call(n,i)=(gmax-1./bin(0,n,i))*bin(0,n,i) > 0. $         ;methode actuarielle
		else begin
			tm=(1./bin(0,n,i)-1.)/t       			 ; methode interbancaire
			tree_call(n,i)=((k-tm)*t/(1.+tm*t))>0. 	 	 ; methode interbancaire
		endelse
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	for i=0,ni do begin
		tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
	endfor
endfor
return,tree_call(0,0)
end


function taux_Acall,k,n,na,tstep,duration,tauxzero,p,delta,method_flag
;option sur FRA garantissant un taux k et payant en precompte
;ammericaine a partir de tta jusqu a l echeance
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le taux-prix d exercice; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis jusqu a l echeance
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=n*tstep
tta=na*tstep
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
if method_flag eq 0 then gmax=(1+k)^duration else gmax=1+k*duration
tree_call=dblarr(n+1,n+1)
; fixation de la valeur a l echeance de l option
for i=0,n do begin
	if method_flag eq 0 then tree_call(n,i)=(1./bin(0,n,i)-gmax)*bin(0,n,i) > 0. $	 		;methode actuarielle
		else begin
			tm=(1./bin(0,n,i)-1.)/t       			 				; methode interbancaire
			tree_call(n,i)=((tm-k)*t/(1.+tm*t))>0. 	 					; methode interbancaire
		endelse
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	tcalcul=ni*tstep
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
	for i=0,ni do begin
			if method_flag eq 0 then val_exe=(1./bin(0,ni,i)-gmax)*bin(0,ni,i) > 0. $       ;methode actuarielle
		else begin
			val_exe=((1./bin(0,ni,i)-1.-k*tcalcul)/(1./bin(0,ni,i)))>0. 	 	 	; methode interbancaire
		endelse
		tree_call(ni,i)=val_exe>(bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i)))
	endfor
	endelse
endfor
return,tree_call(0,0)
end

function taux_Aput,k,n,na,tstep,duration,tauxzero,p,delta,method_flag
;option sur FRA garantissant un taux k et payant en precompte
;ammericaine a partir de tta jusqu a l echeance
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le taux-prix d exercice; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis jusqu a l echeance
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=n*tstep
tta=na*tstep
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
if method_flag eq 0 then gmax=(1+k)^duration else gmax=1+k*duration
tree_call=dblarr(n+1,n+1)
; fixation de la valeur a l echeance de l option
for i=0,n do begin
	if method_flag eq 0 then tree_call(n,i)=(gmax-1./bin(0,n,i))*bin(0,n,i) > 0. $         ;methode actuarielle
		else begin
			tm=(1./bin(0,n,i)-1.)/t       			 ; methode interbancaire
			tree_call(n,i)=((k-tm)*t/(1.+tm*t))>0. 	 	 ; methode interbancaire
		endelse
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	tcalcul=ni*tstep
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
	for i=0,ni do begin
			if method_flag eq 0 then val_exe=(1./bin(0,ni,i)-gmax)*bin(0,ni,i) > 0. $         ;methode actuarielle
		else begin
			val_exe=((k*tcalcul-1./bin(0,ni,i)+1.)/(1./bin(0,ni,i)))>0. 	 	 ; methode interbancaire
		endelse
		tree_call(ni,i)=val_exe>(bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i)))
	endfor
	endelse
endfor
return,tree_call(0,0)
end




function taux_cap,k,n1,tstep,ncrenaux,tauxzero,p,delta,method_flag
;cap garantissant un taux k maximum et payant en precompte sur des periodes
; dont la largeur est t/ncrenaux
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le taux-prix d exercice; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis jusqu a l echeance
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=n1*tstep
duration=t/ncrenaux
nstep=ceilling(n1/double(ncrenaux))
n=nstep*ncrenaux
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
tree_call=dblarr(n+1,n+1)
options=dblarr(ncrenaux)
if method_flag eq 0 then gmax=(1+k)^duration else gmax=1+k*duration
for nj=0,ncrenaux-1 do begin
	; fixation de la valeur a l echeance de l option
	tree_call(nj*nstep+nstep,0:n)=0.
	for i=0,nj*nstep do begin
		;print,'nj=',nj,' bin=',bin(0,nj*nstep,i),1/gmax,(1./bin(0,nj*nstep,i)-gmax)
		if method_flag eq 0 then tree_call(nj*nstep,i)=(1./bin(0,nj*nstep,i)-gmax)*$
				bin(0,nj*nstep,i) > 0. $         					;methode actuarielle
			else begin
				tm=(1./bin(0,nj*nstep,i)-1.)/duration       			; methode interbancaire
				tree_call(nj*nstep,i)=((tm-k)*duration/(1.+tm*duration))>0. 	 	; methode interbancaire
		endelse
	endfor
	;deroulement de l'arbre
	for ni=nj*nstep-1,0,-1 do begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endfor
options(nj)=tree_call(0,0)
endfor
return,total(options)
end

function taux_floor,k,n1,tstep,ncrenaux,tauxzero,p,delta,method_flag
;cap garantissant un taux k maximum et payant en precompte sur des periodes
; dont la largeur est t/ncrenaux
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le taux-prix d exercice; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis jusqu a l echeance
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=n1*tstep
duration=t/ncrenaux
nstep=ceilling(n1/double(ncrenaux))
n=nstep*ncrenaux
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
tree_call=dblarr(n+1,n+1)
options=dblarr(ncrenaux)
if method_flag eq 0 then gmax=(1+k)^duration else gmax=1+k*duration
for nj=0,ncrenaux-1 do begin
	; fixation de la valeur a l echeance de l option
	tree_call(nj*nstep+nstep,0:n)=0.
	for i=0,nj*nstep do begin
		;print,'nj=',nj,' bin=',bin(0,nj*nstep,i),1/gmax,(gmax-1./bin(0,nj*nstep,i))
		if method_flag eq 0 then tree_call(nj*nstep,i)=(gmax-1./bin(0,nj*nstep,i))*$
				bin(0,nj*nstep,i) > 0. $         					;methode actuarielle
			else begin
				tm=(1./bin(0,nj*nstep,i)-1.)/duration       			; methode interbancaire
				tree_call(nj*nstep,i)=((k-tm)*duration/(1.+tm*duration))>0. 	 	; methode interbancaire
		endelse
	endfor
	;deroulement de l'arbre
	for ni=nj*nstep-1,0,-1 do begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endfor
	options(nj)=tree_call(0,0)
endfor
return,total(options)
end



function oblig_call,k,n,tstep,duree,vect_divid,nbdiv,tauxzero,p,delta,method_flag
;option sur obligation nouvellement cree 
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
t=n*tstep
if nbdiv gt 0 then begin
	bin =binom_tauxn(tauxzero,p,delta,vect_divid,t,tstep)		;coupons flottants
	;visubin,bin,4
endif
bind =binom_taux2(tauxzero,p,delta,duree,t,tstep)
		        ;capital
tree_call=dblarr(n+1,n+1)
; fixation de la valeur a l echeance de l option
for i=0,n do begin
			oblig=0.
			if nbdiv gt 0 then begin
				n_div=n_elements(vect_divid(*,0))
				for idiv=0,n_div-1 do begin
					oblig=oblig+bin(idiv,n,i)*vect_divid(idiv,1)	;les coupons
				endfor
			endif
			oblig=oblig+bind(0,n,i)				; le capital
			tree_call(n,i)=oblig-k > 0.         
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	for i=0,ni do begin
		tree_call(ni,i)=bind(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
	endfor
endfor
return,tree_call(0,0)
end


function oblig_put,k,n,tstep,duree,vect_divid,nbdiv,tauxzero,p,delta,method_flag
;option sur obligation nouvellement cree 
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
t=n*tstep
if nbdiv gt 0 then begin
	bin =binom_tauxn(tauxzero,p,delta,vect_divid,t,tstep)		;coupons flottants
	;visubin,bin,4
endif
bind =binom_taux2(tauxzero,p,delta,duree,t,tstep)
		        ;capital
tree_call=dblarr(n+1,n+1)
; fixation de la valeur a l echeance de l option
for i=0,n do begin
			oblig=0.
			if nbdiv gt 0 then begin
				n_div=n_elements(vect_divid(*,0))
				for idiv=0,n_div-1 do begin
					oblig=oblig+bin(idiv,n,i)*vect_divid(idiv,1)	;les coupons
				endfor
			endif
			oblig=oblig+bind(0,n,i)				; le capital
			tree_call(n,i)=k-oblig> 0.         
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	for i=0,ni do begin
		tree_call(ni,i)=bind(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
	endfor
endfor
return,tree_call(0,0)
end



function oblig_Acall,k,n,na,tstep,duree,vect_divid,nbdiv,tauxzero,p,delta,method_flag
;option sur obligation nouvellement cree 
;americaine a partir de tta (distance depuis 0)
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
t=n*tstep
tta=na*tstep
if nbdiv gt 0 then bin =binom_tauxn(tauxzero,p,delta,vect_divid,t,tstep)		;coupons flottants
bind =binom_taux2(tauxzero,p,delta,duree,t,tstep)		        ;capital
tree_call=dblarr(n+1,n+1)
if nbdiv gt 0 then n_div=n_elements(vect_divid(*,0)) else n_div=0
; fixation de la valeur a l echeance de l option
for i=0,n do begin
			oblig=0.
			if n_div gt 0 then begin
				for idiv=0,n_div-1 do begin
					oblig=oblig+bin(idiv,n,i)*vect_divid(idiv,1)	;les coupons
				endfor
			endif
				oblig=oblig+bind(0,n,i)				; le capital
			tree_call(n,i)=oblig-k > 0.         
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	tcalcul=ni*t/float(n)
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bind(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
		for i=0,ni do begin
			oblig=0.
			if nbdiv gt 0 then begin
				for idiv=0,n_div-1 do begin
					oblig=oblig+bin(idiv,ni,i)*vect_divid(idiv,1)	;les coupons
				endfor
			endif
			oblig=oblig+bind(0,ni,i)				; le capital
			tree_call(ni,i)=oblig-k > bind(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endelse
endfor
return,tree_call(0,0)
end

function oblig_Aput,k,n,na,tstep,duree,vect_divid,nbdiv,tauxzero,p,delta,method_flag
;option sur obligation nouvellement cree 
;americaine a partir de tta (distance depuis 0)
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
t=n*tstep
tta=na*tstep
if nbdiv gt 0 then bin =binom_tauxn(tauxzero,p,delta,vect_divid,t,tstep)		;coupons flottants
bind =binom_taux2(tauxzero,p,delta,duree,t,tstep)		        ;capital
tree_call=dblarr(n+1,n+1)
if nbdiv gt 0 then n_div=n_elements(vect_divid(*,0)) else n_div=0
; fixation de la valeur a l echeance de l option
for i=0,n do begin
			oblig=0.
			if n_div gt 0 then begin
				for idiv=0,n_div-1 do begin
					oblig=oblig+bin(idiv,n,i)*vect_divid(idiv,1)	;les coupons
				endfor
			endif
				oblig=oblig+bind(0,n,i)				; le capital
			tree_call(n,i)=k-oblig > 0.         
endfor
;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	tcalcul=ni*t/float(n)
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bind(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
		for i=0,ni do begin
			oblig=0.
			if nbdiv gt 0 then begin
				for idiv=0,n_div-1 do begin
					oblig=oblig+bin(idiv,ni,i)*vect_divid(idiv,1)	;les coupons
				endfor
			endif
			oblig=oblig+bind(0,ni,i)				; le capital
			tree_call(ni,i)=k-oblig > bind(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endelse
endfor
return,tree_call(0,0)
end
   

function oblig_existante_call,k,n,tstep,duree,vect_divid,nbdiv,tauxzero,p,delta,method_flag
;option sur obligation nouvellement cree 
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
tee=n*tstep
if nbdiv gt 0 then begin
	tdiv=vect_divid(*,0)
	ddiv=vect_divid(*,1)
	nn=n_elements(tdiv)
	indt=where(tdiv gt tee,ndiv)
endif else ndiv = 0
if ndiv gt 0 then begin
	new_tdiv=tdiv(indt)-tee
	new_ddiv=ddiv(indt)
	new_nn=n_elements(new_tdiv)
	eche_vect_divid=fltarr(new_nn,2)
	eche_vect_divid(*,0)=new_tdiv(0:new_nn-1)
	eche_vect_divid(*,1)=new_ddiv(0:new_nn-1)
endif else new_nn=0
return,oblig_call(k,n,tstep,duree-tee,eche_vect_divid,new_nn,tauxzero,p,delta,method_flag)
end


function oblig_existante_put,k,n,tstep,duree,vect_divid,nbdiv,tauxzero,p,delta,method_flag
;option sur obligation nouvellement cree 
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
tee=n*tstep
if nbdiv gt 0 then begin
	tdiv=vect_divid(*,0)
	ddiv=vect_divid(*,1)
	nn=n_elements(tdiv)
	indt=where(tdiv gt tee,ndiv)
endif else ndiv = 0
if ndiv gt 0 then begin
	new_tdiv=tdiv(indt)-tee
	new_ddiv=ddiv(indt)
	new_nn=n_elements(new_tdiv)
	eche_vect_divid=fltarr(new_nn,2)
	eche_vect_divid(*,0)=new_tdiv(0:new_nn-1)
	eche_vect_divid(*,1)=new_ddiv(0:new_nn-1)
endif else new_nn=0
return,oblig_put(k,n,tstep,duree-tee,eche_vect_divid,new_nn,tauxzero,p,delta,method_flag)
end



function oblig_existante_Acall,k,n_echeance,tstep,na,duree,vect_divid,nbdiv,nmax,tauxzero,p,delta,method_flag
;option sur obligation deja existante
;americain a partir de tta (distance depuis 0)
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
;n_echeance=n
tte=n_echeance*tstep
tta=na*tstep
tmax=0
if nbdiv gt 0 then begin
	tmax=(max(vect_divid(*,0))>tte)>duree
	bin =binom_tauxnk(tauxzero,p,delta,tmax,nmax,tte,tstep)
	;visubin,bin,4
endif
tree_call=dblarr(n_echeance+1,n_echeance+1)

n1=ceilling(double(tte)/tstep)+1
n2=ceilling(double(duree)/tstep)+1
n_arbre=ceilling(n2/double(nmax))+1
step=tstep*nmax
; fixation de la valeur a l echeance de l option
indt=where(vect_divid(*,0) gt tte,ndiv)
if ndiv gt 0 then new_nn=n_elements(indt) else new_nn=0
	for i=0,n_echeance do begin
			oblig=0.
			if new_nn gt 0 then begin
				for id=0,new_nn-1 do begin
					idiv=indt(id)
					n_arbre1=int((vect_divid(idiv,0)-tte)/step)
					g1=bin(n_arbre1,n_echeance,i)
					if n_arbre1 lt n_arbre then begin
						g2=bin(n_arbre1+1,n_echeance,i)
						g=g1+((vect_divid(idiv,0)-tte)-n_arbre1*step)/step*(g2-g1)
					endif else g=g1
					oblig=oblig+g*vect_divid(idiv,1)	;les coupons
 				endfor
			endif
			n_arbre3=int((duree-tte)/step)
			g3=bin(n_arbre3,n_echeance,i)
			if n_arbre3 lt n_arbre then begin
						g4=bin(n_arbre3+1,n_echeance,i)
						g=g3+((duree-tte)-n_arbre3*step)/step*(g4-g3)
				endif else g=g3
			oblig=oblig+g					; le capital
;print,'  capital bin=',g,'oblig=',oblig
 			tree_call(n_echeance,i)=oblig-k > 0.         
endfor
;print,'val a l eche.:',transpose(tree_call(n_echeance,*))
;deroulement de l'arbre
for ni=n_echeance-1,0,-1 do begin
	tcalcul=ni*tstep
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(0,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
		if ndiv gt 0 then indt=where(vect_divid(*,0) gt tcalcul,new_nn) else new_nn=0
		for i=0,ni do begin
			oblig=0.
			if new_nn gt 0 then begin
				for id=0,new_nn-1 do begin
					idiv=indt(id)
					n_arbre1=int((vect_divid(idiv,0)-tcalcul)/step)
					g1=bin(n_arbre1,ni,i)
					if n_arbre1 lt n_arbre then begin
						g2=bin(n_arbre1+1,ni,i)
						g=g1+((vect_divid(idiv,0)-tcalcul)-n_arbre1*step)/step*(g2-g1)
					endif else g=g1
					oblig=oblig+g*vect_divid(idiv,1)	;les coupons
				endfor
			endif
			n_arbre3=int((duree-tcalcul)/step)
			g3=bin(n_arbre3,ni,i)
			if n_arbre3 lt n_arbre then begin
						g4=bin(n_arbre3+1,ni,i)
						g=g3+((duree-tcalcul)-n_arbre3*step)/step*(g4-g3)
				endif else g=g3

			oblig=oblig+g					; le capital
			tree_call(ni,i)=(oblig-k)>$
				(bin(0,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i)))
		endfor
	endelse
endfor
return,tree_call(0,0)
end

function oblig_existante_Aput,k,n_echeance,tstep,na,duree,vect_divid,nbdiv,nmax,tauxzero,p,delta,method_flag
;option sur obligation deja existante
;americain a partir de tta (distance depuis 0)
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
tte=n_echeance*tstep
tta=na*tstep
tmax=0
if nbdiv gt 0 then begin
	tmax=(max(vect_divid(*,0))>tte)>duree
	bin =binom_tauxnk(tauxzero,p,delta,tmax,nmax,tte,tstep)
	;visubin,bin,4
endif
tree_call=dblarr(n_echeance+1,n_echeance+1)
n1=ceilling(double(tte)/tstep)+1
n2=ceilling(double(duree)/tstep)+1
n_arbre=ceilling(n2/double(nmax))+1
step=tstep*nmax
; fixation de la valeur a l echeance de l option
indt=where(vect_divid(*,0) gt tte,ndiv)
if ndiv gt 0 then new_nn=n_elements(indt) else new_nn=0
	for i=0,n_echeance do begin
			oblig=0.
			if new_nn gt 0 then begin
				for id=0,new_nn-1 do begin
					idiv=indt(id)
					n_arbre1=int((vect_divid(idiv,0)-tte)/step)
					g1=bin(n_arbre1,n_echeance,i)
					if n_arbre1 lt n_arbre then begin
						g2=bin(n_arbre1+1,n_echeance,i)
						g=g1+((vect_divid(idiv,0)-tte)-n_arbre1*step)/step*(g2-g1)
					endif else g=g1
					oblig=oblig+g*vect_divid(idiv,1)	;les coupons
 				endfor
			endif
			n_arbre3=int((duree-tte)/step)
			g3=bin(n_arbre3,n_echeance,i)
			if n_arbre3 lt n_arbre then begin
						g4=bin(n_arbre3+1,n_echeance,i)
						g=g3+((duree-tte)-n_arbre3*step)/step*(g4-g3)
				endif else g=g3
			oblig=oblig+g					; le capital
;print,'  capital bin=',g,'oblig=',oblig
 			tree_call(n_echeance,i)=k-oblig > 0.         
endfor
;print,'val a l eche.:',transpose(tree_call(n_echeance,*))
;deroulement de l'arbre
for ni=n_echeance-1,0,-1 do begin
	tcalcul=ni*tstep
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(0,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
		if ndiv gt 0 then indt=where(vect_divid(*,0) gt tcalcul,new_nn) else new_nn=0
		for i=0,ni do begin
			oblig=0.
			if new_nn gt 0 then begin
				for id=0,new_nn-1 do begin
					idiv=indt(id)
					n_arbre1=int((vect_divid(idiv,0)-tcalcul)/step)
					g1=bin(n_arbre1,ni,i)
					if n_arbre1 lt n_arbre then begin
						g2=bin(n_arbre1+1,ni,i)
						g=g1+((vect_divid(idiv,0)-tcalcul)-n_arbre1*step)/step*(g2-g1)
					endif else g=g1
					oblig=oblig+g*vect_divid(idiv,1)	;les coupons
				endfor
			endif
			n_arbre3=int((duree-tcalcul)/step)
			g3=bin(n_arbre3,ni,i)
			if n_arbre3 lt n_arbre then begin
						g4=bin(n_arbre3+1,ni,i)
						g=g3+((duree-tcalcul)-n_arbre3*step)/step*(g4-g3)
				endif else g=g3

			oblig=oblig+g					; le capital
			tree_call(ni,i)=(k-oblig)>$
				(bin(0,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i)))
		endfor
	endelse
endfor
return,tree_call(0,0)
end



function oblig_existante_AFcall,k,n_echeance,tstep,na,nf,duree,vect_divid,nbdiv,nmax,tauxzero,p,delta,method_flag
;option sur futur d obligation deja existante
;americain a partir de na*tstep (distance depuis 0)
;le contrat de futur est a echeance en  nf*tstep (distance depuis 0)
; duree :duree de l obligation nouvellement cree
; vect_divid ; le vecteur [[date_relative,coupon],...] des coupons  de l obligation
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; t est l echeance de l option
; k est le prix d exercice de l option : il est fixe 
; n est le nombre d'etapes du treillis
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
;vect_divid=format_coupons(duree,coupon,delai_prem_coupon,method_flag)
;n_echeance=n
tte=n_echeance*tstep
tta=na*tstep
ttf=nf*tstep
tmax=0
if nbdiv gt 0 then begin
	tmax=(max(vect_divid(*,0))>ttf)>duree
	bin =binom_tauxnk(tauxzero,p,delta,tmax,nmax,ttf,tstep)
	;visubin,bin,4
endif
tree_call=dblarr(n_echeance+1,n_echeance+1)
tree_future=dblarr(nf+1,nf+1)
n1=ceilling(double(ttf)/tstep)+1
n2=ceilling(double(duree)/tstep)+1
n_arbre=ceilling(n2/double(nmax))+1
step=tstep*nmax
; fixation de la valeur a l echeance du future
indt=where(vect_divid(*,0) gt ttf,ndiv)
print,'div=',vect_divid(*,0)
print,'indt=',indt
if ndiv gt 0 then new_nn=n_elements(indt) else new_nn=0
	for i=0,nf do begin
			oblig=0.
			if new_nn gt 0 then begin
				for id=0,new_nn-1 do begin
					idiv=indt(id)
					n_arbre1=int((vect_divid(idiv,0)-ttf)/step)
					g1=bin(n_arbre1,nf,i)
					if n_arbre1 lt n_arbre then begin
						g2=bin(n_arbre1+1,nf,i)
						g=g1+((vect_divid(idiv,0)-ttf)-n_arbre1*step)/step*(g2-g1)
					endif else g=g1
					oblig=oblig+g*vect_divid(idiv,1)	;les coupons
 				endfor
			endif 
			n_arbre3=int((duree-ttf)/step)
			g3=bin(n_arbre3,nf,i)
			if n_arbre3 lt n_arbre then begin
						g4=bin(n_arbre3+1,nf,i)
						g=g3+((duree-ttf)-n_arbre3*step)/step*(g4-g3)
				endif else g=g3
			oblig=oblig+g					; le capital
 			tree_future(nf,i)=oblig        
endfor
print,'val a l eche. du future de l oblig:',transpose(tree_future(nf,0:nf))
;deroulement de l'arbre du futur jusqu a l echeance de l option au minimum
for ni=nf-1,0,-1 do begin
	for i=0,ni do begin
		tree_future(ni,i)=(p*tree_future(ni+1,i+1)+(1.-p)*tree_future(ni+1,i))
	endfor
endfor
print,'valeur du future :',tree_future(0,0)

;fixation de la valeur de l option a l echeance de celle-ci
for i=0,n_echeance do tree_call(n_echeance,i)=tree_future(n_echeance,i)-k>0
print,'val a l eche. de l option  du future:',transpose(tree_call(n_echeance,0:n_echeance))
;deroulement de l arbre pour l option
for ni=n_echeance-1,0,-1 do begin
	tcalcul=ni*tstep
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(0,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
		if ndiv gt 0 then indt=where(vect_divid(*,0) gt tcalcul,new_nn) else new_nn=0
		for i=0,ni do begin

			tree_call(ni,i)=(tree_future(ni,i)-k)>$
				(bin(0,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i)))
		endfor
	endelse
endfor
return,tree_call(0,0)
end

function zerocoupon_AFcall,k,n,tstep,na,nf,duration,tauxzero,p,delta
;option sur pibor future garantissant un taux k et payant en precompte
;ammericaine a partir de tta=na*tstep jusqu a l echeance tte=n*tstep
; option sur future dont l echeance est en ttf=nf*tstep
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; k est le  prix d exercice exprime par un niveau de prix de zero coupon futur
; le paiement correspond a un financement sur une periode de 'duration'
; 	d un differentiel de taux 
; n est le nombre d'etapes du treillis jusqu a l echeance
; p et delta sont les parametres du modele
; method_flag vaut 1 si la methode d 'actualisation est la methode interbancaire avec un taux 
; 	d interet lineaire dans le temps ou 0 si on utilise la methode d  actualisation exponentielle
; le record reuter des option pibor3mois est  Q0#PIBH3+
t=n*tstep
tta=na*tstep
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
tree_future=dblarr(nf+1,nf+1)
tree_call=dblarr(n+1,n+1)

;fixation de la valeur a l echeance du future
for i=0,nf do tree_future(nf,i)=bin(0,nf,i)

;deroulement de l arbre du future
for ni = nf-1,0,-1 do for i= 0,ni do tree_future(ni,i)=p*tree_future(ni+1,i+1)+(1-p)*tree_future(ni+1,i)
print,'valeur du future :actuariel=',tree_future(0,0)^(1/duration),'monetaire=',1.-4*(1-tree_future(0,0))
	
; fixation de la valeur a l echeance de l option
for i=0,n do tree_call(n,i)=(tree_future(n,i)-k)*bin(0,n,i)>0.

;deroulement de l'arbre
for ni=n-1,0,-1 do begin
	tcalcul=ni*tstep
	if tcalcul lt tta then begin
		for i=0,ni do begin
			tree_call(ni,i)=bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i))
		endfor
	endif else begin
	for i=0,ni do begin
		val_exe=(tree_future(ni,i)-k)*bin(0,ni,i) >0
 		tree_call(ni,i)=val_exe>(bin(1,ni,i)*(p*tree_call(ni+1,i+1)+(1-p)*tree_call(ni+1,i)))
	endfor
	endelse
endfor
return,tree_call(0,0)
end

function zerocoupon_future,tstep,nf,duration,tauxzero,p,delta
;future sur pibor  
; option sur future dont l echeance est en ttf=nf*tstep
;tauxzero est une liste dont chaque element est de la sorte:
;    [temps  , yield ]
; duration est la duration de l instrument de taux qui est l'actif sousjacent d l option
; n est le nombre d'etapes du treillis jusqu a l echeance
t=nf*tstep
bin =binom_taux2(tauxzero,p,delta,duration,t,tstep)
tree_future=dblarr(nf+1,nf+1)
;fixation de la valeur a l echeance du future
for i=0,nf do tree_future(nf,i)=bin(0,nf,i)
print,transpose(tree_future(nf,*))
;deroulement de l arbre du future
for ni = nf-1,0,-1 do for i= 0,ni do tree_future(ni,i)=p*tree_future(ni+1,i+1)+(1-p)*tree_future(ni+1,i)
return,tree_future(0,0)^(1/duration) 	
end



function  taux_recalage_aux,opt_list,tstep,tauxzero,p,delta
;fonction qui calcule la valeur de la distance entre les previsions theoriques et les prix de marche
valmax=0.
slist=size(opt_list)
if slist(0) eq 1 then nb=1 else nb=slist(1)
for i = 0,nb-1 do begin
	param=opt_list(i,*)
	nech=int(param(1)/tstep)
	k=param(2)
	case 1 of
		param(0) eq 11 : begin
			f1=taux_Acall(k,nech,0,tstep,91./360.,tauxzero,p,delta,0)
			f2=taux_Acall(k,nech+1,0,tstep,91./360.,tauxzero,p,delta,0)
			end
		param(0) eq 12 : begin
			f1=taux_Aput(k,nech,0,tstep,91./360.,tauxzero,p,delta,0)
			f2=taux_Aput(k,nech+1,0,tstep,91./360.,tauxzero,p,delta,0)
			end
		else :print,'type non connu:',param(0)
	endcase
	f=f1+(f2-f1)*(param(1)/tstep-nech)
print,'f=',f,'prime recherche=',param(3)
	valmax=max([valmax,abs(f-param(3))])
endfor
print,valmax
return,valmax
end

function  taux_recalage_aux1,opt_list,tstep,tauxzero,p,delta
;fonction qui calcule la valeur de la distance entre les previsions theoriques et les prix de marche
valmax=0.
slist=size(opt_list)
if slist(0) eq 1 then nb=1 else nb=slist(1)
for i = 0,nb-1 do begin
	param=opt_list(i,*)
	nech=int(param(1)/tstep)
	k=param(2)
	case 1 of
		param(0) eq 1 : begin
			f1=taux_Acall(k,nech,0,tstep,91./360.,tauxzero,p,delta,0)
			f2=taux_Acall(k,nech+1,0,tstep,91./360.,tauxzero,p,delta,0)
			end
		param(0) eq 2 : begin
			f1=taux_Aput(k,nech,0,tstep,91./360.,tauxzero,p,delta,0)
			f2=taux_Aput(k,nech+1,0,tstep,91./360.,tauxzero,p,delta,0)
			end
		else :print,'type non connu:',param(0)
	endcase
	f=f1+(f2-f1)*(param(1)/tstep-nech)
print,'f=',f,'prime recherche=',param(3)
	valmax= f-param(3)
endfor
print,valmax
return,valmax
end


pro print_delta_prix,opt_list,tstep,tauxzero,p,delta
;fonction qui affiche les previsions theoriques et les prix de marche
valmax=0.
slist=size(opt_list)
if slist(0) eq 1 then nb=1 else nb=slist(1)
print,'                prix de marche / prix theorique      prime de marche / prime theorique'
for i = 0,nb-1 do begin
	param=opt_list(i,*)
	nech=int(param(1)/tstep)
	k=param(2)
	case 1 of
		param(0) eq 1 : begin
			f1=taux_Acall(k,nech,0,tstep,91./360.,tauxzero,p,delta,0)
			f2=taux_Acall(k,nech+1,0,tstep,91./360.,tauxzero,p,delta,0)
			end
		param(0) eq 2 : begin
			f1=taux_Aput(k,nech,0,tstep,91./360.,tauxzero,p,delta,0)
			f2=taux_Aput(k,nech+1,0,tstep,91./360.,tauxzero,p,delta,0)
			end
		else :print,'type non connu:',param(0)
	endcase
	prime=f1+(f2-f1)*(param(1)/tstep-nech)
	prix_theo=((1+prime)^(1./param(1)) - 1.)*100
	prix_march=((1+param(3))^(1./param(1)) - 1.)*100
	print,'prix ',i,'=',prix_march,' / ',prix_theo,'  ',param(3),prime
endfor
end

function taux_recalage,opt_list,tstep,tauxzero,p0,delta0
;fonction qui trouve le p etr le delta implicite dans les prix de marche fourni
;par les options sur pibor de opt_list dont le format est :
; opt_list=transpose([[type,echeance,strike,prix],...])
; echeance en annee , strike : coef d actualisation en % exemple  92 pour 0.92 ,
; type = 1 pour call, 2 pour put,  prix = 0.60 pour 0.6 % de taux sur la periode jusqu a echeance

;          premierement
;          ------------
;on supose p bon et on minimise la fonction auxiliaire sur delta
p=p0
ax=double(0.99)
bx=double(delta0)
cx=double(0.9999999)
tol=0.0001
r=double(0.61803399)
c=1.-r
x0=ax
x3=cx
if abs(cx-bx) gt abs(bx-ax) then begin
	x1=bx
	x2=bx+c*(cx-bx)
endif else begin
	x2=bx
	x1=bx-c*(bx-ax)
endelse
f1=taux_recalage_aux(opt_list,tstep,tauxzero,p,x1)
f2=taux_recalage_aux(opt_list,tstep,tauxzero,p,x2)
dlist=[0.998,0.999,0.9999,0.99999]
for ii=0,n_elements(dlist)-1 do begin
	fd=taux_recalage_aux1(opt_list,tstep,tauxzero,p,dlist(ii))
	print,'essai',dlist(ii),fd
endfor
	
while abs(x3-x0) gt tol*(abs(x1)+abs(x2)) do begin
print,x1,x2,f1,f2
	if f2 lt f1 then begin
		x0=x1
		x1=x2
		x2=r*x1+c*x3
		f0=f1
		f1=f2
		f2=taux_recalage_aux(opt_list,tstep,tauxzero,p,x2)
	endif else begin
		x3=x2
		x2=x1
		x1=r*x2+c*x0
		f3=f2
		f2=f1
		f1=taux_recalage_aux(opt_list,tstep,tauxzero,p,x1)
	endelse
endwhile
if f1 lt f2 then begin
	golden=f1
	xmin=x1
endif else begin
	golden=f2
	xmin=x2
endelse
print,'minimum=',golden,'   valeur optimale=',xmin
print_delta_prix,opt_list,tstep,tauxzero,p,xmin
delta=xmin
return,[p,delta]
end

;tentative d extraction des paramtre p et delta implicite du marche le 6/1/93
pro test_recalage
; on cherche sur la page RGREL les information milieu de fourchette des bons
taux=transformat_taux([	[1/360.,12+1/8.],$
			[1/12.,13+1/4.],$
			[2/12.,12+3/4.],$
			[3/12.,12+3/4.],$
			[6/12.,11+1/4.],$
			[9/12.,10+1/4.],$
			[1.,10.],$
			[2.,9.25],$
			[3.,9.05],$
			[4.,8.90],$
			[5.,8.90],$
			[7.,9.],$
			[8.,9.05],$
			[9.,9.10],$
			[10.,9.15]])
;on cherche des prix de marche (sur les records  :
; le record reuter des option pibor3mois mars93 est  Q0#PIBH3+
; le record reuter des option pibor3mois juin93 est  Q0#PIBM3+
; le record reuter des option pibor3mois sept93 est  Q0#PIBV3+
; le record reuter des option pibor3mois dec93  est  Q0#PIBZ3+  par exemple)
; prixmarche=transpose([[type,echeance,strike,prix],...])
; echeance en annee , strike : coef d actualisation en % exemple  92 pour 0.92 ,
; type = 11 pour call pibor, 12 pour put pibor 21 pour call matif 22 pour put matif
; ,  prix = 0.60 pour 0.6 % de taux sur la periode jusqu a echeance
t0=jj(100193)
prixmarche=transformat_prix([$
			 [11,ech('mar3',11),90.5,0.555]$
		;	,[11,ech('mar3',11),91.,0.330]$		
		;	,[11,ech('mar3',11),91.4,0.180]$
		;	,[11,ech('mar3',11),92.,0.075]$
		;	,[12,ech('mar3',12),89.,0.135]$
		;	,[12,ech('mar3',12),90.,0.270]$
		;	,[12,ech('mar3',12),91.,0.645]$
		;	,[11,ech('jun3',11),91.8,0.625]$
		;	,[11,ech('jun3',11),92.,0.5]$
		;	,[11,ech('jun3',11),92.3,0.365]$
		;	,[11,ech('jun3',11),92.5,0.285]$
		;	,[11,ech('sep3',11),92.5,0.525]$
		;	,[11,ech('sep3',11),92.6,0.55]$
		;	,[22,ech('mar3',22),109.,0.09]$
		;	,[22,ech('mar3',22),110.,0.16]$
		;	,[21,ech('mar3',21),112.,1.41]$
		;	,[22,ech('mar3',22),112.,0.48]$
		;	,[21,ech('mar3',21),113.,0.80]$
		;	,[21,ech('mar3',21),114.,0.35]$
		;	,[21,ech('mar3',21),115.,0.14]$
		;	,[22,ech('jun3',22),112.,0.82]$
		;	,[22,ech('jun3',22),113.,1.15]$
		;	,[21,ech('jun3',21),114.,1.36]$
		;	,[21,ech('jun3',21),115.,0.95]$
		;	,[21,ech('jun3',21),116.,0.62]$
		;	,[22,ech('dec3',22),111.,0.055]$
		;	,[21,ech('dec3',21),113.,0.37]$
		;	,[21,ech('dec3',21),114.,0.07]$

			],t0)

tauxzero=couponszero(taux)		;on calcule les coupons zero
p0=0.2
delta0=0.9995
tstep=0.33/12.
resu= taux_recalage(prixmarche,tstep,tauxzero,p0,delta0)
end


;tentative d extraction des paramtre p et delta implicite du marche le 6/1/93
pro test_recalage1
; on cherche sur la page RGREL les information milieu de fourchette des bons
t0=jj(100193)
taux=transformat_taux([	[1/360.,12+1/8.],$
			[1/12.,13+1/4.],$
			[2/12.,12+3/4.],$
			[3/12.,12+3/4.],$
			[6/12.,11+1/4.],$
			[9/12.,10+1/4.],$
			[1.,10.],$
			[2.,9.25],$
			[3.,9.05],$
			[4.,8.90],$
			[5.,8.90],$
			[7.,9.],$
			[8.,9.05],$
			[9.,9.10],$
			[10.,9.15]])
 			param=[1,ech('mar3',1),90.5,0.555]
	type=param(0)			; le type d option
	tech=(param(1)-t0)/365.			; temps restant a courir jusque a l echeance
	k=1./(param(2)/100.)-1.	; taux annuel sous entendu
	prime=(1+double(param(3)/100.))^tech-1. ; calcul de la prime
	 tauxzero=couponszero(taux)		;on calcule les coupons zero
p=0.2
delta0=0.9999999
tstep=0.1/12.
nech=int((tech)/tstep)
na=10
print,'nech=',nech
print,'delta=',delta0
print,'k=',k
	case 1 of
		type eq 1 : begin
			f1=taux_Acall(k,nech,na,tstep,91./360.,tauxzero,p,delta0,0)
			f2=taux_Acall(k,nech+1,na,tstep,91./360.,tauxzero,p,delta0,0)
			end
		type eq 2 : begin
			f1=taux_Aput(k,nech,na,tstep,91./360.,tauxzero,p,delta0,0)
			f2=taux_Aput(k,nech+1,na,tstep,91./360.,tauxzero,p,delta0,0)
			end
		else :print,'type non connu:',type
	endcase
	f=f1+(f2-f1)*(param(1)/tstep-nech)
print,' prime theorique=',f,'  prime de marche=',prime
end





pro testtaux_echeance
taux=transpose([$
[1./12.,0.992],$
[2./12.,0.984],$
[3./12.,0.976],$
[4./12.,0.969],$
[5./12.,0.961],$
[6./12.,0.953]$
;,[7./12,0.945],[8./12.,0.937],[9./12,0.93],[10./12.,0.922]$
])
staux=size(taux)
ntaux=staux(1)
taux1=taux
for i=0,ntaux-1 do taux1(i,1)=taux(i,1)^(-1/taux(i,0))-1.
vol=0.07
q=0.5
r=0.1
duration=3./12.
t=4./12.
k=0.09
tstep=0.5/12.
p=0.3
delta=0.98
for nt=4,15 do begin
print,'n=',nt,'  t=',nt*tstep,'euro=',taux_call(k,nt,tstep,duration,taux1,p,delta ,0)*100,' % DU NOMINAL',$
	'amer=',taux_acall(k,nt,0,tstep,duration,taux1,p,delta ,0)*100,' % DU NOMINAL'
endfor
end

pro testtaux_strike
taux=transpose([[1./12.,0.992],[2./12.,0.984],[3./12.,0.976],[4./12.,0.969],[5./12.,0.961],[6./12.,0.953]])
taux1=taux
for i=0,5 do taux1(i,1)=taux(i,1)^(-1/taux(i,0))-1.
vol=0.07
q=0.5
r=0.1
duration=3./12.
nt=40
k=0.09
tstep=1./120.
w =calcule_p_delta(tstep*100,100,vol,q,r)
p=w(0)
delta=w(1)
for strike=0.08,0.11,0.002 do begin
print,'europeenne actuarielle  strike=',strike,'call=',taux_call(strike,nt,tstep,duration,taux1,p,delta ,0),$
' put=',taux_put(strike,nt,tstep,duration,taux1,p,delta ,0)
endfor
end

pro testv,duration
taux=transpose([[1./12.,0.992],[2./12.,0.984],[3./12.,0.976],[4./12.,0.969],[5./12.,0.961],[6./12.,0.953]] )
tauxzero=taux
for i=0,5 do tauxzero(i,1)=taux(i,1)^(-1/taux(i,0))-1.
p=0.5
delta=0.998
duration=1./12.
t1=6./12.
n=6
tstep=t1/n
ss=size(tauxzero)
nn=ss(1)
dt=tauxzero(*,0)
cp1=tauxzero(*,1)
cp=dblarr(nn)
for i=0,nn-1 do cp(i)=1./(1.+double(cp1(i)))^dt(i)
n=ceilling((double(t1+duration)/tstep))+2
n1=ceilling(double(t1)/tstep)+1
tmps=(1.+dindgen(n))*tstep
cps=dblarr(n)
for i=0,n-1 do begin
	if tmps(i) le dt(0) then cps(i)=(1.+cp1(0))^(-tmps(i))
	;if tmps(i) ge dt(nn-1) then cps(i)=cp(nn-1)
	if tmps(i) ge dt(nn-1) then cps(i)=(1.+cp1(nn-1))^(-tmps(i))
	if (tmps(i) gt dt(0)) and (tmps(i) lt dt(nn-1)) then begin
		indless=where(dt lt tmps(i))
		indt1=max(indless)
		indsup=where(dt ge tmps(i))
		indt2=min(indsup)
		;cps(i)=cp(indt1)+((cp(indt2)-cp(indt1))/(dt(indt2)-dt(indt1)))*(tmps(i)-dt(indt1))
		tt1=double(dt(indt1))
		tt2=double(dt(indt2))
			 ; on interpolle les taux
		cps(i)=(1.+(cp1(indt1)*(1.-(tmps(i)-tt1)/(tt2-tt1))+cp1(indt2)*((tmps(i)-tt1)/(tt2-tt1))))^(-tmps(i)) 
	endif
endfor
bin1=binom_taux1(cps,p,delta,1)
display,bin1
end

pro testmatif_call
tauxbon=transformat_taux([	[1/360.,12+1/8.],$
			[1/12.,11+11/16.],$
			[2/12.,11+9/16.],$
			[3/12.,11+9/16.],$
			[6/12.,10+9/16.],$
			[9/12.,9+7/8.],$
			[1.,9.68],$
			[2.,8.90],$
			[3.,8.70],$
			[4.,8.55],$
			[5.,8.55],$
			[7.,8.65],$
			[8.,8.70],$
			[9.,8.65],$
			[10.,8.65]])
tauxzero=couponszero(tauxbon)
tauxzero1=transformat_taux([	[1/360.,12+1/8.]$
			,[3/12.,12.]$
			,[6/12.,11.]$
			,[1.,9.85]$
			,[2.,8.15]$
			,[3.,7.95]$
			,[4.,7.85]$
			,[5.,7.85]$
			,[6.,7.95]$
			,[7.,8.]$
			,[8.,8.15]$
			,[9.,8.2]$
			,[10.,8.3]$
			,[15.,8.5]$
			,[20.,8.7]$
			])
printcourbetaux,tauxzero,2
j0=jj(120193)

; caracterisitique de lq cheapest
tc=(jj(270593)-j0)/365.
tech_support=7.+tc
vect_divid=transpose([[tc,0.1],[1.+tc,0.1],[2.+tc,0.1],[3.+tc,0.1],[4.+tc,0.1],[5.+tc,0.1],[6.+tc,0.1],[7.+tc,0.1]])
nbdiv=8

tstep=.05
tech_future=(jj(3100393)-j0)/365.
n_future=int(tech_future/tstep)
tech_option=(jj(280293)-j0)/365.
n_option=int(tech_option/tstep)
na=0
; calcul du niveau d exercice
niv1=112.
k1= (niv1/100.)
niv2=113.
k2= (niv2/100.)
nmax=3
p=0.3
delta=0.99965
method_flag=0

res=oblig_existante_AFcall(k1,n_option,tstep,na,n_future,tech_support,vect_divid,nbdiv,nmax,tauxzero,p,delta,method_flag)
print,'option matif au niveau ',niv1,'  prime=',res
res=oblig_existante_AFcall(k2,n_option,tstep,na,n_future,tech_support,vect_divid,nbdiv,nmax,tauxzero,p,delta,method_flag)
print,'option matif au niveau ',niv2,'  prime=',res
end 

pro testpibor_call
tauxbon1=transformat_taux([	[1/360.,12+1/2.],$
			[1/12.,12+1/8.],$
			[2/12.,11+3/4.],$
			[3/12.,11+5/8.],$
			[6/12.,10+3/8.],$
			[9/12.,10+5/8.],$
			[1.,9.50],$
			[2.,8.70],$
			[3.,8.50],$
			[4.,8.35],$
			[5.,8.35],$
			[7.,8.45],$
			[8.,8.50],$
			[9.,8.55],$
			[10.,8.55]])
tauxbon2=transformat_taux([	[1/360.,12+1/2.],$
			[1/12.,12+1/8.],$
			[2/12.,12.],$
			[3/12.,11+7/8.],$
			[6/12.,10+3/4.],$
			[9/12.,9+1/8.],$
			[1.,9.87],$
			[2.,9.10],$
			[3.,8.90],$
			[4.,8.75],$
			[5.,8.75],$
			[7.,8.85],$
			[8.,8.90],$
			[9.,8.95],$
			[10.,8.95]])
tauxbon=tauxbon1*0.5+tauxbon2*0.5
tauxzero=couponszero(tauxbon)
;printcourbetaux,tauxzero,2
j0=jj(120193)
tstep=.003
tech_future=(jj(150393)-j0)/365.
n_future=int(tech_future/tstep)
tech_option=(jj(150393)-j0)/365.
print,'echeance de l option=',tech_option
n_option=int(tech_option/tstep)
na=0
duration=3/12.
; calcul du niveau d exercice
niv1=90.5
k1= 1.-(1.-(niv1/100.))*duration
niv2=91.5
k2= 1.-(1.-(niv2/100.))*duration
nmax=3
p=0.5
delta=0.99993
method_flag=0
taux_court=0.1
vol1=sqrt((p*(1-p)*((alog(delta))^2))/(taux_court^2*(tstep^3)))
vol2=sqrt(p*(1-p)*((alog(delta))^2)*n_option^2)
print,'vol taux=',vol1,vol2
vol=0.1
res1=zerocoupon_AFcall(k1,n_option,tstep,na,n_future,duration,tauxzero,p,delta)
print,'option pibor au niveau ',niv1,'  prime=',res1
f1=zerocoupon_future(tstep,n_future,duration,tauxzero,p,delta)
f=0.903
vol=0.0344
print,'pour le modele de black: futur=',f,'  volatilite=',vol,'  alors que ho et lee future=',f1
print,'modele de black=',callF(f,niv1/100.,n_option*tstep,taux_court,vol)
res2=zerocoupon_AFcall(k2,n_option,tstep,na,n_future,duration,tauxzero,p,delta)
print,'option pibor au niveau ',niv2,'  prime=',res2
print,'modele de black=',callF(f,niv2/100.,n_option*tstep,taux_court,vol)
end


pro testpibor_future,delta
tauxbon1=transformat_taux([	[1/360.,12+1/2.],$
			[1/12.,12+1/8.],$
			[2/12.,11+3/4.],$
			[3/12.,11+5/8.],$
			[6/12.,10+3/8.],$
			[9/12.,10+5/8.],$
			[1.,9.50],$
			[2.,8.70],$
			[3.,8.50],$
			[4.,8.35],$
			[5.,8.35],$
			[7.,8.45],$
			[8.,8.50],$
			[9.,8.55],$
			[10.,8.55]])
tauxbon2=transformat_taux([	[1/360.,12+1/2.],$
			[1/12.,12+1/8.],$
			[2/12.,12.],$
			[3/12.,11+7/8.],$
			[6/12.,10+3/4.],$
			[9/12.,9+1/8.],$
			[1.,9.87],$
			[2.,9.10],$
			[3.,8.90],$
			[4.,8.75],$
			[5.,8.75],$
			[7.,8.85],$
			[8.,8.90],$
			[9.,8.95],$
			[10.,8.95]])
tauxbon=tauxbon1*0.5+tauxbon2*0.5
tauxzero=couponszero(tauxbon)
;printcourbetaux,tauxzero,2
j0=jj(120193)
tstep=.01
tech_future=(jj(150393)-j0)/365.
n_future=int(tech_future/tstep)
duration=3/12.
; calcul du niveau d exercice
niv1=90.5
k1= 1.-(1.-(niv1/100.))*duration
niv2=91.5
k2= 1.-(1.-(niv2/100.))*duration
nmax=3
p=0.5
f1=zerocoupon_future(tstep,n_future,duration,tauxzero,p,delta)
print,'  alors que ho et lee future=',f1
end
