 


;                 ***************************** S W A P S **********************


;calcule la courbe te taux forward  de terme tf a partir d une courbe zero coupons

function taux_forward,curvezero,tf
x=curvezero(*,0)
y=curvezero(*,1)
zsize=size(curvezero)
n=zsize(1)
nt=int(curvezero(n-1,0)/float(tf))
xnew=(findgen(nt)+1)*tf
ynew=spline(x,y,xnew)
res=fltarr(nt-1,2)
res(0,0)=tf*0.0001
res(0,1)=ynew(0)
for i=1,nt-2 do begin
	res(i,0)=i*tf
	res(i,1)=((1.+ynew(i+1))^(i+1))/((1.+ynew(i))^(i))-1.
endfor
return,res
end


;valeur d un bond
; spread est un spread de taux caracterisant la signature
;les taux sont tous des taux discrets
function bond,notionnel,maturity,periodpaiements,fixedrated,curvezero,spreadd
fixedrate=alog(1.+fixedrated)
spread=alog(1.+spreadd)
if maturity lt 0 then return,0
nbp=int((maturity-0.000001)/periodpaiements)+1
paymentdate=findgen(nbp)*periodpaiements+(maturity-(nbp-1)*periodpaiements)
sum=0
for i=0,nbp-1 do begin
	rzd=interpolation(curvezero,paymentdate(i))
	rz=alog(1.+rzd)
	cash_fixe=(exp(fixedrate*periodpaiements)-1.)*notionnel
	sum=sum+cash_fixe*exp(-(rz+spread)*paymentdate(i))
endfor
sum=sum+notionnel*exp(-(rz+spread)*maturity)
return,sum
end

;valeur d une rente (bond sans echange de capital)
function rente,notionnel,maturity,periodpaiements,fixedrated,curvezero,spreadd
fixedrate=alog(1.+fixedrated)
spread=alog(1.+spreadd)
if maturity lt 0 then return,0
nbp=int((maturity-0.000001)/periodpaiements)+1
paymentdate=findgen(nbp)*periodpaiements+(maturity-(nbp-1)*periodpaiements)
sum=0
for i=0,nbp-1 do begin
	rzd=interpolation(curvezero,paymentdate(i))
	rz=alog(1.+rzd)
	cash_fixe=(exp(fixedrate*periodpaiements)-1.)*notionnel
	sum=sum+cash_fixe*exp(-(rz+spread)*paymentdate(i))
endfor
return,sum
end

;valeur d un swap

; le side 1 est le taux float donne en interpolant la courbe curveF (forward rate)
; le side 2 est le taux fixe k
; paymntdate est la liste des dates de payment
; curveZ est la coube des zero coupons
; curveN est la liste des notionnels applicables (meme longueur que paymntdate)
; les taux sont suposes s appliquer exponentiellement sur les durees exactes
; k est un vecteur qui doit avoir la meme dim que paymentdate et que curveN (notionnel)
; les taux entres sont tous des taux discrets.
function swap,k,paymntdate,curveF,curveZ,curveN,spread
	sum=0
	for i=0,n_elements(paymntdate)-1 do begin
		if i eq 0 then debut=double(0.000001) else debut=paymntdate(i-1)
		rf=alog(1.+interpolation(curveF,paymntdate(i)))
		rz=alog(1.+interpolation(curveZ,paymntdate(i)))
		cash_float=(exp(rf*(paymntdate(i)-debut))-1.)*curveN(i)
		cash_fixe=(exp(alog((1+k(i))*(1+double(spread)))*(paymntdate(i)-debut))-1.)*curveN(i)
		sum=sum+(cash_fixe-cash_float)*exp(-rz*paymntdate(i))
	endfor
	return,sum
end

; valeur future d un swap
function Swap_futureValue,k,paymntdate,curveF,curveZ,curveN,spread,t
Return,swap(shift_dt2(k,paymntdate,t),shift_dt1(paymntdate,t), curveF,curveZ,shift_dt2(curveN,paymntdate,t),spread)
end

    
; *******************  M O D E L E   D E   H U L L   A N D    W H I T E  **********************

; fonction qui calcule un arbre de HW  qui matche la courbe de zero coupons rz(i)
; rz contient n+2 element  c est a dire une information sur
; 2 periodes de plus que le theta qui sortira
; afin de matcher les notations de l article de J.Hull(93)
; le modele dynamique des taux court r(i) est:
;  d r(t)=[ theta(t) - a r(t) ] dt + sigma dw(t)
; on suppose que rz1 contient des taux discret. [ rdiscret=exp(rcontinu) - 1 ]
function HWbuild_tree,rz1,a,sigma,delta_t
rz=alog(1.+rz1)
n=n_elements(rz)
m=n+1
  delta_r=sigma*sqrt(3*delta_t)
r=fltarr(n,2*m)
nodes=fltarr(n,2*m)
k=fltarr(n,2*m)
p1=fltarr(n,2*m)
p2=fltarr(n,2*m)
p3=fltarr(n,2*m)
q=fltarr(n,2*m)		;Q(i,j) de l'  article
theta=fltarr(n-1)
r(0,m)=rz(0)
nodes(0,m)=1
q(0,m)=1
theta(0)=(2./delta_t)*rz(1)+sigma^2*delta_t/2.-2*rz(0)/delta_t+a*rz(0) 
mu=theta(0)-a*r(0,m)
k_brut=mu*delta_t/delta_r+m
k(0,m)=int(k_brut+0.5)
nodes(1,k(0,m)-1)=1
nodes(1,k(0,m))=1
nodes(1,k(0,m)+1)=1
eta=mu*delta_t+(m-k(0,m))*delta_r
p1(0,m)=(sigma^2*delta_t+eta^2)/(2*delta_r^2)+eta/(2*delta_r)
p2(0,m)=1.-(sigma^2*delta_t+eta^2)/(delta_r^2)
p3(0,m)=1.-p1(0,m)-p2(0,m)
for i=1,n-2 do begin
	for j=0,2*m-1 do if nodes(i,j) eq 1 then begin
			; calcul de q(i,j) en renversant les fleches k(i-1,*)
			r(i,j)=r(0,m)+(j-m)*delta_r
			sum=0
			for jstar=0,2*m-1 do begin
				case 1 of 
					j eq k(i-1,jstar)+1 : sum=sum+p1(i-1,jstar)*q(i-1,jstar)*exp(-r(i-1,jstar)*delta_t)
					j eq k(i-1,jstar) : sum=sum+p2(i-1,jstar)*q(i-1,jstar)*exp(-r(i-1,jstar)*delta_t)
					j eq k(i-1,jstar)-1 : sum=sum+p3(i-1,jstar)*q(i-1,jstar)*exp(-r(i-1,jstar)*delta_t)
					else:
				endcase
 			endfor
			q(i,j)=sum
	endif
	sum=0
	for j=0,2*m-1 do if nodes(i,j) eq 1 then sum=sum+q(i,j)*exp(-2*r(i,j)*delta_t+a*r(i,j)*delta_t^2)
	theta(i)=(i+2)*rz(i+1)/delta_t+sigma^2*delta_t/2.+alog(sum)/(delta_t^2)
;if i eq 1 then theta(i)=0.0213
	for j=0,2*m-1 do if nodes(i,j) eq 1 then begin
		mu=theta(i)-a*r(i,j)
		k_brut=mu*delta_t/delta_r+j
		k(i,j)=int(k_brut+0.5)
		nodes(i+1,k(i,j)-1)=1
		nodes(i+1,k(i,j))=1
		nodes(i+1,k(i,j)+1)=1
		eta=mu*delta_t+(j-k(i,j))*delta_r 
		p1(i,j)=(sigma^2*delta_t+eta^2)/(2*delta_r^2)+eta/(2*delta_r) 
		p2(i,j)=1.-(sigma^2*delta_t+eta^2)/(delta_r^2)
		p3(i,j)=1.-p1(i,j)-p2(i,j)
	endif
endfor
resultat={,r:r,nodes:nodes,k:k,q:q,p1:p1,p2:p2,p3:p3,theta:theta,delta_t:[delta_t],delta_r:[delta_r],nt:[n],nr:[m]}
;print,nodes
return,resultat
end


; fonction quiretourne les differentes courbe de coupon zero de maturite, 
;  ceci dans l'arbre calcule par la fonction HWbuild_tree
; n2 est le numero de la maturite du zero coupon
; n1 est numero de la date d'evaluation
; le resultat est un tenseur d ordre 3 
function HWcompute_crb,arg
nodes=arg.nodes
p1=arg.p1
p2=arg.p2
p3=arg.p3
k=arg.k
r=arg.r
n=arg.nt
m=arg.nr
delta_t=arg.delta_t
delta_r=arg.delta_r
couponZ=fltarr(n-2,n-2,2*m)  ; 1iere dim : forwardage du cp (il faut rajouter 1 et multiplier par delta_t)
                             ; 2ieme dim : maturite relative du cp (il faut rajouter 1 et multiplier par delta_t)
                             ; 3ieme dim : etat stochastique (meme convention que pour l arbre) 
qq=fltarr(n,2*m)
for n1a=0,n-3 do begin
  n1=n1a+1
  for n2a=0,n-n1-2 do begin
	n2=n1+n2a+1
	qq(*,*)=0.
	qq(n2,*)=1.
	for i=n2-1.,n1,-1 do for j=0,2*m-1 do if nodes(i,j) eq 1 then $
	   qq(i,j)=(exp(-r(i,j)*delta_t))*(p1(i,j)*qq(i+1,k(i,j)+1)+p2(i,j)*qq(i+1,k(i,j))+p3(i,j)*qq(i+1,k(i,j)-1))
        couponZ(n1a,n2a,*)=qq(n1,*)
  endfor
endfor
return,couponZ
end


;construit une structure comprenant l arbre et les courbes de taux en interpollant la courbe de taux
function HWbuild_crbtree,curverz,a,sigma,n
rz=fltarr(n)
zsize=size(curverz)
nz=zsize(1)
delta_t=curverz(nz-1,0)/float(n)
for i=0,n-1 do rz(i)=interpolation(curverz,delta_t*(i+1))
arbre=HWbuild_tree(rz,a,sigma,delta_t)
virtual_discount_crb=HWcompute_crb(arbre)
crbtree={,arbre:arbre,discount_tree:virtual_discount_crb}
return,crbtree
end



;*************************************************************************************************************

;                   F O N C T I O N S  Q U I   S P E C I F I E N T  L A   S T R U C T U R E 

;                                 D U   S W A P    E T    D U  S W A P T I O N 


;*****************************************************************************************************************

;   STRUCTURE DU SWAP

function EVAL_SWAP,curvezero,t
paymntdate=[0.5,1.,1.5,2.,2.5,3.]
tf=0.5
curvez=curvezero
curvef=taux_forward(curvezero,tf)
curveN=[1000,1000,1000,1000,1000,1000]
curvek=[0.12,0.12,0.12,0.12,0.12,0.12]  ;taux fixe
return,swap(shift_dt2(curvek,paymntdate,t),shift_dt1(paymntdate,t), curveF,curveZ,shift_dt2(curveN,paymntdate,t),0)
; on met shift_dt1 et shift_dt2 si on veut que le swap ne glisse pas sinon rien 
end


;  STRUCTURE DU PAYBACK DE L OPTION

function EVAL_PAYBACK,val_swap,t            ;call
if val_swap gt 0 then return,val_swap else return,0
end


;  DESCRIPTION DE  L'EXERCABILITE (CARACTERE AMERICAIN)

function AMERICANP,t
if t gt 1. then return,1 else return,0
end

;******************************************************************************************************************


;***********************       C A L C U L      D E       L ' O P T I O N     *************************************


; fonction qui calcule en tout point d un arbre une valeur de swap sur la courbe zero coupons correspondant au neud
; calcule alors un valeur d exercice grace a cette courbe zero coupons
; et evalue l arbre en supposant l' exercice possible quand AMERICAINP est = 1  jusqua EXPIRATION
; la fonction qui evalue la valeur de swap a partir d une courbe de zero coupons est: EVAL_SWAP
; la fonction qui evalue la valeur de l'option a partir de la valeur du swap  est : EVAL_PAYBACK
;   le mot cle prob_exer stipule que l option est une option dont l exercice est aleatoire conditione par une 
;   probabilte annuelle de prob_exer
;   le mot cle comp_exer pour compulsory exercise  si present signifie que lorsque l'exercice est possible 
;   il est obligatoire.
function HWSwapOption,crbtree,payment_dates,notionals,fixed_rates,type,expiration,american_start,slidingflag
r=crbtree.arbre.r
nodes=crbtree.arbre.nodes
k=crbtree.arbre.k
p1=crbtree.arbre.p1
p2=crbtree.arbre.p2
p3=crbtree.arbre.p3
delta_t=crbtree.arbre.delta_t
delta_r=crbtree.arbre.delta_r
n=crbtree.arbre.nt
m=crbtree.arbre.nr
virtual_discount_crb=crbtree.discount_tree
if keyword_set(spread) eq 1 then spr=1 else spr=0.
imax=int(EXPIRATION/delta_t)+1         ;index temporel (date=imax*delta_t)
support=fltarr(imax,2*m)
opt=fltarr(imax,2*m)
for i=imax-1,0,-1 do begin             ;index pointant dans virtual_taux_crb(i,*,*)
	for j=0,2*m-1 do begin         ; index d etat
		if nodes(i,j) eq 1 then begin
			maturity_max=n-i-2

			curveZ=fltarr(maturity_max,2)
			curveZ(*,0)=(findgen(maturity_max)+1)*delta_t
			for km=0,maturity_max-1 do curveZ(km,1)=virtual_discount_crb(i,km,j)^(-1./curveZ(km,0))-1. 
								  ;index pointant dans virtual_taux_crb(i,km,*)
			if type eq 0 then $
			support(i,j)=max([0,EVAL_SWAP(curveZ,(i+1.)*delta_t)]) $
				else $
			support(i,j)=-min([0,EVAL_SWAP(curveZ,(i+1.)*delta_t)]) 


			if i eq imax-1 then $
					 opt(i,j)=support(i,j)  $
				else  $
					opt(i,j)=exp(-r(i,j)*delta_t)*(p1(i,j)*opt(i+1,k(i,j)+1)+p2(i,j)*opt(i+1,k(i,j))+p3(i,j)*opt(i+1,k(i,j)-1))

			if (american_start ge (i+1)*delta_t)   and (i ne imax-1)  then $
			        opt(i,j)=max([opt(i,j),support(i,j)]) 
			   
			;print,'crb(,',int(i+1),',',j+1,') at time',(i+1.)*delta_t,' k=',k(i,j)
			;print,curvez
			;pprint,'sup=',support(i,j),' opt=',opt(i,j) 

		endif
	endfor
endfor
return,opt(0,m)
end

function HWswapValue,crbtree
virtual_discount_crb=crbtree.discount_tree
maturity_max=crbtree.arbre.nt-2
m=crbtree.arbre.nr
curveZ=fltarr(maturity_max,2)
curveZ(*,0)=(findgen(maturity_max)+1)*crbtree.arbre.delta_t
for km=0,maturity_max-1 do curveZ(km,1)=virtual_discount_crb(0,km,m)^(-1./curveZ(km,0))-1.
return,EVAL_SWAP(curveZ,0.)
end


;*******************************   T  E  S  T  S  *********************************************************


;test de swap
pro test01
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115]])
print,EVAL_SWAP(curverz,1.,0.)
end


;test de HWbuild_tree
pro test02
rz=[0.1,0.105,0.11,0.1125,0.115]
a=0.1
sigma=0.014
delta_t=1.
resultat=HWbuild_tree(rz,a,sigma,delta_t)
print,'resultat=',resultat.theta
end

;test de HWcompute_crb
pro test03
rz=[0.1,0.105,0.11,0.1125,0.115,0.117,0.119,0.122]
a=0.1
sigma=0.014
delta_t=1.
tree=HWbuild_tree(rz,a,sigma,delta_t)
z_etats=HWcompute_crb(tree)
print,z_etats
end


;test de HWSwapOption
pro test04
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115],[6.,0.117],[7.,0.123]])
a=0.1
sigma=0.015
n=15
payment_dates=0
notionals=0
fixed_rates=0
type=0
expiration=2.
american_start=0.
slidingflag=0
crbarbre=HWbuild_crbtree(curverz,a,sigma,n)
print,'prix swaption= ',HWSwapOption(crbarbre,payment_dates,notionals,fixed_rates,type,expiration,american_start,slidingflag)
end



