  ;                     ***************************** C A P S *******************
;calcul d un cap  
; paymntdate = liste des date de paiement 
; k = strike
; curveF= courbe de taux forward
; curveZ= courbe de taux zero coupons


function cap,k,paymntdate,curveF,curveZ,vol
	sum=0
	for i=0,n_elements(paymntdate)-1 do begin
		if i eq 0 then debut=0 else debut=paymntdate(i-1)
		s0=interpolation(curveF,debut)
		rz=interpolation(curveZ,paymntdate(i))
		caplet=(paymntdate(i)-debut)*callF(s0,k,paymntdate(i)-debut,0,vol)$
						*exp(-rz*paymntdate(i))
;print,'let=,discount=',callF(s0,k,paymntdate(i)-debut,0,vol),exp(-rz*paymntdate(i))
		sum=sum+caplet
	endfor
	return,sum
	end

function floor,k,paymntdate,curveF,curveZ,vol
	sum=0
	for i=0,n_elements(paymntdate)-1 do begin
		if i eq 0 then debut=0 else debut=paymntdate(i-1)
		s0=interpolation(curveF,debut)
		rz=interpolation(curveZ,paymntdate(i))
		caplet=(paymntdate(i)-debut)*putF(s0,k,paymntdate(i)-debut,0,vol)$
						*exp(-rz*paymntdate(i))
;print,'let=,discount=',putF(s0,k,paymntdate(i)-debut,0,vol),exp(-rz*paymntdate(i))
		sum=sum+caplet
	endfor
	return,sum
	end




;                 ***************************** S W A P S **********************
;valeur d un swap
; le side 1 est le taux float donne en interpolant la courbe curveF (forward rate)
; le side 2 est le taux fixe k
; paymntdate est la liste des dates de payment
; curveZ est la coube des zero coupons
; curveN est la liste des notionnels applicables (meme longueur que paymntdate)
; les taux sont suposes s appliquer exponentiellement sur les durees exactes

function swap,k,paymntdate,curveF,curveZ,curveN
	sum=0
	for i=0,n_elements(paymntdate)-1 do begin
		if i eq 0 then debut=double(0.000001) else debut=paymntdate(i-1)
		rf=interpolation(curveF,paymntdate(i))
		;rf=moyenne_courbe2(curveF,debut,paymntdate(i))     ;si on veut moyenner le taux forw. sur la periode de paiement
		rz=interpolation(curveZ,paymntdate(i))*paymntdate(i)
		cash_float=(exp(rf*(paymntdate(i)-debut))-1.)*curveN(i)
		cash_fixe=(exp(k*(paymntdate(i)-debut))-1.)*curveN(i)
;print,'fixe,flt,disc fct',cash_fixe,cash_float,exp(-rz*paymntdate(i))
		sum=sum+(cash_fixe-cash_float)*exp(-rz*paymntdate(i))
	endfor
;print,'prix du swap',sum
	return,sum
end

;calcul du prix futur d un swap
function forward_swap_price,t,k,paymntdate,curveF,curveZ,curveN
;print,'py-dt',shift_dt1(paymntdate,t)
;print,'f',shift_courbe(curveF,t)
;print,'z',shift_courbe(curveZ,t)
;print,'n',shift_dt2(curveN,paymntdate,t)
return,swap(k,shift_dt1(paymntdate,t), $
		shift_courbe(curveF,t), $
		shift_courbe(curveZ,t), $
		shift_dt2(curveN,paymntdate,t))
end



 
;               ******************* S W A P T I O N S   E U R O P E E N S **********************
; vol est la volatilite du bond forward
function swaptioncall,te,ks,k,paymntdate,curveF,curveZ,curveN,vol
	forward_swap=forward_swap_price(te,k,paymntdate,curveF,curveZ,curveN)
	rz=interpolation(curveZ,te)
	Notionnal=interpolation1(paymntdate,curveN,te)
	forward_bond=Notionnal+forward_swap
	zero_float=transpose([[0,0],[100,0]])
	delt_k=forward_swap_price(te,k-ks,paymntdate,zero_float,curveZ,curveN)
	forward_k=Notionnal+delt_k
;print,'f,k,t,r,v=',forward_bond,forward_k,te,rz,vol
	return,callF(forward_bond,forward_k,te,rz,vol)
end

function swaptionput,te,ks,k,paymntdate,curveF,curveZ,curveN,vol
	forward_swap=forward_swap_price(te,k,paymntdate,curveF,curveZ,curveN)
	rz=interpolation(curveZ,te)
	Notionnal=interpolation1(paymntdate,curveN,te)
	forward_bond=Notionnal+forward_swap
	zero_float=transpose([[0,0],[100,0]])
	delt_k=forward_swap_price(te,k-ks,paymntdate,zero_float,curveZ,curveN)
	forward_k=Notionnal+delt_k
;print,'f,k,t,r,v=',forward_bond,forward_k,te,rz,vol
	return,putF(forward_bond,forward_k,te,rz,vol)
end


    
;               ******************* S W A P T I O N S   A M E R I C A I N S **********************

; fonction qui calcule la courbe implicite theta(i)  qui matche la courbe de zero coupons rz(i)
; rz contient n+2 element  c est a dire une information sur
; 2 periodes de plus que le theta qui sortira
; afin de matcher les notations de l article de J.Hull(93)
; le modele dynamique des taux court r(i) est:
;  d r(t)=[ theta(t) - a r(t) ] dt + sigma dw(t)
function HWbuild_tree,rz,a,sigma,delta_t
n=n_elements(rz)
m=n+1
  delta_r=sigma*sqrt(3*delta_t)
r=fltarr(n,2*m)
nodes=fltarr(n,2*m)
k=fltarr(n,2*m)
p1=fltarr(n,2*m)
p2=fltarr(n,2*m)
p3=fltarr(n,2*m)
q=fltarr(n,2*m)
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
resultat={,r:r,nodes:nodes,k:k,p1:p1,p2:p2,p3:p3,theta:theta,delta_t:[delta_t],delta_r:[delta_r],nt:[n],nr:[m]}
;print,nodes
return,resultat
end

;
; delta_t variable !!!
; mais cette fonction conduit a des probabilite negatives...
; il faut mieux utiliser le truc de HW en choisissant deltat2 un multiple de deltat1 ...
; fonction qui calcule la courbe implicite theta(i)  qui matche la courbe de zero coupons rz(i)
; rz contient n+2 element  c est a dire une information sur
; 2 periodes de plus que le theta qui sortira
; afin de matcher les notations de l article de J.Hull(93)
; le modele dynamique des taux court r(i) est:
;  d r(t)=[ theta(t) - a r(t) ] dt + sigma dw(t)
function VHWbuild_tree,curverz,a,sigma
timez=curverz(*,0)
rz=curverz(*,1)
n=n_elements(timez)
delta_t=fltarr(n)
delta_r=fltarr(n)
for i=1,n-1  do delta_t(i)=timez(i)-timez(i-1)
delta_t(0)=timez(0)
m=n+5
delta_r(0)=sigma*sqrt(3*delta_t(0))
r=fltarr(n,2*m)
p1=fltarr(n,2*m)
p2=fltarr(n,2*m)
p3=fltarr(n,2*m)
q=fltarr(n,2*m)
rm=fltarr(n)
theta=fltarr(n-1)
r(0,m)=rz(0)
q(0,m)=1
rm(0)=rz(0)
theta(0)=(2./delta_t(0))*rz(1)+sigma^2*delta_t(0)/2.-2*rz(0)/delta_t(0)+a*rz(0)
rm(1)=rm(0)*(1-a*delta_t(0))+theta(0)*delta_t(0)
delta_r(1)=delta_r(0)*(1-a*delta_t(0)) 
mu=theta(0)-a*r(0,m)
p1(0,m)=(sigma^2*delta_t(0))/(2*delta_r(0)^2)
p2(0,m)=1.-(sigma^2*delta_t(0))/(delta_r(0)^2)
p3(0,m)=1.-p1(0,m)-p2(0,m)
for i=1,n-2 do begin
	for j=m-i,m+i do  begin
			r(i,j)=rm(i)+(j-m)*delta_r(i) 
			; calcul de q(i,j) en renversant les fleches k(i-1,*)
			q(i,j)=p1(i-1,j-1)*q(i-1,j-1)*exp(-r(i-1,j-1)*delta_t(i))+$
				p2(i-1,j)*q(i-1,j)*exp(-r(i-1,j)*delta_t(i))+$
				p3(i-1,j+1)*q(i-1,j+1)*exp(-r(i-1,j+1)*delta_t(i))
	endfor
	sum=0
	for j=m-i,m+i do  sum=sum+q(i,j)*exp(-2*r(i,j)*delta_t(i)+a*r(i,j)*delta_t(i)^2)
	theta(i)=(i+2)*rz(i+1)/delta_t(i)+sigma^2*delta_t(i)/2.+alog(sum)/(delta_t(i)^2) 
	rm(i+1)=rm(i)*(1-a*delta_t(i))+theta(i)*delta_t(i)
	delta_r(i+1)=delta_r(i)*(1-a*delta_t(i)) 
print,'delta+r=',delta_r(i+1),' inclu dans [',sigma*sqrt(3*delta_t(i+1))/2.,',',2*sigma*sqrt(delta_t(i+1)),']'
	for j=m-i,m+i do begin
		mu=theta(i)-a*r(i,j)
		p1(i,j)=(sigma^2*delta_t(i))/(2*delta_r(i)^2) 
		p2(i,j)=1.-(sigma^2*delta_t(i))/(delta_r(i)^2)
		p3(i,j)=1.-p1(i,j)-p2(i,j)
print,p1(i,j),p2(i,j),p3(i,j)
	endfor
endfor
resultat={,r:r,p1:p1,p2:p2,p3:p3,theta:theta,delta_t:delta_t,delta_r:delta_r,nt:[n],nr:[m]}
return,resultat
end

;construit un arbre en interpollant la courbe de taux
function NHWbuild_tree,curverz,a,sigma,n
rz=fltarr(n)
zsize=size(curverz)
nz=zsize(1)
delta_t=curverz(nz-1,0)/float(n)
for i=0,n-1 do rz(i)=interpolation(curverz,delta_t*(i+1))
arbre=HWbuild_tree(rz,a,sigma,delta_t)
return,arbre
end





; fonction qui calcule les differents etat possible du coupon zero de maturite fixe n2, 
; a la date n1, ceci dans l'arbre calcule par la fonction HWbuild_tree
; n2 est le numero de la maturite du zero coupon
; n1 est numero de la date d'evaluation
; le resultat est un vecteur 
function HWcompute_Z,arg,n1,n2
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
qq=fltarr(n,2*m)
qq(n2,*)=1.
for i=n2-1.,n1,-1 do for j=0,2*m-1 do if nodes(i,j) eq 1 then $
	 qq(i,j)=(exp(-r(i,j)*delta_t))*(p1(i,j)*qq(i+1,k(i,j)+1)+p2(i,j)*qq(i+1,k(i,j))+p3(i,j)*qq(i+1,k(i,j)-1))
  return,qq(n1,*)
end



;test de la construction des coupons zer0

pro test02a
rz=[0.1,0.105,0.11,0.1125,0.115]
a=0.1
sigma=0.014
delta_t=1.
arbre=HWbuild_tree(rz,a,sigma,delta_t)
z_etats=HWcompute_Z(arbre,2,5)
print,z_etats
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


pro test02bb
rz=[0.1,0.105,0.11,0.1125,0.115,0.117,0.119,0.122]
a=0.1
sigma=0.014
delta_t=1.
arbre=HWbuild_tree(rz,a,sigma,delta_t)
z_etats=HWcompute_crb(arbre)
print,z_etats
end




;calcule la courbe te taux forward  de terme tf a partie d une courbe zero coupons
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



function EVAL_SWAP,curvezero,t
paymntdate=[0.5,1.,1.5,2.,2.5,3.]
tf=0.5
curvez=curvezero
curvef=taux_forward(curvezero,tf)
curveN=[1000,1000,1000,1000,1000,1000]
k=0.12   ;taux fixe
return,swap(k,shift_dt1(paymntdate,t), curveF,curveZ,shift_dt2(curveN,paymntdate,t))
; on met shift_dt1 et shift_dt2 si on veut que le swap glisse sinon rien 
end

function EVAL_PAYBACK,val_swap,t            ;call
if val_swap gt 0 then return,val_swap else return,0
end

function AMERICANP,t
if t gt 1. then return,1 else return,0
end

pro test_uo
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115]])
print,EVAL_SWAP(curverz,1.)
end


;fonction qui calcule en tout point d un arbre une valeur de swap sur la courbe zero coupons correspondant au neud
; calcule alors un valeur d exercice grace a cette courbe zero coupons
; et evalue l arbre en supposant l' exercice possible a partir de DEBUTAMERICAIN jusqua EXPIRATION
;la fonction qui evalue la valeur de swap a partir d une courbe de zero coupons est: EVAL_SWAP
; la fonction qui evalue la valeur de l'option a partir de la valeur du swap  est : EVAL_PAYBACK
; 
function HWoption,curvezero,EXPIRATION,a,sigma,n
arbre=NHWbuild_tree(curvezero,a,sigma,n)
; r:r,nodes:nodes,k:k,p1:p1,p2:p2,p3:p3,theta:theta,delta_t:[delta_t],delta_r:[delta_r],nt:[n],nr:[m]
r=arbre.r
nodes=arbre.nodes
k=arbre.k
p1=arbre.p1
p2=arbre.p2
p3=arbre.p3
delta_t=arbre.delta_t
delta_r=arbre.delta_r
n=arbre.nt
m=arbre.nr
imax=int(EXPIRATION/delta_t)+1         ;index temporel (date=imax*delta_t)
virtual_discount_crb=HWcompute_crb(arbre)
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
					support(i,j)=EVAL_PAYBACK(EVAL_SWAP(curveZ,(i+1.)*delta_t),(i+1.)*delta_t)
			if i eq imax-1 then opt(i,j)=support(i,j) else  opt(i,j)=exp(-r(i,j)*delta_t)*($
							p1(i,j)*opt(i+1,k(i,j)+1)+$
							p2(i,j)*opt(i+1,k(i,j))+$
							p3(i,j)*opt(i+1,k(i,j)-1))
			if AMERICANP((i+1)*delta_t) eq 1 then opt(i,j)=max([opt(i,j),support(i,j)])

			;print,'crb(,',int(i+1),',',j+1,') at time',(i+1.)*delta_t,' k=',k(i,j)
			;print,curvez
			;pprint,'sup=',support(i,j),' opt=',opt(i,j) 

		endif
	endfor
endfor

val=opt(0,m)
return,val
end
 

pro testopt
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115],[6.,0.117],[7.,0.123]])
EXPIRATION=2.
a=0.1
sigma=0.015
n=20
op=HWoption(curverz,EXPIRATION,a,sigma,n)
print,op
end

pro test0f
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115]])
m=EVAL_SWAP(curverz,0.6)
print,m
end

;test de HWbuild_tree
pro test0
rz=[0.1,0.105,0.11,0.1125,0.115]
a=0.1
sigma=0.014
delta_t=1.
resultat=HWbuild_tree(rz,a,sigma,delta_t)
print,'resultat=',resultat.theta
end

;test de VHWbuild_tree
pro test00
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115]])
a=0.1
sigma=0.014
resultat=VHWbuild_tree(curverz,a,sigma)
print,'resultat=',resultat.theta
end
 
;test de VHWbuild_tree
pro test01
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[5.,0.1125],[6,0.115]])
a=0.1
sigma=0.014
resultat=VHWbuild_tree(curverz,a,sigma)
print,'resultat=',resultat.theta
end

;test de NHWbuild_tree
pro test02
curverz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5,0.115]])
a=0.1
sigma=0.014
resultat=NHWbuild_tree(curverz,a,sigma,20)
print,'resultat=',resultat.theta
end



;test de HWbuild_tree
pro testa
rz=[0.1,0.105,0.11,0.1125,0.115]
a=0.1
sigma=0.014
delta_t=1.
result=HWbuild_tree(rz,a,sigma,delta_t)
theta=result.theta
x=(findgen(n_elements(rz)-1)+1)*delta_t
rz11=rz
rzb=[0.1,0.105,0.11,0.1125,0.115]
total_time=5.
nb=7
a=0.1
sigma=0.014
n=n_elements(rzb)
rzb1=fltarr(n+1)
rzb1(1:n)=rzb
rzb1(0)=rzb(0)-(rzb(1)-rzb(0))
rz1=interpol(rzb1,nb+1)
rz=rz1(1:nb)
delta_t=total_time/float(nb)
result=HWbuild_tree(rz,a,sigma,delta_t)
theta2=result.theta
x2=(findgen(n_elements(rz)-1)+1)*delta_t
window,1,title='theta'
plot,x,theta
oplot,x2,theta2,linestyle=1
window,0,title='yield curve'
plot,x,rz11
oplot,x2,rz,linestyle=1
end

pro test03
curverz=transpose([[0.001,0.04],[1.,0.041],[5.,0.042]])
timez=curverz(*,0)
rz=curverz(*,1)
n=n_elements(timez)
delta_t=fltarr(n)
for i=1,n-1 do delta_t(i)=timez(i)-timez(i-1)
delta_t(0)=timez(0)
print,timez
print,delta_t
end


;test du swap (forward price)
pro test1
paymntdate=[0.5,1.,1.5,2.]
curvef=transpose([[0,0.04],[1.,0.041],[5.,0.042]])
curvez=transpose([[0,0.04],[1.,0.041],[5.,0.042]])
curveN=[100,100,100,100,100,100,100,100]
k=0.07
print,'forward swap price =',forward_swap_price(1.,k,paymntdate,curveF,curveZ,curveN)
end

;test des caps et des floors
pro test2
paymntdate=[0.5,1.,1.5,2.,2.5,3.,3.5,4.]
curvef=transpose([[0,0.04],[1.,0.041],[5.,0.042]])
curvez=transpose([[0,0.04],[1.,0.041],[5.,0.042]])
k=0.08
vol=0.15
print,'cap=',cap(k,paymntdate,curveF,curveZ,vol)
print,'floor=' ,floor(k,paymntdate,curveF,curveZ,vol)
end

;test des swaptions
pro test3
paymntdate=[0.5,1.,1.5,2.,2.5,3.,3.5,4.]
curvef=transpose([[0,0.04],[1.,0.041],[5.,0.042]])
curvez=transpose([[0,0.04],[1.,0.041],[5.,0.042]])
curveN=[100,100,100,100,100,100,100,100]
te=.6
ks=0.042
k=ks
vol=0.05
print,'swaptioncall=',swaptioncall(te,ks,k,paymntdate,curveF,curveZ,curveN,vol)
print,'swaptionput=',swaptionput(te,ks,k,paymntdate,curveF,curveZ,curveN,vol)
end

