; DEMO   EQUITY SWAPS

 function call,s,k,t,r,sigm
if (t le 0)  then return,0>(s-k) else begin
	if (sigm eq 0.0) then return,(s-k*exp(-r*t))>0 else begin
		d1=(alog(s/float(k))+(r+sigm*sigm/2)*t)/(sigm*sqrt(t))
		d2=d1-(sigm*sqrt(t))
		return,s*gaussint(d1)-k*gaussint(d2)*exp(-r*t)
	endelse
endelse
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

;********************************Auxilliary calculations***************************************
;calcule la courbe te taux forward  de terme tf a partir d une courbe zero coupons

function taux_forward_swap,curvezero,tf
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


function forward_sliding_swaps,curvezero,t,tf,paymntdate,curveN,curvek
curvef=taux_forward_swap(curvezero,tf)
spread=0
return,swap(curvek,paymntdate, curveF,curvezero,curveN,spread)
; on met shift_dt1 et shift_dt2 si on veut que le swap ne glisse pas sinon rien 
end



function forward_non_sliding_swaps,curvezero,t,tf,paymntdate,curveN,curvek
curvef=taux_forward_swap(curvezero,tf)
spread=0
return,swap(shift_dt2(curvek,paymntdate,t),shift_dt1(paymntdate,t), curveF,curvezero,shift_dt2(curveN,paymntdate,t),spread)
; on met shift_dt1 et shift_dt2 si on veut que le swap ne glisse pas sinon rien 
end


; *******************  M O D E L E   D E   H U L L   A N D    W H I T E  **********************

; fonction qui calcule un arbre de HW  qui matche la courbe de zero coupons rz(i)
; rz contient n+2 element  c est a dire une information sur
; 2 periodes de plus que le theta qui sortira
; afin de matcher les notations de l article de J.Hull(93)
; le modele dynamique des taux court r(i) est:
;  d r(t)=[ theta(t) - a r(t) ] dt + sigma dw(t)
; on suppose que rz1 contient des taux discret. [ rdiscret=exp(rcontinu) - 1 ]
; et que les taux de rz sont regulierement espace en temps separe par dela_t
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
for j=0,2*m-1 do if nodes(n-1,j) eq 1 then begin
r(n-1,j)=r(0,m)+(j-m)*delta_r
endif
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
crbtree={,arbre:arbre,discount_tree:virtual_discount_crb,rvol:sigma,a:a}
return,crbtree
end
;*************************************************************************************************************

;                   Model a deux facteurs  HW  X  LN

;*****************************************************************************************************************



;construit un arbre a deux facteurs , en calculant les probabilite de branchement dans un modele a six points de branchement
;vol : volatility of the index
;sigma : sigma parameter of the HW factor
function HWLNtree,crbtree,vol,sigma,corr,s
r=crbtree.arbre.r
nodes=crbtree.arbre.nodes
k=crbtree.arbre.k
p1=crbtree.arbre.p1
p2=crbtree.arbre.p2
p3=crbtree.arbre.p3
delta_t=crbtree.arbre.delta_t
delta_r=crbtree.arbre.delta_r
theta=crbtree.arbre.theta
n=crbtree.arbre.nt
m=crbtree.arbre.nr
pu1=fltarr(n,2*m)
pd1=pu1
pu2=pu1
pd2=pu1
pu3=pu1
pd3=pu1
for i=0,n-2 do begin
	for j=0,2*m-1 do if nodes(i,j) eq 1 then begin
		a=exp(r(i,j)*delta_t)
		u=exp(vol*sqrt(delta_t))
		d=1/u
		pu=(a-d)/(u-d)
	;problem lie a la positivite de pu et au caractere arbitrage-free du model qui nest pas encore demontre
	;peut etre meme un probleme de convergence vers un model continu quand n -> infini
		if pu lt 0. then pu =0.
		pd=1.-pu
		rk=r(i+1,k(i,j))
		;eta=(p1(i,j)-p3(i,j))*delta_r
		; determination des bornes max de covariances :r1 et r2
		Vcorr=transpose([(rk+delta_r)*S*u,(rk+delta_r)*S*d,(rk-delta_r)*S*u,(rk-delta_r)*S*d,$
		      rk*S*u,rk*S*d])
		alpha=pu
		x=p1(i,j)
		y=p2(i,j)

		px= [p1(i,j)*pu,p1(i,j)*pd,p2(i,j)*pu,p2(i,j)*pd,p3(i,j)*pu,p3(i,j)*pd]
 		v1=[1,-1,0,0,-1,1]
		v2=[0,0,1,-1,-1,1]

;determination du maximum de correlation
	if (1-alpha)*x -alpha*y lt -(1-alpha)*(1-x-y) then begin
		eta1max=(1-alpha)*x
		eta2max=-(1-alpha)*(1-y)
		endif else if alpha*(1-x-y) lt (1-alpha)*x-alpha*y then begin
		eta1max=alpha*(1-x)
		eta2max=-alpha*y
	endif else begin
		eta1max=(1-alpha)*x
		eta2max=-alpha*y
	endelse
	pxmax=px+eta1max*v1+eta2max*v2
;determination du minimum de correlation
	if (1-alpha)*y -alpha*x lt -(1-alpha)*(1-x-y) then begin
		eta1min=-(1-alpha)*(1-x)
		eta2min=(1-alpha)*y
		endif else if alpha*(1-x-y) lt (1-alpha)*y -alpha*x then begin
	eta1min=-alpha*x
	eta2min=alpha*(1-y)
	endif else begin
	eta1min=-alpha*x
	eta2min=(1-alpha)*y
	endelse
	pxmin=px+eta1min*v1+eta2min*v2
	xcorr=[vcorr#pxmin,vcorr#pxmax]
			r1=min(xcorr)
			r2=max(xcorr)
			Prod_Expect=sqrt(r1*r2)
			prod_sigma=alog(r2/r1)/(2*delta_t)
;calcul du jeu de probabilite qui rend la correlation donnee
			RR=Prod_Expect*exp(corr*prod_sigma*delta_t)
	krk=(RR-r1)/(r2-r1)
	eta1=eta1min+krk*(eta1max-eta1min)
	eta2=eta2min+krk*(eta2max-eta2min)
	prob=px+eta1*v1+eta2*v2
	prob=prob>0.
	aa=total(prob)
;	prob=prob/aa
	pu1(i,j)=prob(0)
	pd1(i,j)=prob(1)
	pu2(i,j)=prob(2)
	pd2(i,j)=prob(3)
	pu3(i,j)=prob(4)
	pd3(i,j)=prob(5)
        endif
endfor
resultat={,pu1:pu1,pd1:pd1,pu2:pu2,pd2:pd2,pu3:pu3,pd3:pd3,svol:vol,corr:corr,s:s}
return,resultat
end



;portef est un portefeuille soit une matrice de 80 * n  reels

function HWLNOption,crbtree,hwlnprob,portef,nbinstrument,yieldlist,t_divlist
r=crbtree.arbre.r
nodes=crbtree.arbre.nodes
k=crbtree.arbre.k
p1r=crbtree.arbre.p1
p2r=crbtree.arbre.p2
p3r=crbtree.arbre.p3
delta_t=crbtree.arbre.delta_t
delta_r=crbtree.arbre.delta_r
n=crbtree.arbre.nt
m=crbtree.arbre.nr
rvol=crbtree.rvol
pu1=hwlnprob.pu1
pd1=hwlnprob.pd1
pu2=hwlnprob.pu2
pd2=hwlnprob.pd2
pu3=hwlnprob.pu3
pd3=hwlnprob.pd3
svol=hwlnprob.svol
corr=hwlnprob.corr
s=hwlnprob.s
;sizeportef=size(portef)
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
virtual_discount_crb=crbtree.discount_tree
imaxinst=int(expiration_date/delta_t-0.0000001)+1         ;index temporel (date=imax*delta_t)
imax=max(imaxinst)
if n-imax+1 lt 2 then print,'expiration of option to far compared to the initial interest rate zero coupons curve'
index=fltarr(imax,imax+1)
u=exp(svol*sqrt(delta_t))
d=1/u
for i=0,imax-1 do begin
	indexlist=where(t_divlist lt i*delta_t,count)
	if count eq 0 then pp =1 else begin
		ylist=yieldlist(indexlist)
		pp=1.
		for ii=0,count-1 do pp=pp*(1-ylist(ii))
	endelse
	for j=0,i do index(i,j)=s*u^j*d^(i-j)*pp
endfor
sumopt=fltarr(imax,2*m,imax+1)
opt=fltarr(imax,2*m,imax+1,ninst)
support=fltarr(imax,2*m)
for i=imax-1,0,-1 do begin             ;index pointant dans virtual_taux_crb(i,*,*)

		for jtaux=0,2*m-1 do begin         ; index d etat
		if nodes(i,jtaux) eq 1 then begin
			maturity_max=n-i-2                ;maturite max de la courbe de zero coupon
			curveZ=fltarr(maturity_max,2)
			curveZ(*,0)=(findgen(maturity_max)+1)*delta_t
			for km=0,maturity_max-1 do curveZ(km,1)=virtual_discount_crb(i,km,jtaux)^(-1./curveZ(km,0))-1. 
								  ;index pointant dans virtual_taux_crb(i,km,*)
			for jindex=0,i do begin
			   for iinst=0,ninst-1 do begin
				if i eq imaxinst(iinst)-1 then begin
; ici les conditions terminales pour chaque instruments qui peuvent etre une fonction 
; de la courbe des taux curveZ, du temps : (i+1)*delta_t, de l index : index(i,jindex) suivant le type d'instrument : typeinst(iinst)
; et de la maturite effecive de l instrument expirationdate(iinst)
				case instrument_type(iinst) of 
				 100 : opt(i,jtaux,jindex,iinst)=0
				 200 : begin
					nbp=nb_payments(iinst)-1
					opt(i,jtaux,jindex,iinst)=max($
		[forward_non_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst)$
					,curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end
				 201 : begin
					nbp=nb_payments(iinst)-1
					opt(i,jtaux,jindex,iinst)=max($
		[forward_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst)$
					,curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end

				 250 :  begin
					nbp=nb_payments(iinst)-1
					opt(i,jtaux,jindex,iinst)=-min($
		[forward_non_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst),$
					curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end

				 251 : begin
					nbp=nb_payments(iinst)-1
					opt(i,jtaux,jindex,iinst)=-min($
		[forward_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst),$
					curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end

				 300 : opt(i,jtaux,jindex,iinst)=0
				 400 : opt(i,jtaux,jindex,iinst)=0
				 401 : opt(i,jtaux,jindex,iinst)=0
				 402 : opt(i,jtaux,jindex,iinst)=0
				 403 : opt(i,jtaux,jindex,iinst)=0
				 450 : opt(i,jtaux,jindex,iinst)=0
				 451 : opt(i,jtaux,jindex,iinst)=0
				 452 : opt(i,jtaux,jindex,iinst)=0
				 453 : opt(i,jtaux,jindex,iinst)=0
				 500 : opt(i,jtaux,jindex,iinst)=0

				 600 : opt(i,jtaux,jindex,iinst)=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 601 : opt(i,jtaux,jindex,iinst)=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 602 : opt(i,jtaux,jindex,iinst)=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 603 : opt(i,jtaux,jindex,iinst)=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 650 : opt(i,jtaux,jindex,iinst)=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 651 : opt(i,jtaux,jindex,iinst)=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 652 : opt(i,jtaux,jindex,iinst)=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 653 : opt(i,jtaux,jindex,iinst)=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 700 : opt(i,jtaux,jindex,iinst)=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)

				 701 : opt(i,jtaux,jindex,iinst)=0
				 702 : opt(i,jtaux,jindex,iinst)=0
				 703 : opt(i,jtaux,jindex,iinst)=0
				 750 : opt(i,jtaux,jindex,iinst)=0
				 751 : opt(i,jtaux,jindex,iinst)=0
				 752 : opt(i,jtaux,jindex,iinst)=0
				 753 : opt(i,jtaux,jindex,iinst)=0
				endcase
				endif else begin
				if i lt imax-1 then begin
					opt(i,jtaux,jindex,iinst)=exp(-r(i,jtaux)*delta_t)*($
					pu1(i,jtaux)*opt(i+1,k(i,jtaux)+1,jindex+1,iinst)+$
					pd1(i,jtaux)*opt(i+1,k(i,jtaux)+1,jindex,iinst)+$
					pu2(i,jtaux)*opt(i+1,k(i,jtaux),jindex+1,iinst)+$
					pd2(i,jtaux)*opt(i+1,k(i,jtaux),jindex,iinst)+$
					pu3(i,jtaux)*opt(i+1,k(i,jtaux)-1,jindex+1,iinst)+$
					pd3(i,jtaux)*opt(i+1,k(i,jtaux)-1,jindex,iinst)) 
				endif else opt(i,jtaux,jindex,iinst)=0        ;si on passe pas par une condition terminale mais qu on est encore > intrument
				endelse

				if ((i+1)*delta_t ge american_starting_date(iinst)) and (i le (imaxinst(iinst)-1)) then begin

;correction du a l exercabilite de linstrument , la valeur d exercise peut dependre de 
; de la courbe des taux curveZ, du temps : (i+1)*delta_t, de l index : index(i,jindex) suivant le type d'instrument : typeinst(iinst)
		case instrument_type(iinst) of 
				 100 : exercise=0
				 200 : begin
					nbp=nb_payments(iinst)-1
					exercise=max($
		[forward_non_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst),curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end
				 201 : begin
					nbp=nb_payments(iinst)-1
					exercise=max($
		[forward_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst),curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end

				 250 :  begin
					nbp=nb_payments(iinst)-1
					exercise=-min($
		[forward_non_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst),curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end

				 251 : begin
					nbp=nb_payments(iinst)-1
					exercise=-min($
		[forward_sliding_swaps(curveZ,(i+1.)*delta_t,tenor(iinst),paymntdate(0:nbp,iinst),curveN(0:nbp,iinst),curvek(0:nbp,iinst)),0])
					end

				 300 : exercise=0
				 400 : exercise=0
				 401 : exercise=0
				 402 : exercise=0
				 403 : exercise=0
				 450 : exercise=0
				 451 : exercise=0
				 452 : exercise=0
				 453 : exercise=0
				 500 : exercise=0

				 600 : exercise=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 601 : exercise=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 602 : exercise=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 603 : exercise=max([index(i,jindex)-strike(iinst),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 650 : exercise=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 651 : exercise=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 652 : exercise=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 653 : exercise=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)
				 700 : exercise=max([strike(iinst)-index(i,jindex),0])*contract_size(iinst)*nb_of_contracts(iinst)

				 700 : exercise=0
				 701 : exercise=0
				 702 : exercise=0
				 703 : exercise=0
				 750 : exercise=0
				 751 : exercise=0
				 752 : exercise=0
				 753 : exercise=0
				endcase
				opt(i,jtaux,jindex,iinst)=max([opt(i,jtaux,jindex,iinst),exercise])
				endif
				sumopt(i,jtaux,jindex)=total(opt(i,jtaux,jindex,0:ninst-1))
			   endfor
			endfor

		endif
	endfor
endfor
return,sumopt
end

;test de HWbuild_tree
pro test03
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2      ; volatility of the index
corr=+0.     ;  correlation  interest rate / correlation
s=100.
nbinstrument=5
portef=fltarr(80,nbinstrument)
;portef(1:17,0)=[201,0,2.5,0,0,0,0,0,0,0,0,0,0,0,0,3,1]
;portef(20:22,0)=[1,1.5,2]
;portef(40:42,0)=10000000*[1,1,1]
;portef(60:62,0)=[0.05,0.05,0.5]
portef(1:17,0)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,1)=[600,0,5.,0,0,0,0,200,-1,0,110.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,2)=[600,0,5.,0,0,0,0,200,-1,0,120.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,3)=[650,0,5.,0,0,0,0,200,5,0,80.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,4)=[650,0,5.,0,0,0,0,200,5,0,90.,0,0,0,0,0,0]  ;un call 5 an

yieldlist=[0.00,0.00]
t_divlist=[1.5,2.5]
crbtree=HWbuild_crbtree(rz,a,sigma,nb)
hwlnprob=HWLNtree(crbtree,vol,sigma,corr,s)
print,'arbre construit'
resulta=HWLNOption(crbtree,hwlnprob,portef,nbinstrument,yieldlist,t_divlist)
end


function HWLNTOption,rz,a,sigma,vol,corr,nb,yieldlist,t_divlist,s,option
if n_elements(option) ne 80 then begin
	print,'error in the dimension of the option:',n_elements(option)
	return,0
	endif
crbtree=HWbuild_crbtree(rz,a,sigma,nb)
hwlnprob=HWLNtree(crbtree,vol,sigma,corr,s)
m=crbtree.arbre.nr
portef=fltarr(80,1)
portef(0:79,0)=option
nbinstrument=1
resulta=HWLNOption(crbtree,hwlnprob,portef,nbinstrument,yieldlist,t_divlist)
return,resulta(0,m,0)
end


pro test05
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2
corr=+0.
s=100.
option=fltarr(80)
option(1:17)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
yieldlist=[0.00,0.00]
t_divlist=[1.5,2.5]
for nb=8,30 do print,HWLNTOption(rz,a,sigma,vol,corr,nb,yieldlist,t_divlist,s,option)
end


function impvol_HWLNTOption,rz,a,sigma,vol,corr,nb,yieldlist,t_divlist,s,option
instrument_type=option(1)
option_type=option(2)
expiration_date=option(3)
american_starting_date=option(4)
sliding_flag=option(5)
actual_price=option(6)
maturity=option(7)
contract_size=option(8)
nb_of_contracts=option(9)
american_flag=option(10)
strike=option(11)
underlying_type=option(12)
barrier_type=option(13)
barrier=option(14)
notional=option(15)
nb_payments=option(16)
tenor=option(17)
prime=HWLNTOption(rz,a,sigma,vol,corr,nb,yieldlist,t_divlist,s,option)/(contract_size*nb_of_contracts)
print,'prix bimodel:',prime,'/ prix normal:',$
		ACRYLoption(option_type,underlying_type,s,strike,expiration_date,american_starting_date,0.1,vol,50,yieldlist,t_divlist,ctaux=rz)
voll=index_impvol_spot(prime,option_type,underlying_type,american_flag,s,strike,expiration_date,rz,yieldlist,t_divlist)
return,voll
end


pro test06
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2
corr=+0.
s=100.
option=fltarr(80)
option(1:17)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
yieldlist=[0.00,0.00]
t_divlist=[1.5,2.5]
print,impvol_HWLNTOption(rz,a,sigma,vol,corr,nb,yieldlist,t_divlist,s,option)
end

pro test07
ctaux=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
yieldlist=[0.00,0.00]
t_yieldlist=[1.5,2.5]
optype=0
cashtype=0
s=100.
k=100.
t=5.
tta=0.
vol=0.01
n=50
cvol=transpose([[0.0001,vol],[10.,vol]])
print,'prix=',ACRYLoption(optype,cashtype,s,k,t,tta,0.1,vol,n,yieldlist,t_yieldlist,ctaux=ctaux,cvol=cvol)
end

; *********************bucket analysis dans le model HWLN *************


;on se donne une option , une courbe de taux, une coube de vol 
; on rend un vecteur de derivee pour les taux et 
; un vecteur de derivee pour les vol
;  le grid pour ces vecteur est naturellement le grid de calcul binomial

function bucket_taux_HWLNoption, option,a,sigma,vol,corr,n,yieldlist,t_yieldlist,nb,ctaux,bumptaux=bumptaux
;creation d une courbe de taux equivalente composee de nb buckets
btaux=fltarr(nb,2)
delta_b=t/float(nb)
if keyword_set(bumptaux) eq 0 then bumptaux=0.005
for i=0,nb-1 do begin
	btaux(i,0)=(i+1)*delta_b
	btaux(i,1)=integrale_courbec(ctaux,i*delta_b,(i+1)*delta_b)
end
;calcul des buckets
opt_base=HWLNTOption(btaux,a,sigma,vol,corr,n,yieldlist,t_yieldlist,s,option)
bucket_rate=fltarr(nb,2)
for i=0,nb-1 do begin
	bucket_rate(i,0)=(i+1)*delta_b
	btaux(i,1)=btaux(i,1)+bumptaux
	bucket_rate(i,1)=HWLNTOption(btaux,a,sigma,vol,corr,n,yieldlist,t_yieldlist,s,option)-opt_base
	btaux(i,1)=btaux(i,1)-bumptaux
end
return,bucket_rate
end


pro test09
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2
corr=+0.
s=100.
option=fltarr(80)
option(1:17)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
yieldlist=[0.00,0.00]
t_yieldlist=[1.5,2.5]
print,bucket_taux_HWLNoption, option,a,sigma,vol,corr,n,yieldlist,t_yieldlist,nb,ctaux,bumptaux=bumptaux
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



;******************************************************************************************************************


;***********************  visualisation du P&L du portefeuille          *************************************


pro PL_HWLN_visualisation,resulta,crbtree,hwlnprob,rz,s,t

; structure de crbtree
; {,arbre:arbre,discount_tree:virtual_discount_crb,rvol:sigma}
; structure de arbre
; {,r:r,nodes:nodes,k:k,q:q,p1:p1,p2:p2,p3:p3,theta:theta,delta_t:[delta_t],delta_r:[delta_r],nt:[n],nr:[m]}
; structure de hwlnprob
; {,pu1:pu1,pd1:pd1,pu2:pu2,pd2:pd2,pu3:pu3,pd3:pd3,svol:vol,corr:corr,s:s}
; n nombre de pas de la'arbre
; maturite maximale de tout les instruments: (n-2)*deltat_t
delta_t=crbtree.arbre.delta_t
nf=int((t+0.00001)/(delta_t))
rsize=size(resulta)
ntaux=rsize(3)
nindex=rsize(2)
portf_value=fltarr(nindex,ntaux)
for i=0,nindex-1 do portf_value(i,*)=resulta(nf,i,0:ntaux-1)
window,0,title='P&L on date = '+string(nf*delta_t)
surface,portf_value,ytitle='index',xtitle='short rate'
end

 function HWLN_PortefValue,resulta,crbtree,t
delta_t=crbtree.arbre.delta_t
nf=int((t+0.00001)/(delta_t))
rsize=size(resulta)
ntaux=rsize(3)
nindex=rsize(2)
portf_value=fltarr(nindex,ntaux)
for i=0,nindex-1 do portf_value(i,*)=resulta(nf,i,0:ntaux-1)
return,portf_value
end

pro test04
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2
corr=+0.
s=100.
nbinstrument=5
portef=fltarr(80,nbinstrument)
;portef(1:17,0)=[201,0,2.5,0,0,0,0,0,0,0,0,0,0,0,0,3,1]
;portef(20:22,0)=[1,1.5,2]
;portef(40:42,0)=10000000*[1,1,1]
;portef(60:62,0)=[0.05,0.05,0.5]
portef(1:17,0)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,1)=[600,0,5.,0,0,0,0,200,-1,0,110.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,2)=[600,0,5.,0,0,0,0,200,-1,0,120.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,3)=[650,0,5.,0,0,0,0,200,5,0,80.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,4)=[650,0,5.,0,0,0,0,200,5,0,90.,0,0,0,0,0,0]  ;un call 5 an
yieldlist=[0.00,0.00]
t_divlist=[1.5,2.5]
crbtree=HWbuild_crbtree(rz,a,sigma,nb)
hwlnprob=HWLNtree(crbtree,vol,sigma,corr,s)
print,'arbre construit'
t=4.
resulta=HWLNOption(crbtree,hwlnprob,portef,nbinstrument,yieldlist,t_divlist)
PL_HWLN_visualisation,resulta,crbtree,hwlnprob,rz,s,t
end
 


;************************************* visualisation du P&L probabiliste


;calcul de la variance d'une variable lognormale

common aux_lognormal_integration,mu,sigma

function aux1_lognormal_variance,x
common aux_lognormal_integration,mu,sigma
return,x^2*exp(-((x-mu)/sigma)^2/2.-2*x)/(2.5066282746*sigma)
end

function aux2_lognormal_variance,x
common aux_lognormal_integration,mu,sigma
return,x^3*exp(-((x-mu)/sigma)^2/2.-2*x)/(2.5066282746*sigma)
end

function lognormal_mean_variance,m,s
;return  [mean,variance]
n=7
common aux_lognormal_integration,mu,sigma
mu=m
sigma=s
ey=intfcn("aux1_lognormal_variance",mu-n*sigma,mu+n*sigma,rule=6)
ey2=intfcn("aux2_lognormal_variance",mu-n*sigma,mu+n*sigma,rule=6)
return,[ey,ey2-ey^2]
end

pro testint,n
print,lognormal_mean_variance(0.1,0.2,n)
end


; calcul de la densite binormale
common aux_binormal_density,sigma1,sigma2,rho,mu1,mu2

function binormal_density,sigma1,sigma2,rho,mu1,mu2,x1,x2
f=1/(6.28318530718*sigma1*sigma2*sqrt(1-rho^2))*$
	exp(-1/(2*(1-rho^2))*(((x1-mu1)/sigma1)^2-2*rho*(x1-mu1)/sigma1*(x2-mu2)/sigma2+((x2-mu2)/sigma2)^2))
return,f
	
end


function aux1_binormal_density,x
common aux_binormal_density,sigma1,sigma2,rho,mu1,mu2
return,binormal_density(sigma1,sigma2,rho,mu1,mu2,x(0),x(1))
end

;  calcul de la densite binormale associee au pavage pour une matrice
 function  binormal_density_tree,n1,delta_x1,n2,delta_x2,sigma1z,sigma2z,rhoz,mu1z,mu2z,x01,x02
common aux_binormal_density,sigma1,sigma2,rho,mu1,mu2
sigma1=sigma1z
sigma2=sigma2z
rho=rhoz
mu1=mu1z
mu2=mu2z
;return an n1 X n2 arrays
; les domaines associe au couple (i,j) sont  :
;              [x01+(i-0.5)*delta_x1,x01+(i+0.5)*delta_x1] X [x02+(j-0.5)*delta_x2,x01+(j+0.5)*delta_x2]
;
result=fltarr(n1,n2)
for i=0,n1-1 do for j=0,n2-1 do begin
ip=i-n1/2.+0.5
jp=j-n2/2.+0.5
result(i,j)= $
  intfcnhyper("aux1_binormal_density", $
              [x01+(ip-0.5)*delta_x1, x02+(jp-0.5)*delta_x2], $
              [x01+(ip+0.5)*delta_x1, x02+(jp+0.5)*delta_x2], $
              Max_Evals = 10000)
endfor
return,result
end

pro testbin
n1=20
delta_x1=0.022
n2=20
delta_x2=0.022
sigma1z=0.1
sigma2z=0.1
rhoz=0.7
mu1z=0.2
mu2z=0.2
x01=0.2
x02=0.2
f=binormal_density_tree(n1,delta_x1,n2,delta_x2,sigma1z,sigma2z,rhoz,mu1z,mu2z,x01,x02)
shade_surf,f
end


; calcul des probabilite forward absolue associe a un model HW X LN  sur un arbre 
function HWLN_probability,crbtree,hwlnprob,valorisation,rz,t

; structure de crbtree
; {,arbre:arbre,discount_tree:virtual_discount_crb,rvol:sigma,a:a}
; structure de arbre
; {,r:r,nodes:nodes,k:k,q:q,p1:p1,p2:p2,p3:p3,theta:theta,delta_t:[delta_t],delta_r:[delta_r],nt:[n],nr:[m]}
; structure de hwlnprob
; {,pu1:pu1,pd1:pd1,pu2:pu2,pd2:pd2,pu3:pu3,pd3:pd3,svol:vol,corr:corr,s:s}
; n nombre de pas de la'arbre
; maturite maximale de tout les instruments: (n-2)*deltat_t
delta_t=crbtree.arbre.delta_t
nf=int((t+0.00001)/(delta_t))
ntaux=2*crbtree.arbre.nr
rsize=size(valorisation)
nindex=rsize(3)
sigmar=crbtree.rvol
r=crbtree.arbre.r
nodes=crbtree.arbre.nodes
rmoy=0.
nnodes=0
for j=0,crbtree.arbre.nr-1 do if nodes (nf,j) eq 1 then begin
	nnodes=nnodes+1
	rmoy=rmoy+r(nf,j)
endif
x0taux=rmoy/nnodes
theta=crbtree.arbre.theta
index0=hwlnprob.s
volindex=hwlnprob.svol
a=crbtree.a
x0index=index0
delta_taux=crbtree.arbre.delta_r
delta_index=(exp(volindex*sqrt(nf*delta_t))-1.)*index0
sigmataux=sigmar*sqrt((exp(a*t)-1)/(2*a))
res=lognormal_mean_variance((1+interpolation(rz,nf*delta_t))^t,volindex)
sigmaindex=res(1)*index0
rho=hwlnprob.corr
mutaux=x0taux
muindex=res(0)*index0
;print,'ntaux,delta_taux,nindex,delta_index,sigmataux,sigmaindex,rho,mutaux,muindex,x0taux,x0index'
;print,ntaux,delta_taux,nindex,delta_index,sigmataux,sigmaindex,rho,mutaux,muindex,x0taux,x0index
f=binormal_density_tree(ntaux,delta_taux,nindex,delta_index,sigmataux,sigmaindex,rho,mutaux,muindex,x0taux,x0index)
return,f
end





;il faut traiter la dimension de de lindex en integrant entre des x(j) distant de s0[u^(j+1)-u^j] dans le calcul de la matrice de probabilite
;on ne tient pas compte du mean reversion pour la determination du mutaux ,??
;il faut verifier l influence du yield_list dans la determination de l'arbre croise hwlnprob .


; obtention de la courbe de distribution du P&L
function profit_loss_probability,portef_value,probability
sizez=size(portef_value)
n=sizez(1)
m=sizez(2)
p=m*n
;linearisation des matrices
val=fltarr(p)
prob=fltarr(p)
for i=0,n-1 do for j=0,m-1 do begin
	val(i*m+j)=portef_value(i,j)
	prob(i*m+j)=probability(i,j)
endfor
indexval=sort(val)
val1=val(indexval)
prob1=prob(indexval)
y=fltarr(p+1)
x=fltarr(p+1)
x(0)=0
y(0)=0
for i=1,p do begin
	x(i)=x(i-1)+prob1(i-1)
	y(i)=y(i-1)+val1(i-1)
endfor
return,{,probability:x,profit:y}
end

function distribution_to_density,distrib
prob=distrib.probability
profit=distrib.profit
dens=deriv(profit,prob)
return,{,density:dens,profit:profit}
end

function distribution_to_density1,distrib
prob=distrib.probability
profit=distrib.profit
n=n_elements(prob)
dens=fltarr(n)
for i=1,n-1 do dens(i)=prob(i)-prob(i-1)
return,{,density:dens,profit:profit}
end


function distribution_to_risk,distrib
prob=distrib.probability
profit=distrib.profit
n=n_elements(prob)
dens=fltarr(n)
for i=1,n-1 do dens(i)=prob(i)-prob(i-1)
avgprice=fltarr(n)
for i=0,n-1 do avgprice(i)=dens(i)*profit(i)
price=total(avgprice)
risk=fltarr(n)
for i=0,n-1 do risk(i)=price-total(avgprice(0:i))
return,{,risk:risk,probability:prob}
end


pro testpl
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2
corr=+0.
s=100.
nbinstrument=5
portef=fltarr(80,nbinstrument)
;portef(1:17,0)=[201,0,2.5,0,0,0,0,0,0,0,0,0,0,0,0,3,1]
;portef(20:22,0)=[1,1.5,2]
;portef(40:42,0)=10000000*[1,1,1]
;portef(60:62,0)=[0.05,0.05,0.5]
portef(1:17,0)=[600,0,5.,0,0,0,0,200,-1,0,110.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,1)=[600,0,5.,0,0,0,0,200,-1,0,120.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,2)=[650,1,5.,0,0,0,0,200,1,0,80.,0,0,0,0,0,0]  ;un put 5 an
portef(1:17,3)=[650,1,5.,0,0,0,0,200,5,0,90.,0,0,0,0,0,0]  ;un put 5 an
portef(1:17,4)=[600,1,5.,0,0,0,0,200,5,0,110.,0,0,0,0,0,0]  ;un call 5 an
yieldlist=[0.00,0.00]
t_divlist=[1.5,2.5]
t=4.

;#######construction des arbres
crbtree=HWbuild_crbtree(rz,a,sigma,nb)			; arbre de taux
hwlnprob=HWLNtree(crbtree,vol,sigma,corr,s)			; arbre unifie index-taux
valorisation=HWLNOption(crbtree,hwlnprob,portef,nbinstrument,yieldlist,t_divlist)  ;valorisation d un portefeuille
print,'arbres construit'

;#######construction du risque probabiliste
portef_value=HWLN_portefvalue(valorisation,crbtree,t)               ; future dependnace du portefeuill aux facteurs
probability=HWLN_probability(crbtree,hwlnprob,valorisation,rz,t)     ;  implicit forward probability
pnl=profit_loss_probability(portef_value,probability)			;  distribution du P & L
;density=distribution_to_density(pnl)			;  density of the P & L
density1=distribution_to_density1(pnl)			;  density of the P & L  other definition
risk=distribution_to_risk(pnl)			; risk of the P & L

; #### display
window,0,title='P&L on date = '+string(t)
surface,portef_value,ytitle='index',xtitle='short rate'
window,1,title='probability'
surface,probability,ytitle='index',xtitle='short rate'
window,2,title='P&L probability distribution'
plot,pnl.profit,pnl.probability
;window,3,title='P&L probability density'
;plot,density.profit,density.density
window,5,title='P&L probability density 1'
plot,density1.profit,density1.density
window,4,title='P&L Probabilistic Risk'
plot,risk.probability,risk.risk

end
