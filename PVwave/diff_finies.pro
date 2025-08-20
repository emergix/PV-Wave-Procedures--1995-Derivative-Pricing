
; on resoud le systeme 
;   T X = G
; T=T(a,b,c)
; a,b,c sont des reels positifs 
; G(i) sert a stocker un valeur intermediaire de dimension egale a celle de G ;
;  la fonction retourne X

function solve_matrice,xa,xb,xc,F,G
m=n_elements(G)
if m ne N_elements(F) then print,'dimension de F et G non correspondante' else begin
	b=replicate(xb,m)
	b(0)=b(0)+xa
	b(m-1)=b(m-1)+xc
	gp=fltarr(m)
	bp=fltarr(m)
;remontee
	bp(m-1)=b(m-1)
	gp(m-1)=g(m-1)
	for i=m-2,0,-1 do begin
		bp(i)=b(i)-(xc*xa/bp(i+1))
		gp(i)=g(i)-(xc*gp(i+1)/bp(i+1))
	endfor
	x=fltarr(m)
;descente
	x(0)=gp(0)/bp(0)
	for i=1,m-1 do begin
		x(i)=(gp(i)-(xa*x(i-1)))/bp(i)
	endfor
	return,x
endelse
end


; on resoud le systeme 
;   T X >= G
;   X >= 0
;   (T X - G , X - F) = 0
; T=T(a,b,c)
; a,b,c sont des reels positifs 
; F(i) est la valeur intrinseque de l'option si le support vaut S(i)
; G(i) sert a stocker un valeur intermediaire de dimension egale a celle de G ;
;  la fonction retourne X

function solve_matrice_sup,a,b,c,F,G
m=n_elements(G)
if m ne N_elements(F) then print,'dimension de F et G non correspondante' else begin
	b2=replicate(b,m)
	b2(0)=b2(0)+a
	b2(m-1)=b2(m-1)+c
	g2=g
	b2p=b2
	g2p=g2
	for i=m-2,0,-1 do begin
		b2p(i)=b2(i)-(c*a/b2p(i+1))
		g2p(i)=g2(i)-(c*g2p(i+1)/b2p(i+1))
	endfor
	x=fltarr(m)
	x(0)=g2p(0)/b2p(0)
	for i=1,m-1 do begin
		xp=(g2p(i)-(a*x(i-1)))/b2p(i)
		x(i)=xp>f(i)
	endfor
	return,x
endelse
end

function call_diff_finies,s,strike,t,r,sigma,l,n,m,theta,alpha
; l : parametre de localisation du log du support
; n nombre de pas temporels jusqu a t
; m : nombre de pas pour le log du support entre -l et + l
;theta ; parametre du theta schema theta =1 pour un schema totalement implicite , theta=0 pour un schema explicite, theta=0.5 pour krank-nicholson
; alpha parametre de ..
;calcule la valeur du put americain par les differences finies

k=t/float(n)
h=2*l/float(m)
St = alog(s)+(findgen(m)/float(m) - 0.5)*2*l  ; support final
f=fltarr(m)
for i=0,m-1 do f(i)=(strike-exp(St(i)))>0.   ; vecteur de l option a l echeance
a=theta*k*(-sigma^2/(2*h^2)+(r-sigma^2/2.)/(2*h))
b=1+theta*k*(sigma^2/h^2+r)
c=-theta*k*(sigma^2/(2*h^2)+(r-sigma^2/2.)/(2*h))
x=f
g=fltarr(m)
for i=n-2,0,-1 do begin
print,'pas ',i,x
	for j=1,m-2 do g(i)=x(i)+k*(1-theta)*((sigma^2/(2*h^2))*(x(j+1)-2*x(j)+x(j-1))+(alpha/(2*h))*(x(j+1)-x(j-1))-r*x(j))
	x=solve_matrice_sup(a,b,c,f,g)
endfor
if m mod 2 eq 0 then res=(x(m/2)+x(m/2-1))/2. else res=x((m+1)/2)
return,res
end

pro test
s=100.
strike=100.
t=1.
r=0.1
sigma=0.2
l=2
n=11
m=11
theta=1
alpha=0.
print,'resultat=',call_diff_finies(s,strike,t,r,sigma,l,n,m,theta,alpha)
end
