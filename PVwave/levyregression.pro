 
function uauxillevy1,phi,theta,alpha 
   	pi=double(3.1415926535)
	denom=cos(pi*phi/2)
	s=sin(pi*alpha/2*(phi+theta))
 	if (s ne 0.) then return,(s/denom)^(alpha/(1-alpha))* $
		cos(pi/2*((alpha-1)*phi+alpha*theta))/denom else return,0.		
end

function s1,x,alpha,nb
	pi=double(3.1415926535)
	if (x ne 0.) then begin 
 	du=dblarr(nb)
	k=1.-abs(1.-alpha)	
	if x ge 0 then thetastar=k/alpha else thetastar=-k/alpha
 	for i=0,nb-1 do begin
	phi=(1+thetastar)*double(i+0.5)/nb-thetastar
	ux=uauxillevy1(phi,thetastar,alpha)
 	du(i)=ux*exp(-abs(x)^(alpha/(alpha-1))*ux)
	endfor
 ;plot,du*(1+thetastar)*alpha*abs(x)^(1/(alpha-1))/(2*abs(1-alpha)*nb)
  	return,(total(du)*(1+thetastar)*alpha*abs(x)^(1/(alpha-1))/(2*abs(1-alpha)*nb))
	endif else return,gamma(1+1/alpha)/pi*cos(pi/2*(1-abs(1-alpha))/alpha)
end

function sd1,x,alpha,nb
;distribution dans le cas beta =-1 et alpha > 1
	pi=double(3.1415926535)
 	du=dblarr(nb)
	if (x ge 0.) then  sgnx=1. else sgnx=-1.
	k=(1-abs(1-alpha))
  	for i=0,nb-1 do begin
	phi=((1+sgnx*k/alpha)*double(i+0.5)/nb-k*sgnx/alpha)*pi/2
 	valpha=(sin(alpha*phi+pi/2*sgnx*k)/abs(x)/cos(phi))^(alpha/(1-alpha))* $
				cos(phi*(alpha-1)+pi/2*k)/cos(phi)
 	du(i)=exp(-valpha)
	endfor
  	return,(1+sgnx)/2.-sgnx*(total(du)*(1+sgnx*k/alpha))/2/nb 
end

function sd22,x,alpha,beta,nb
;distribution dans le cas beta quelconque x>0 mais alpha > 1
 	pi=double(3.1415926535)
 	du=dblarr(nb)
 	k=(1-abs(1-alpha))
  	for i=1,nb-1 do begin
	phi=((1+k*beta/alpha)*double(i+0.5)/nb-beta*k/alpha)*pi/2.
 	valpha=(sin(alpha*phi+beta*pi/2*k)/abs(x)/cos(phi))^(alpha/(1-alpha))* $
				cos(phi*(alpha-1)+pi/2*beta*k)/cos(phi)
  	du(i)=exp(-valpha)
	endfor
	dphidi=(1+k*beta/alpha)/2/nb*pi
  	return,(1.-total(du)*dphidi/pi)
end



function sd2,x,alpha,beta,nb
;distribution dans le cas beta quelconque mais alpha > 1
if x lt 0 then return,1.-sd22(-x,alpha,-beta,nb)  
if x gt 0 then  return,sd22(x,alpha,beta,nb)
if x eq 0 then return,0.5*(1.-beta*(1-abs(1-alpha))/alpha)
end


function sdg,z
	return,gaussint(z/sqrt(2.))
end

;calcul de l integrale I1 intervenant dans la formule de d evaluation de Mc Culloch
function ki1,alpha,c,z1,c1,c2,b1,b2,nb
	dv=dblarr(nb)
	for i=0,nb-1 do begin
	u=double(i)/nb
	z=u/(1-u)
	h1=-c2*z+alog(s1(z,alpha,nb))+alog((1-sd1(z*b1-b2,alpha,nb))>exp(double(-710)))
	if(abs(h1) lt 710.) then dv(i)=exp(h1)
 	h2=c2*z+alog(s1(-z,alpha,nb))+ alog((1-sd1(-z*b1-b2,alpha,nb))> $
						exp(double(-710)))-2*alog(1-u) 
	if(abs(h2) lt 710.) then dv(i)=dv(i)+exp(h2)
 	endfor
	;plot,dv
	return,total(dv)/nb>0

end


;calcul de l integrale I2 intervenant dans la formule de d evaluation de Mc Culloch
function ki2,alpha,c,z1,c1,c2,b1,b2,nb
 	som=double(0.)
	for i=0,nb-1 do begin
	u=double(i)/nb
	z=u/(1-u)
	h1=-c1*z+alog(s1(z,alpha,nb))+alog(sd1(z*b1+b2,alpha,nb))
	if(abs(h1) lt 710.) then som=som+exp(h1)
	h2=c1*z+alog(s1(-z,alpha,nb))+ alog(sd1(-z*b1+b2,alpha,nb))-2*alog(1-u) 
	if(abs(h2) lt 710.) then som=som+exp(h2)
	endfor
	return,som/nb>0

end

;evaluation du call selon la formule de Mc Culloch
; nb est un nombre de pas d integration dans les calculs
; usuellement pris entre 50 et 100
function LVcall,s,k,t,r,alpha,beta,c,nb
	r1=r
	F=s*exp(r*t)
	pi=double(3.1415926535)
 	alpha1=1/alpha
	c1=((1+beta)/2)^alpha1*c
	c2=((1-beta)/2)^alpha1*c
	theta=pi*alpha/2
	z=1/Cos(theta)
	z1=(alog(F/k)-beta*c^alpha*z)
	b1=c2/c1
	b2=z1/c1
	i1=ki1(alpha,c,z1,c1,c2,b1,b2,nb)

	b1=c1/c2
	b2=z1/c2
 	i2=ki2(alpha,c,z1,c1,c2,b1,b2,nb)
 	return,F*exp(-r1*t+c2^alpha*z)*i1-k*exp(-r1*t+c1^alpha*z)*i2
end

function Cetoile,x0F,alpha,beta,c,nb
	return,LVcall(x0F,1.,1.,0.,alpha,beta,c,0,1,nb)
end

function k_tabulation,alpha,n
case 1 of
	alpha gt 1.9: return,10
	alpha le 1.9 and alpha gt 1.5 : begin
	case 1 of
	n lt 200:return,10
	n ge 200 and n lt 800 :return,closestint(interabcd(1.5,1.9,200,800,11,11,9,9,alpha,n))
	n ge 800 and n lt 1600 :return,closestint(interabcd(1.5,1.9,800,1600,11,11,9,10,alpha,n))
	n ge 1600:return,11
	endcase
	end
	alpha le 1.5 and alpha gt 1.3 : begin
	case 1 of
		n lt 200:return,17
		n ge 200 and n lt 800 :return,closestint(interabcd(1.3,1.5,200,800,22,16,11,11,alpha,n))
		n ge 800 and n lt 1600 :return,closestint(interabcd(1.3,1.5,800,1600,16,14,11,11,alpha,n))
		n ge 1600:return,13
	endcase
	end
	alpha le 1.3 and alpha gt 1.1 : begin
	case 1 of
		n lt 200:return,23
		n ge 200 and n lt 800 :return,closestint(interabcd(1.1,1.3,200,800,24,18,22,16,alpha,n))
		n ge 800 and n lt 1600 :return,closestint(interabcd(1.1,1.3,800,1600,18,15,16,14,alpha,n))
		n ge 1600:return,15
	endcase
	end
	alpha lt 1.1: return,23
endcase
end

function l_tabulation,alpha,n
case 1 of 
	alpha gt 1.9 : return,11
	alpha le 1.9 and alpha gt 1.5 : begin
	case 1 of
		n lt 200:return,11
		n ge 200 and n lt 800 :return,closestint(interabcd(1.5,1.9,200,800,12,14,9,10,alpha,n))
		n ge 800 and n lt 1600 :return,closestint(interabcd(1.5,1.9,800,1600,14,15,10,11,alpha,n))
		n ge 1600:return,13
	endcase
	end
	alpha le 1.5 and alpha gt 1.1 : begin
	case 1 of
		n lt 200:return,14
		n ge 200 and n lt 800 :return,closestint(interabcd(1.1,1.5,200,800,16,18,12,14,alpha,n))
		n ge 800 and n lt 1600 :return,closestint(interabcd(1.1,1.5,800,1600,18,17,14,15,alpha,n))
		n ge 1600:return,16
	endcase
	end
	alpha le 1.1 :return,16
 endcase
end
 
function regresstable,sq
;suppose que sq est un vecteur  stable
;retourne un vecteur [alpha,beta,c,delta] de coefficients stables associes aux rendements
n=n_elements(sq)
s=fltarr(20,n)
zdelta=fltarr(20)
c=fltarr(20)
alpha=fltarr(20)
cestim=fltarr(20)
beta=fltarr(20)
zdeltaestim=fltarr(20)
c(0)=(sq(int(n*0.72))-sq(int(n*0.28)))/1.654
cestim(0)=c(0)
zdelta(0)=avg(sq)
zdeltaestim(0)=zdelta(0)
s(0,*)=(sq-zdelta(0))/c(0)
ck=11
cl=14
for iloop=0,6  do begin
tk=!pi*(findgen(ck)+1)/25
saux=transpose(s(iloop,*))
realphiw1=cos(saux#tk)
imphiw1=sin(saux#tk)
realphi1=fltarr(ck)
imphi1=fltarr(ck)
for i=0,ck-1 do realphi1(i)=avg(realphiw1(*,i))
for i=0,ck-1 do imphi1(i)=avg(imphiw1(*,i))
modphi1=realphi1^2+imphi1^2
window,20,title='regression->alpha,c'
plot,alog(tk),alog(-alog(modphi1)),psym=4
reg=regression(alog(tk),alog(-alog(modphi1)))
alpha(iloop+1)=reg(0)
ck=k_tabulation(alpha(iloop+1),n)
cl=l_tabulation(alpha(iloop+1),n)
c(iloop+1)=exp((reg(1)-alog(2))/alpha(iloop+1))
sprim0=transpose(s(iloop,*)/c(iloop+1))
cestim(iloop+1)=cestim(iloop)*c(iloop+1)
ul=!pi*(findgen(cl)+1)/50
realphiw2=cos(sprim0#ul)
imphiw2=sin(sprim0#ul)
realphi2=fltarr(cl)
imphi2=fltarr(cl)
for i=0,cl-1 do realphi2(i)=avg(realphiw2(*,i))
for i=0,cl-1 do imphi2(i)=avg(imphiw2(*,i))
argphi2=atan(imphi2/realphi2)
vl=tan(!pi*alpha(iloop+1)/2)*abs(ul)^alpha(iloop+1)
for i=0,cl-1 do vl(i)=vl(i)*sgn(ul(i))
reg2=regression2(argphi2,ul,vl)
window,21,title='regression->beta,delta'
plot,argphi2,psym=2
oplot,zdelta(iloop+1)*ul+beta(iloop+1)*vl,color=112,psym=0
beta(iloop+1)=-reg2(1)/c(iloop+1)^alpha(iloop+1)
zdelta(iloop+1)=reg2(0)
zdeltaestim(iloop+1)=zdeltaestim(iloop)+zdelta(iloop+1)*cestim(iloop+1)
print,'pas '+string(iloop+1)+'alpha='+string(alpha(iloop+1))+' cestim='+string(cestim(iloop+1)*n^(1/alpha(iloop+1)))
print,'            beta='+string(beta(iloop+1))+'  zdeltaestim='+string(zdeltaestim(iloop+1))
s(iloop+1,*)=(sq-zdeltaestim(iloop+1))/cestim(iloop+1)
endfor
return,[alpha(iloop),beta(iloop),cestim(iloop),zdeltaestim(iloop)]
end



;procedure chambers ,mallows ,...
function genestable,nbjours,alpha,beta,c,ik
;retourne une liste de prix dont les rendements sont stable
common rando,table,rflag
if rflag eq 0 then begin
	randa=randomu(seed,nbjours)
	wait=randomu(seed,nbjours)
	randb=randomu(seed,nbjours)
endif else begin
	randa=table(2*ik*nbjours:(2*ik+1)*nbjours-1)
	randb=table((2*ik+1)*nbjours*2:(2*ik+2)*nbjours-1)
endelse
kalpha=1-abs(1-alpha)
phi0=-0.5*!pi*beta*kalpha/alpha
phi=double(!pi)*(randa-0.4999)
seed=987345
E=-alog(1.-randb)
q3=fltarr(nbjours)
q4=fltarr(nbjours)
rendement=fltarr(nbjours)
for i=0,nbjours-1 do begin
q3(i)=sin(alpha*(phi(i)-phi0)) / cos( phi(i) - alpha*(phi(i)-phi0) )
q4(i)=Cos( phi(i)-alpha*(phi(i)-phi0) ) / cos(phi(i))
if q4(i) le 0 then begin
print,'*** phi='+string(phi(i))+' alpha='+string(alpha)+' beta='+string(beta)+$
 ' c='+string(c)+' q4='+string(q4(i))+'phi_norm='+string(phi(i)/double(!pi))
q4(i)=0
endif
rendement(i)=q3(i)*q4(i)^(1/alpha)*E(i)^((alpha-1)/alpha)
endfor
s=1
f=fltarr(nbjours)
for i=0,nbjours-1 do begin
s=s*exp(c*rendement(i))
f(i)=s
endfor
return,f
end

 
function teststableKS,s,alpha,beta,cc,delta
;effectue un test de kolmogorov-smirnov sur la structure stable donnee en arguments
;cc est donne par jour
nb=n_elements(s)+1
c=calculdistrib(s)
ncsize=size(c)
nc=ncsize(1)
distribstable=fltarr(nc)
for i=0,nc-1 do distribstable(i)=sd2((c(i,0)-delta)/cc,alpha,beta,100)
window,30,title='distribution empirique/theorique pour le cas stable'
plot,c(*,0),c(*,1)
oplot,c(*,0),distribstable,color=112
distancestable=abs(c(*,1)-distribstable)
Dns=max(distancestable)
nps=n_elements(distancestable)
seuil1=0.01
seuil5=0.05
dn1s=1.6276/sqrt(nps)
dn5s=1.3581/sqrt(nps)
KS1s=Dns-dn1s
KS5s=Dns-dn5s
return,[KS1s,KS5s]
end


function regresscalculistar,y,alpha,beta,gamma,delta
n=n_elements(y)
poids=dblarr(15)
abcisses=dblarr(15)
poids(0:14)=[5.544336631e-2,1.240277389e-1,1.752909438e-1,1.914883407e-1,1.634737971e-1,$
	1.059376372e-1,5.002702115e-2,1.644296900e-2,3.573204214e-3,4.828965093e-4,$
	3.749086502e-5,1.493684115e-6,2.552704969e-8,1.342176791e-10,9.562274467e-14]
abcisses(0:14)=[2.168694746e-2,1.126842203e-1,2.704926714e-1,4.869023703e-1,7.530436830e-1,$
	1.060931003,1.404254958,1.778646379,2.181708131,2.613060845,$
	3.074618113,3.571408151,4.113736089,4.723513062,5.460488935]
istar=double(0)
for j=0,14 do begin
u=abcisses(j)
if alpha ne 1 then omega=double(tan(0.5*!pi*alpha)) else omega=double(2/!pi*alog(abs(u)))
g=double(gamma*abs(u)^alpha*beta*sgn(u)*omega)
c=double(0)
for i=0,n-1 do c=c+double(cos(u*y(i)))
c=double(c/n)
s=double(0)
for i=0,n-1 do s=s+double(sin(u*y(i)))
s=double(s/n)
lambda=double((c-exp(-gamma*abs(u)^alpha)*cos(delta*u+g))^2+(s-exp(-gamma*abs(u)^alpha)*sin(delta*u+g))^2)
istar=istar+poids(j)*lambda
endfor
print,alpha,beta,gamma,delta,istar
return,istar
end

function minimumparabolique,a,b,c,fa,fb,fc
;xx=double(((b-a)^2*(fb-fc)-(b-c)^2*(fb-fa))/((b-a)*(fb-fc)-(b-c)*(fb-fa)))
xx=double(((fa-fc)*(b^2-c^2)+(fb-fc)*(c^2-a^2))/((fa-fc)*(b-c)+(fb-fc)*(c-a)))
return,xx/2.
end


function regresstabstep,x,alpha,beta,gamma,delta,nvar,diff
case 1 of 
(nvar eq 1): begin
	a=alpha+diff & b=alpha & c=alpha-diff
	fa=regresscalculistar(x,alpha+diff,beta,gamma,delta)
	fb=regresscalculistar(x,alpha,beta,gamma,delta)
	fc=regresscalculistar(x,alpha-diff,beta,gamma,delta)
	return,minimumparabolique(a,b,c,fa,fb,fc)
	end
(nvar eq 2): begin
	a=beta+diff & b=beta & c=beta-diff
	fa=regresscalculistar(x,alpha,beta+diff,gamma,delta)
	fb=regresscalculistar(x,alpha,beta,gamma,delta)
	fc=regresscalculistar(x,alpha,beta-diff,gamma,delta)
	return,minimumparabolique(a,b,c,fa,fb,fc)
	end
(nvar eq 3): begin
	a=gamma+diff & b=gamma & c=gamma-diff
	fa=regresscalculistar(x,alpha,beta,gamma+diff,delta)
	fb=regresscalculistar(x,alpha,beta,gamma,delta)
	fc=regresscalculistar(x,alpha,beta,gamma-diff,delta)
	return,minimumparabolique(a,b,c,fa,fb,fc)
	end
(nvar eq 4): begin
	a=delta+diff & b=delta & c=delta-diff
	fa=regresscalculistar(x,alpha,beta,gamma,delta+diff)
	fb=regresscalculistar(x,alpha,beta,gamma,delta)
	fc=regresscalculistar(x,alpha,beta,gamma,delta-diff)
	return,minimumparabolique(a,b,c,fa,fb,fc)
	end
else:return,0
endcase
end

function regresstable1,x
;suppose que x est un vecteur de rendements croissants
n=n_elements(x)
alphabest=1.7
betabest=0.3
gammabest=((x(int(n*0.72))-x(int(n*0.28)))/1.654)^alphabest
deltabest=avg(x)
alphabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,1,0.2)
gammabest=((x(int(n*0.72))-x(int(n*0.28)))/1.654)^alphabest
print,[alphabest,betabest,gammabest,deltabest]
alphabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,1,0.1)
gammabest=((x(int(n*0.72))-x(int(n*0.28)))/1.654)^alphabest
print,[alphabest,betabest,gammabest,deltabest]
gammabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,3,gammabest/2.)
print,[alphabest,betabest,gammabest,deltabest]
gammabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,3,gammabest/4.)
print,[alphabest,betabest,gammabest,deltabest]
betabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,2,0.9)
print,[alphabest,betabest,gammabest,deltabest]
betabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,2,0.3)
print,[alphabest,betabest,gammabest,deltabest]
deltabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,4,max([abs(deltabest)/2,0.0001]))
print,[alphabest,betabest,gammabest,deltabest]
deltabest=regresstabstep(x,alphabest,betabest,gammabest,deltabest,4,max([abs(deltabest)/4,0.00001]))
print,[alphabest,betabest,gammabest,deltabest]
return,[alphabest,betabest,gammabest^(1/alphabest),deltabest]
end


pro testgenestable1,alpha,beta,c,iflg,rflag1
common rando,rtable,rflag
rnumber=1000
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.table',2*nbj*nbs)
x=genestable(rnumber,alpha,beta,c*(1./rnumber)^(1./alpha),0)
window,8,title='cours stable  /  alpha='+string(alpha)+'  beta='+string(beta)+'  c='+string(c)
plot,x

alpha1=findgen(10)/20.+1.45
beta1=findgen(10)/5.-1.
gamma=(c)^alpha/rnumber
res=fltarr(10,10)
for j=0,9 do for i=0,9 do res(i,j)=regresscalculistar(x,alpha1(j),beta1(i),gamma,0)
surface,res,beta1,alpha1

print,'------ si on utilise les valeurs trouvees par regression: on teste le generateur + la regression'
vec=regresstable1(logprix(x))
ras=teststableKS(logprix(x),vec(0),vec(1),vec(2),vec(3))
KS1s=ras(0)
KS5s=ras(1)
if KS1s le 0 then rep1=' = acceptation' else rep1=' = rejet'
if KS5s le 0 then rep5=' = acceptation' else rep5=' = rejet'
print,'kolmogorov-smirnov stable au niveau 1% :'+string(KS1s)+rep1
print,'kolmogorov-smirnov stable au niveau 5% :'+string(KS5s)+rep5
end


pro testgenestable,alpha,beta,c,iflg,rflag1
common rando,rtable,rflag
rnumber=1000
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.table',2*nbj*nbs)
x=genestable(rnumber,alpha,beta,c*(1./rnumber)^(1./alpha),0)
window,8,title='cours stable  /  alpha='+string(alpha)+'  beta='+string(beta)+'  c='+string(c)
plot,x
print,'-----si on utilise les bonne valeurs : on teste le generateur'
ras=teststableKS(logprix(x),alpha,beta,c/rnumber^(1./alpha),0.)
KS1s=ras(0)
KS5s=ras(1)
if KS1s le 0 then rep1=' = acceptation' else rep1=' = rejet'
if KS5s le 0 then rep5=' = acceptation' else rep5=' = rejet'
print,'kolmogorov-smirnov stable au niveau 1% :'+string(KS1s)+rep1
print,'kolmogorov-smirnov stable au niveau 5% :'+string(KS5s)+rep5
print,'------ si on utilise les valeurs trouvees par regression: on teste le generateur + la regression'
vec=regresstable(logprix(x))
ras=teststableKS(logprix(x),vec(0),vec(1),vec(2),vec(3))
KS1s=ras(0)
KS5s=ras(1)
if KS1s le 0 then rep1=' = acceptation' else rep1=' = rejet'
if KS5s le 0 then rep5=' = acceptation' else rep5=' = rejet'
print,'kolmogorov-smirnov stable au niveau 1% :'+string(KS1s)+rep1
print,'kolmogorov-smirnov stable au niveau 5% :'+string(KS5s)+rep5
print,'---- si on regarde l hypothese gaussienne'
vec2=regresgauss(logprix(x))
mean=vec2(0)
sigma=vec2(1)
ra=testgaussKS(logprix(x),mean,sigma)
KS1=ra(0)
KS5=ra(1)
if KS1 le 0 then rep1=' = acceptation' else rep1=' = rejet'
if KS5 le 0 then rep5=' = acceptation' else rep5=' = rejet'
print,'kolmogorov-smirnov gaussien au niveau 1% :'+string(KS1)+rep1
print,'kolmogorov-smirnov gaussien au niveau 5% :'+string(KS5)+rep5
nb=n_elements(logprix(x))+1
c=calculdistrib(logprix(x))
distribgauss=gaussint((c(*,0)-mean)/sigma)
window,9,title='distribution empirique/theorique pour le cas gaussien'
plot,c(*,0),c(*,1)
oplot,c(*,0),distribgauss,color=112
end


function LVcall_delta,s,k,t,r,alpha,beta,c
	nb=50
	pi=double(3.1415926535)
 	alpha1=1/alpha
	c1=((1+beta)/2)^alpha1*c
	c2=((1-beta)/2)^alpha1*c
	theta=pi*alpha/2
	z=1/Cos(theta)
	z1=(alog(s*exp(r*t)/k)-beta*c^alpha*z)
	b1=c2/c1
	b2=z1/c1
	i1=ki1(alpha,c,z1,c1,c2,b1,b2,nb)
print,'delta pour t='+string(t)
 	return,exp(c2^alpha*z)*i1
end



function GESLVcall,f,df,k0,taux,alpha,beta,c,debut,fin,tauxfrais,deltaflag
nb=fin-debut+1
cf=f(debut:fin)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
if deltaflag eq 1 then begin
 for i=0,nb-1 do deltaf(i)=LVcall_delta(cf(i),k0,bdates(i),taux,alpha,beta,c)
endif else begin
vol=c*sqrt(2)
 for i=0,nb-1 do deltaf(i)=call_delta(cf(i),k0,bdates(i),taux,vol)
endelse
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
end
if cf(nb-1) gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else $
				coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)
return,result
end

function MCLVcall1,s01,k01,t01,taux,alpha,beta,c,nbj,nbs,iflg,rflag1
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.table',2*nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
cj=c*(t0/nbj)^(1/alpha)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genestable(nbj+1,alpha,beta,cj,j)
df=findgen(nbj+1)
cout=GEScall(f,df,k0,rj,volj,0,nbj,0.)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))
endif
endfor
mean=avg(vec)
display_distrib,vec,nbs/10,titre,' ',iflg
bs=call(s0,k0,t0,taux,c*sqrt(2.))
print,'alors que black and scholes donne :'+string(bs)
return,mean
end


;calcul de la valeur d un call  dans le cas stable par la methode de Monte Carlo
;le calcul avec le un delta stable de Mc Culloch -> deltaflag=1
;le calcul avec un delta de black and sholes  -> deltaflag=0
;impression = oui -->> iflag=1
;utilisation de table aleatoire predefinit  --> rlfag1 = 1
function MCLVcall,s01,k01,t01,taux,alpha,beta,c,nbj,nbs,deltaflag,iflg,rflag1
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
cj=c*(t0/nbj)^(1/alpha)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genestable(nbj+1,alpha,beta,c,j)
df=findgen(nbj+1)
cout=GESLVcall(f,df,k0,r,alpha,beta,c,0,nbj-1,0.,deltaflag)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))
endif
endfor
mean=avg(vec)
display_distrib,vec,nbs/10,titre,' ',iflg
bs=call(s0,k0,t0,taux,c*sqrt(2.))
print,'alors que black and scholes donne :'+string(bs)
return,mean
end

;teste les differentes procedures d'evaluation d'options en comparant les resultats avec le
;cas gaussien : black and sholes
pro testLVcall,s,k,t,r,alpha,beta,c
 print,'calcul de black and sholes'
 sigm=c*sqrt(2.*365.)
 bs=call(s,k,t,r,sigm)
print,'prix BS = '+string(bs)
print,'calcul de Mc Culloch'
nb=50
mcc=LVcall(s,k,t,r,alpha,beta,c,nb)
print,'prix McC = '+string(mcc)
print,'calcul de monte carlo / loi stable'
nbj=500.
taux=(1+r)^(1/365.)-1
nbs=50
mclv=MCLVcall(s,k,taux,alpha,beta,c,nbj,nbs,0,0,0)
print,'prix MCLV = '+string(mclv)
end


function ni1,X0,F,c
	pi=3.1415926535
	d1=alog(F/X0)/(c*sqrt(2))+c/sqrt(2)
 	i1=Gaussint(d1) 
  	return,i1
	end

function ni2,X0,F,c
	pi=3.1415926535
	d2=alog(F/X0)/(c*sqrt(2))-c/sqrt(2)
 	i2=Gaussint(d2) 
  	return,i2
	end

function callBS,X0,F,c,r1,T
	pi=double(3.1415926535)
   	i1=ni1(X0,F,c)
	print,"x0=",x0,"  F=",F,"  c=",c
	print,"---->i1=",i1
	i2=ni2(X0,F,c)
	print,"x0=",x0,"  F=",F,"  c=",c
	print,"---->i2=",i2
	return,F*exp(-r1*T)*i1-X0*exp(-r1*T)*i2
end 

;fonction qui reproduit 1-gaussint(y) pour y>0
function ni,y,nb
  	som=double(0.)
	pi=double(3.1415926535)
	m=double(y)/(1+y)
	for i=0,nb-1 do begin
	u=m+double(i)*(1-m)/nb
	z=u/(1-u)
	h=-z^2/2-2*alog(1-u)
	if(abs(h) lt 87) then som=som+ exp(h)/sqrt(2*pi)
	endfor
	return,som*(1-m)/nb
end


;fonction de simulation d un processus stable tire de l acticle
function rstab,alpha,bprime,u,w
da=double(0.)
db=double(0.)
piby2=double(1.57079633)
piby4=double(0.785398163)
thr1=double(0.99)
eps=1.0-alpha
phiby2=piby2*(u-0.5)
a=phiby2*tan(phiby2)
bb=tan(eps*phiby2)
b=eps*phiby2*bb
if eps gt -0.99 then tau=bprime/(tan(eps*piby2)*piby2)
if eps le -0.99 then tau=bprime*piby2*eps*(1.-eps)*tan((1.-eps)*piby2)
if a le thr1 then begin
	a2=a^2
	a2p=1.+a2
	a2=1.-a2
	b2=b^2
	b2p=1.+b2
	b2=1.-b2
	endif else begin
da=double(a)^2
db=double(b)^2
a2=1.-da
a2p=1.+da
b2=1.-db
b2p=1.+db
endelse
z=a2p*(b2+2.*phiby2*bb*tau)/(w*a2*b2p)
alogz=alog(z)
d=d2(eps*alogz/(1.-eps))*(alogz/(1.-eps))
rstab=(1.+eps*d)*2.*((a-b)*(1.+a*b)-phiby2*tau*bb*(b*a2-2.*a))/(a2*b2p)+tau*d
return,rstab
end
