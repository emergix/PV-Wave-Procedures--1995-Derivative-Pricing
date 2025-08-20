function multinormale,ya,yb,sigma,n1
;y=[y1...,.ym] et yap <= ybp et y>0
;sigma la matrice de covariance
;n nombre d'essai montecarlo
n=long(n1)
m=n_elements(ya)
up=fltarr(m)
upa=fltarr(m)
ux=fltarr(m)
yap=fltarr(m)
ybp=fltarr(m)
nh=0
for i=0,m-1 do begin
	yap(i)=ya(i)/sqrt(1+ya(i)^2)
	ybp(i)=yb(i)/sqrt(1+yb(i)^2)
endfor
invsigma=invert(sigma)
k=(2*!pi)^(m/2.)*sqrt(determ(sigma))
kp0=1.
for i=0,m-1 do kp0=kp0/(1-yap(i)^2)^1.5
for i=0,m-1 do upa(i)=yap(i)/sqrt(1-yap(i)^2)
kp0=kp0*exp(-0.5*total(upa*(invsigma#upa)))
ur=randomu(seed,(m+1)*n)
for ii=0,n-1 do begin
	for i=0,m-1 do ux(i)=ur(ii*(m+1)+i)*(ybp(i)-yap(i))+yap(i)
	uf=kp0/k*ur(ii*(m+1)+m)
	for i=0,m-1 do up(i)=ux(i)/sqrt(1-ux(i)^2)
	f=exp(-0.5*total(up*(invsigma#up)))/k
	for i=0,m-1 do f=f/(1-ux(i)^2)^1.5
	if uf le f then nh=nh+1
endfor
pp=1.
for i=0,m-1 do pp=pp*(ybp(i)-yap(i))
in=kp0/k*pp*nh/float(n)
return,in
end


pro testq2
ya=[0.5,0.5]
yb=[1.,1.]
sigma=[[0.45^2,0.2*0.45^2],[0.2*0.45^2 ,0.45^2]]
ff=norm2d(1./0.45,1./0.45,0.2)-(norm2d(0.5/0.45,1./0.45,0.2) +$
norm2d(1./0.45,0.5/0.45,0.2))+norm2d(0.5/0.45,0.5/0.45,0.2)
print,multinormale(ya,yb,sigma,3000),' vraie valeur ;',ff
end

  
pro testq3
for a=0.1,0.8,0.1 do begin
	ya=[a,a]
	yb=[1.,1.]
	sigma=[[0.45^2,0.2*0.45^2],[0.2*0.45^2 ,0.45^2]]
	ff=norm2d(1./0.45,1./0.45,0.2)-(norm2d(a/0.45,1./0.45,0.2) +$
	norm2d(1./0.45,a/0.45,0.2))+norm2d(a/0.45,a/0.45,0.2)
	print,multinormale(ya,yb,sigma,3000),' vraie valeur ;',ff
endfor
end


pro testq4
ya=[0.5, 0.3,0.6,0.3,0.6]
yb=[0.7, 0.6,0.9,0.6,1.0]
sigma=[	[0.3^2,0,0,0,0],$
	[0,0.4^2,0,0,0],$
	[0,0,0.5^2,0,0],$
	[0,0,0,0.4^2,0],$
	[0,0,0,0,0.3^2]]
ff1=gaussint(0.7/0.3)-gaussint(0.5/0.3)
ff2=gaussint(0.6/0.4)-gaussint(0.3/0.4)
ff3=gaussint(0.9/0.5)-gaussint(0.6/0.5)
ff4=gaussint(0.6/0.4)-gaussint(0.3/0.4)
ff5=gaussint(1.0/0.3)-gaussint(0.6/0.3)
ff=ff1*ff2*ff3*ff4*ff5
print,multinormale(ya,yb,sigma,10000),' vraie valeur ;',ff
end

pro testq5
ya=[0.5, 0.3,0.6,0.3,0.6,0.5, 0.3,0.6,0.3,0.6]
yb=[0.7, 0.6,0.9,0.6,1.0,0.7, 0.6,0.9,0.6,1.0]
sigma=[	[0.3^2,0,0,0,0,0,0,0,0,0],$
	[0,0.4^2,0,0,0,0,0,0,0,0],$
	[0,0,0.5^2,0,0,0,0,0,0,0],$
	[0,0,0,0.4^2,0,0,0,0,0,0],$
	[0,0,0,0,0.3^2,0,0,0,0,0],$
	[0,0,0,0,0,0.3^2,0,0,0,0],$
	[0,0,0,0,0,0,0.4^2,0,0,0],$
	[0,0,0,0,0,0,0,0.5^2,0,0],$
	[0,0,0,0,0,0,0,0,0.4^2,0],$
	[0,0,0,0,0,0,0,0,0,0.3^2]]
ff1=gaussint(0.7/0.3)-gaussint(0.5/0.3)
ff2=gaussint(0.6/0.4)-gaussint(0.3/0.4)
ff3=gaussint(0.9/0.5)-gaussint(0.6/0.5)
ff4=gaussint(0.6/0.4)-gaussint(0.3/0.4)
ff5=gaussint(1.0/0.3)-gaussint(0.6/0.3)
ff6=gaussint(0.7/0.3)-gaussint(0.5/0.3)
ff7=gaussint(0.6/0.4)-gaussint(0.3/0.4)
ff8=gaussint(0.9/0.5)-gaussint(0.6/0.5)
ff9=gaussint(0.6/0.4)-gaussint(0.3/0.4)
ff10=gaussint(1.0/0.3)-gaussint(0.6/0.3)
ff=ff1*ff2*ff3*ff4*ff5*ff6*ff7*ff8*ff9*ff10
print,multinormale(ya,yb,sigma,20000),' vraie valeur ;',ff
end





pro tq
n=long(50)
nh=0
nh2=0
sh=0.
sh2=0.
kn=6
a=0.3
b=0.5
sigma=0.3
c0=1/(sqrt(2*!pi)*sigma)
print,'c0=',c0
c=c0
cste=0.7
c2=2*cste
ur=randomu(seed,2.*n)
ur1=fltarr(n)
ur2=fltarr(n)
ur22=fltarr(n,2)
for i=long(0),n-1 do begin
ur1(i)=ur(2*i)
ur2(i)=ur(2*i+1)
ur22(i,0)=ur(2*i)
ur22(i,1)=ur(2*i+1)
endfor
print,'avg=',total(ur)/(2*float(n))
print,'cov=',covariances(ur22)
ws=sort(ur1)
px=ur1(ws)
py=ur2(ws)
pxx=fltarr(n)
pyy=fltarr(n)
pyf=fltarr(n)
for i=long(0),n-1 do begin
	pxx(i)=px(i)*(b-a)+a
	pyy(i)=py(i)*c
	pyf(i)=(exp(-pxx(i)^2./(2.*sigma^2.)))/(sqrt(2.*!pi)*sigma)
endfor
wy=where(pyy gt pyf)
nfu=n_elements(wy)
print,'nfu=',nfu

Plot,pxx,pyy,psym=3
oplot,pxx,pyf,psym=-1,thick=0.2

for i=long(0),n-1 do begin
	u=[ur(2*i),ur(2*i+1)]
	xu=(u(0)*(b-a))+a
	fu=(1/(sqrt(2*!pi)*sigma))*exp(double(-xu^2/(2*sigma^2)))
	fu2=(xu-a)^kn
	yu=u(1)*c
	yu2=u(1)*c2
	if fu ge yu then nh=nh+1
	if fu2 ge yu2 then nh2=nh2+1
	sh=sh+fu
	sh2=sh2+fu2
endfor
i=(b-a)*c*nh/float(n)
i2=(b-a)*c2*nh2/float(n)
is=(b-a)/float(n)*sh
is2=(b-a)/float(n)*sh2
print,'result MC=',i,is,'  vraie valeur=',gaussint(b/sigma)-gaussint(a/sigma),nh
print,'result MC=',i2,is2,'  vraie valeur=',((b-a)^(kn+1.))/(kn+1.),nh2
end
