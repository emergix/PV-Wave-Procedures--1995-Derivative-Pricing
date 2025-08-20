function f2,x,y
return,abs((x*x-1.)/(y*y))
end
function f3,x,y
return,abs(-x*(x*x-3.)/(y*y*y))
end
function f4,x,y
return,abs(3.+x*x*(-6.+x*x))/(y^4))
end
function f5,x,y
return,abs(x*(-5.+x*x*(10.-x*x*)))/y^5
end
function f6,x,y
return,abs(-15.+x*x*(45.-x*x*(-15.+x*x)))/y^6
end



function mulnor,a,b,sig,eps,n,inf,prob,bound,ifault
coef=[[0.311111111111,1.4222222222222,0.5333333333333,1.42222222222222,0.31111111111111111],$
      [0.333333333333,0.,1.333333333333,0.,0.33333333333333],$
      [0.5,0.,0.,0.,0.5]]
bcs=[1.,4.,6.,4.,1.]
bcn=[1.,6.,15.,20.,15.,6.,1.]
do=[-3.6097169,-3.324574,-2.85697,-2.7969244,-2.3344142,-1.8891759,-1.7320508,-1.3556262,-1.0,0.74196378,$
-0.61670659,-0.38361163,0.0,0.38361163,0.61670659,0.74196378,1.0,1.3556262,1.7320508,1.8891759,2.3344142,2.7969244,$
2.85697,3.324574,3.6097169]
eo=[-2.85697,-2.3344144,-1.7320508,-1.3556262,-1.0,-0.74196378,0.0,0.74196378,1.0,1.3556262,1.7320508,2.3344142,2.85697]
cons2=6.879833
cons3=4.517004
cons4=76.371214
epsmin=1.e-8
epssim=6.e-5
true=1
false=0
;checking faulty data
ifault=0
if eps le epsmin then ifault=2
if n le 0 or n gt 7 then ifault =3
if ifault ne 0 then return,ifault
for i=1,n do if inf(i) eq 2 and a(i) lt b(i) then ifault=100+i
if ifault ne 0 then return,ifault
bound=eps
ept=eps*p15/float(n)
z=-gauin(ept,ifault)+epsmin
if ifault ne 0 then return
cup=alnorm(z,true)
; inverting sig and the n-2 lower right hand principal minors
ik=0
ij=0
for i=1,n do for j=1,i do begin
	ik=ik+1
	if i eq j then goto,ll2
	ij=ij+1
	sigma(ik)=sig(ij)
	goto,ll3
	ll2:sigma(ik)=1.
ll3:endfor
if n le 2. then goto,ll4
det=determ(sigma)
sinv=inv(sigma)
simps=true
if det lt p05 or eps le epssim then simps=false
prob=0
