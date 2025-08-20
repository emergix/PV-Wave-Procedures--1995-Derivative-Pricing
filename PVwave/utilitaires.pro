
function sgn,x
if x gt 0 then return,1.
if x lt 0 then return,-1.
if x eq 0 then return,0.
end

function int,x
return,x-(x mod double(1.))
end

function ceilling,x
a=int(x)
if x eq a then return,x else return,a+1
end



function closestint,x
z=x mod 1.
if z ge 0.5 then return,x-z+1 else return,x-z
end

function lecture,z
p='/home/p6379j/'+z+'.histo'
openr,1,p
date=fltarr(2000)
open=date & high=date & low=date
close=date & vol=date & interest=date
i=0
while not eof(1) do begin
readf,1,dates,opens,highs,lows,closes,vols,interests
date(i)=dates & open(i)=opens & high(i)=highs & low(i)=lows
close(i)=closes & vol(i)=vols & interest(i)=interests
i=i+1
endwhile
print,i-1,'record'
close,1
r=fltarr(7,i)
r(0,*)=date(0:i-1)
r(1,*)=open(0:i-1)
r(2,*)=high(0:i-1)
r(3,*)=low(0:i-1)
r(4,*)=close(0:i-1)
r(5,*)=vol(0:i-1)
r(6,*)=interest(0:i-1)
return,r
end

;lecture des table de nombre aleatoire
function lectureR,z,nbmax
p='/vol/mathematica/appli/'+z
openr,1,p
rand=fltarr(200001)
i=double(0)
while not eof(1) do begin
readf,1,nbr
rand(i)=nbr
i=i+1
if i ge nbmax then  begin
	print,i,'record'
	close,1
	r=fltarr(i)
	r(*)=rand(0:i-1)
	return,r
endif
endwhile
print,i-1,'record'
close,1
r=fltarr(i)
r(*)=rand(0:i-1)
return,r
end

function lecture2R,z,nbmax1
p='/vol/mathematica/appli/'+z
openr,1,p
rand=dblarr(1000000)
nbr=double(0.0)
i=double(0)
nbmax=double(nbmax1)
while not( eof(1) or i gt nbmax) do begin
readf,1,nbr
rand(i)=nbr
i=i+1
if i ge nbmax then  begin
	print,i,'record'
	close,1
	r=fltarr(i)
	r(*)=rand(0:i-1)
	return,r
endif
endwhile
print,i-1,'record'
close,1
r=dblarr(i)
r(*)=rand(0:i-1)
return,r
end

function mapstring,str
n=n_elements(str)
a=''
for i=0,n-1 do a=a+string(str(i))
return,'['+a+']'
end


function interunion,a1,a2
a3=fltarr(13,2000)
i1=0 & i2=0 & i3=0
l1=n_elements(a1)/7-1 & l2=n_elements(a2)/7-1
while (i1 le l1) and (i2 le l2) do begin
if a1(0,i1) eq a2(0,i2) then begin
a3(0:6,i3)=a1(0:6,i1) & a3(7:12,i3)=a2(1:6,i2)
i1=i1+1 & i2=i2+1 & i3=i3+1
                         endif else begin
if a1(0,i1) lt a2(0,i2) then i1=i1+1 else i2=i2+1
                                    endelse
   end
a4=fltarr(13,i3)
a4=a3(*,0:i3-1)
return,a4
end

function nbjdate,d
j=d mod 100
m=((d-j) mod 10000)/100
a=(d-j-100*m)/10000+1900
return,julday(m,j,a)
end

function extractdates,f
si=size(f)
n=si(2)
dates=f(0,0:n-1)
df=fltarr(n)
b0=nbjdate(dates(0))
for i=0,n-1 do df(i)=nbjdate(dates(i))
return,df
end



function jj,d
;calcule le jour julien a partir du jour jjmmaa , inverse de caldat
a=d mod 100
m=((d-a) mod 10000)/100
j=(d-a-100*m)/10000
if a lt 10 then a=a+100
return,julday(m,j,a+1900)
end

function jj_dt,d
;calcule le jour julien a partir du jour jjmmaa , inverse de caldat
a=d mod 100
m=((d-a) mod 10000)/100
j=(d-a-100*m)/10000
if a lt 10 then a=a+100
return,var_to_dt(a+1900,m,j,12.,0.,0.)
end

function jj2,d
;calcule le jour julien a partir du jour aammjj , inverse de caldat
j=d mod 100
m=((d-j) mod 10000)/100
a=(d-j-100*m)/10000
if a lt 10 then a=a+100
return,julday(m,j,a+1900)
end

function jj2_dt,d
;calcule le jour julien a partir du jour aammjj , inverse de caldat
j=d mod 100
m=((d-j) mod 10000)/100
a=(d-j-100*m)/10000
if a lt 10 then a=a+100
return,var_to_dt(a+1900,m,j,12.,0.,0.)
end


function jj3,d
;calcule le jour julien a partir du jour mmjjaa , inverse de caldat
a=d mod 100
j=((d-a) mod 10000)/100
m=(d-a-100*j)/10000
if a lt 10 then a=a+100
return,julday(m,j,a+1900)
end

function jj3_dt,d
;calcule le jour julien a partir du jour mmjjaa , inverse de caldat
a=d mod 100
j=((d-a) mod 10000)/100
m=(d-a-100*j)/10000
if a lt 10 then a=a+100
return,var_to_dt(a+1900,m,j,12.,0.,0.)
end


;fonction qui convertit une designation de mois/annee en jours juilien permettant par difference
;les calculs de delai entre dates d est egal pra exemple a 0392 ce qui signifie fin mars 1992
function je_option,d
a=d mod 100
m=(d-a)/100
if m eq 1 then j=31
if m eq 2 then j=28
if m eq 3 then j=31
if m eq 4 then j=30
if m eq 5 then j=31
if m eq 6 then j=30
if m eq 7 then j=31
if m eq 8 then j=31
if m eq 9 then j=30
if m eq 10 then j=31
if m eq 11 then j=30
if m eq 12 then j=31
return,julday(m,j,a+1900)
end

;transforme l entier jjmmaa en aammjj
function inverse_date,z
a=z-int(z/100.)*100
j=int(z/10000.)
mstar=z-a-j*10000
return,mstar+j+a*10000.
end

function caldat,julian
;calcule le jour jjmmaa a partir du jour julien : fonction inverse de jj
igreg=double(2299161.)
if julian ge igreg then begin
jalpha=int(((julian-double(1867216.))-0.25)/36524.25)
ja=julian+1.+jalpha-int(0.25*jalpha)
endif else begin
ja=julian
endelse
jb=ja+1524.
jc=int(6680.+((jb-double(2439870.))-122.1)/365.25)
jd=365.*jc+int(0.25*jc)
je=int((jb-jd)/30.6001)
id=jb-jd-int(30.6001*je)
mm=je-1.
if mm gt 12. then mm=mm-12.
iyyy=jc-4715.
if mm gt 2. then iyyy=iyyy-1.
if iyyy le 0. then iyyy=iyyy-1.
yy=iyyy-1900.
result=id*10000.+mm*100.+yy
return,result
end

function caldat2,julian
;calcule le jour aammjj a partir du jour julien : fonction inverse de jj2
igreg=double(2299161.)
if julian ge igreg then begin
jalpha=int(((julian-double(1867216.))-0.25)/36524.25)
ja=julian+1.+jalpha-int(0.25*jalpha)
endif else begin
ja=julian
endelse
jb=ja+1524.
jc=int(6680.+((jb-double(2439870.))-122.1)/365.25)
jd=365.*jc+int(0.25*jc)
je=int((jb-jd)/30.6001)
id=jb-jd-int(30.6001*je)
mm=je-1.
if mm gt 12. then mm=mm-12.
iyyy=jc-4715.
if mm gt 2. then iyyy=iyyy-1.
if iyyy le 0. then iyyy=iyyy-1.
yy=iyyy-1900.
result=yy*10000.+mm*100.+id
return,result
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

function intervalue,courbe,delai    ;calcule la valeur exacte  ou zero
zsize=size(courbe)
n=zsize(1)
if n eq 1 then return,courbe(0,1)
win=where(courbe(*,0) ge delai,ct)
if ct eq 0 then return,0
if delai lt 0 then return,0
winmin=min(win)
if winmin eq 0 then return,courbe(0,1)
return,courbe(winmin,1)
end


function interpolation1,tvector,svector,t
r=tvector#[1,1]
r(*,1)=svector
return,interpolation(r,t)
end


function intervalue1,courbe,delai1,delai2   	; calcule la valeur estime  ou zero
   						; pour avoir une valeur integree exacte
						; (important pour les dividendes)
						; entre les bornes delai1 et delai2

zsize=size(courbe)
n=zsize(1)
t1=courbe(0,0)
t2=courbe(n-1,0)
sum=0.
delai1a=delai1>0
delai2a=delai2<t2
win1=where(courbe(*,0) ge delai1a,ct1)
n1=min(win1)
win2=where(courbe(*,0) ge delai2a,ct2)
n2=min(win2)
sum=sum+((courbe(n1,0)-delai1a)*courbe(n1,1))
sum=sum-((courbe(n2,0)-delai2a)*courbe(n2,1))
if n1 eq n2 then return,sum
for i=n1,n2-1 do sum=sum+((courbe(i+1,0)-courbe(i,0))*courbe(i+1,1))
return,sum/float(delai2-delai1)
end

function moyenne_courbe,courbe,delai
if delai eq 0 then print,'moyenne_courbe:erreur delai =0 '
zsize=size(courbe)
n=zsize(1)
t1=courbe(0,0)
f1=courbe(0,1)
if delai lt t1 then return,f1
sigma_s=f1*t1
i=1
while i lt n  do  begin
	if courbe(i,0) eq delai then return,sigma_s/delai
	if courbe(i,0) gt delai then begin
		sigma_s=sigma_s+(delai-t1)*(2*f1+(courbe(i,1)-f1)*(delai-t1)/(courbe(i,0)-t1))/2.
		return,sigma_s/delai
		endif
	sigma_s=sigma_s+(courbe(i,0)-t1)*(courbe(i,1)+f1)/2.
	t1=courbe(i,0)
	f1=courbe(i,1)
	i=i+1	
endwhile
sigma_s=sigma_s+(delai-t1)*f1
return,sigma_s/delai
end

function moyenne_courbe2,courbe,delai1,delai2
f1=moyenne_courbe(courbe,delai1)*delai1
f2=moyenne_courbe(courbe,delai2)*delai2
return,(f2-f1)/(delai2-delai1)
end

 ; fonction utilise pour calcule somme de sigma^2 *dt
function integrale_carre_courbe,courbe,delai    ;calcule l integrale jusqu a delai de courbe^2
zsize=size(courbe)
n=zsize(1)
courbe2=fltarr(n,2)
for i=0,n-1 do begin
	courbe2(i,0)=courbe(i,0)
	courbe2(i,1)=courbe(i,1)^2
endfor
return,delai*moyenne_courbe(courbe2,delai)
end


;fonction supose constante entre t(i) et t(i+1) et valant f(i+1)
; la fonction est supose valoir f(n-1) pour x>t(n-1)
; elle est suppose valoir f(0) pour x<t(0)
function integrale_courbec,f,tt1,tt2
if tt1 eq tt2 then return,0
if tt1 gt tt2 then begin 
	t1=tt2
	t2=tt1
	signes=-1
endif else begin
	t1=tt1
	t2=tt2
	signes=1
endelse
zsize=size(f)
n=zsize(1)
t0=f(0,0)
f0=f(0,1)
sigma_s=0
if t1 le t0 and t2 le t0 then return,(t2-t1)*f0*signes
if t1 ge f(n-1,0) then return,f(n-1,1)*(t2-t1)*signes
nt1=min(where(f(*,*) ge t1))
if nt1 lt 0 then return,f(n-1,1)*(t2-t1)*signes
nt2=max(where(f(*,0) le t2))
if nt2 lt 0 then begin
	nt2=n-1
	plusfin=(t2-f(n-1,0))*f(n-1,1)
	t2=f(n-1,0)
endif else plusfin=0.
;intgrale entre t1 et nt1
plus1=f(nt1,1)*(f(nt1,0)-t1)
;integrale entre nt2 et t2
if nt2 ne n-1 then ft2=f(nt2+1,1) else ft2=f(n-1,1)
plus2=ft2*(t2-f(nt2,0))
if nt1 eq nt2 then return,(plus1+plus2+plusfin)*signes
if nt1 gt nt2 then begin
	h=f(nt1,1)*(f(nt1,0)-f(nt2,0))
	return,(plus1+plus2-h)*signes
	endif
;cas standard nt1<nt2
h=0.
for ideb=nt1,nt2-1 do h=h+f(ideb+1,1)*(f(ideb+1,0)-f(ideb,0))
return,(h+plus1+plus2+plusfin)*signes
end


; fonction qui supose continue la fonction
; et constante en dehors de ]t(0),t(n-1)[
function integrale_courbe,f,tt1,tt2
if tt1 eq tt2 then return,0
if tt1 gt tt2 then begin 
	t1=tt2
	t2=tt1
	signes=-1
endif else begin
	t1=tt1
	t2=tt2
	signes=1
endelse
zsize=size(f)
n=zsize(1)
t0=f(0,0)
f0=f(0,1)
sigma_s=0
if t1 le t0 and t2 le t0 then return,(t2-t1)*f0*signes
if t1 ge f(n-1,0) then return,f(n-1,1)*(t2-t1)*signes
nt1=min(where(f(*,*) ge t1))
if nt1 lt 0 then return,f(n-1,1)*(t2-t1)*signes
nt2=max(where(f(*,0) le t2))
if nt2 lt 0 then begin
	nt2=n-1
	plusfin=(t2-f(n-1,0))*f(n-1,1)
	t2=f(n-1,0)
endif else plusfin=0.
if nt1 ne 0 then ft1=(f(nt1,1)-f(nt1-1,1))/(f(nt1,0)-f(nt1-1,0))*(t1-f(nt1-1,0))+f(nt1-1,1) else ft1=f(0,1)
plus1=(f(nt1,1)+ft1)*(f(nt1,0)-t1)/2.
if nt2 ne n-1 then ft2=(f(nt2+1,1)-f(nt2,1))/(f(nt2+1,0)-f(nt2,0))*(t2-f(nt2,0))+f(nt2,1) else ft2=f(n-1,1)
plus2=(f(nt2,1)+ft2)*(t2-f(nt2,0))/2.
if nt1 eq nt2 then return,(plus1+plus2+plusfin)*signes
if nt1 gt nt2 then begin
	h=(f(nt2,1)+f(nt1,1))*(f(nt1,0)-f(nt2,0))/2.
	return,(plus1+plus2-h)*signes
	endif
;cas standard nt1<nt2
h=0.
for ideb=nt1,nt2-1 do h=h+(f(ideb+1,1)+f(ideb,1))*(f(ideb+1,0)-f(ideb,0))/2.
return,(h+plus1+plus2+plusfin)*signes
end

;fonction qui calcule le y tel que somme(x a y de f(t)dt) =k
;en suposant que f est constante entre t(i) et t(i+1) et vaut f(i+1)
; et que f est toujours positive
;on suppose donc que f(x) =f(n-1) ou t(n-1) est le dernier point de la coube
;
function integral_solve_y,f,t1,kk
k=kk
zsize=size(f)
n=zsize(1)
t0=f(0,0)
f0=f(0,1)
sigma_s=0
if k eq 0 then return,t2
if t1 le t0 then sigma_s=(t0-t1)*f0
if k le sigma_s then return,k/f0+t1
if t1 ge f(n-1,0) then return, k/float(f(n-1,1))+t1
nt1=min(where(f(*,0) gt t1))
plus1=f(nt1,1)*(f(nt1,0)-t1)
if k lt plus1 then return,t1+k/f(nt1,1)
k=k-plus1
nt1=nt1+1
if nt1 ge n then return,f(n-1,0)+k/float(f(n-1,1))
plus1=f(nt1,1)*(f(nt1,0)-f(nt1-1,0))
while 1 do begin
if k lt plus1 then return,f(nt1-1,0)+k/f(nt1,1)
k=k-plus1
nt1=nt1+1
if nt1 ge n then return,f(n-1,0)+k/float(f(n-1,1))
plus1=f(nt1,1)*(f(nt1,0)-f(nt1-1,0))
endwhile

end


function shift_courbe,taux,nyears    
			;shifte une courbe en simulant une date future 
ns=size(taux)
n=ns(1)
nt=0
flag=0
for i=0,n-1 do begin
	if taux(i,0) gt nyears then nt=nt+1
	if taux(i,0) eq nyears then flag =1
endfor
if flag eq 0 then begin 
	if nt gt 0 then res=fltarr(nt+1,2) else return,0
	j=0
	for i=0,n-1 do begin
  	  if taux(i,0) gt nyears then begin
		res(j+1,0)=taux(i,0)-nyears
		res(j+1,1)=taux(i,1)
		j=j+1
 	  endif
	endfor
	res(0,0)=0
	res(0,1)=interpolation(taux,nyears)
endif else begin 
	if nt gt 0 then res=fltarr(nt,2) else return,0
	j=0
	for i=0,n-1 do begin
  	  if taux(i,0) gt nyears then begin
		res(j,0)=taux(i,0)-nyears
		res(j,1)=taux(i,1)
		j=j+1
 	  endif
	endfor
endelse
return,res
end
pro tt
a=transpose([[1,2],[2,3],[3,4]])
end

function shift_dt1,dt,nyears   
			;shifte une liste de date en simulant une date future 
n=n_elements(dt)
nt=0
for i=0,n-1 do begin
	if dt(i) gt nyears then nt=nt+1
endfor
if nt gt 0 then res=fltarr(nt) else return,0
j=0
for i=0,n-1 do begin
  if dt(i,0) gt nyears then begin
	res(j)=dt(i)-nyears
	j=j+1
  endif
endfor
return,res
end


function shift_dt2,dataelements,dt,nyears    
			;shifte une liste de nombre en simulant 
			;une date future en s aidant d une liste de date 
n=n_elements(dt)
nt=0
for i=0,n-1 do begin
	if dt(i) gt nyears then nt=nt+1
endfor
if nt gt 0 then res=fltarr(nt) else return,0
j=0
for i=0,n-1 do begin
  if dt(i) gt nyears then begin
	res(j)=dataelements(i)
	j=j+1
  endif
endfor
return,res
end



; regression : y=a*x+b --> a,b
function regression,x,y
sx=total(x)
n=n_elements(x)
a=(sx*total(y)-n*total(x*y))/(sx^2-total(x*x)*n)
b=(total(y)-a*sx)/n
return,[a,b]
end



; regression : z=a*x+b*y --> a,b
function regression2,z,x,y
yz=total(y*z)
y2=total(y^2)
xy=total(x*y)
a=(total(x*z)*y2-yz*xy)/(y2*total(x^2)-xy^2)
b=(yz-a*total(x*y))/y2
return,[a,b]
end

;regression  :z=a*x^2+b*x+c -->a,b,c
function regression4,z,x
xi4=total(x^4)
xi3=total(x^3)
xi2=total(x^2)
xi=total(x)
sn=n_elements(x)
zixi2=total(z*x^2)
zixi=total(z*x)
zi=total(z)
matrice=double([[xi4,xi3,xi2],[xi3,xi2,xi],[xi2,xi,sn]])
coeff=transpose(invert(matrice))#[zixi2,zixi,zi]
return,coeff
end


;regression  :z=a*x^3+b*x^2+c*x+d -->a,b,c,d
function regression5,z,x
xi6=total(x^6)
xi5=total(x^5)
xi4=total(x^4)
xi3=total(x^3)
xi2=total(x^2)
xi=total(x)
sn=n_elements(x)
zixi3=total(z*x^3)
zixi2=total(z*x^2)
zixi=total(z*x)
zi=total(z)
matrice=double([[xi6,xi5,xi4,xi3],[xi5,xi4,xi3,xi2],[xi4,xi3,xi2,xi],[xi3,xi2,xi,sn]])
coeff=transpose(invert(matrice))#[zixi3,zixi2,zixi,zi]
return,coeff
end

;regression : z=a*x+b*y+c --> a,b,c
function regression3,z,x,y
sx=total(x)
sy=total(y)
xy=total(x*y)
x2=total(x^2)
y2=total(y^2)
m=[[sx,sy,1],[x2,xy,sx],[xy,y2,sy]]
vz=[total(z),total(z*x),total(z*y)]
return,transpose(invert(m))#vz
end


; regression : z=a*y^2+b*y+c*x+d --> a,b,c,d
function regression7,z,x,y
sx=double(total(x))
sy=double(total(y))
xy=double(total(x*y))
x2=double(total(x^2))
y2=double(total(y^2))
x3=double(total(x^3))
x2y=double(total((x^2)*y))
x4=double(total(x^4))
m=[	[x2,  sx,  sy, 1],$
	[x2y, xy , y2, sy],$
	[x3,  x2,  xy, sx],$
	[x4,  x3,  x2y,x2]]
vz=[total(z),total(z*y),total(z*x),total(z*x^2)]
return,transpose(invert(m))#vz
end


; regression : z=a*x^2+b*x+c*y^2+d*y+e*xy+f --> a,b,c,d,e,f
function regression6,z,x,y
sx=double(total(x))
sy=double(total(y))
xy=double(total(x*y))
x2=double(total(x^2))
y2=double(total(y^2))
x3=double(total(x^3))
y3=double(total(y^3))
x2y=double(total((x^2)*y))
xy2=double(total((y^2)*x))
x4=double(total(x^4))
y4=double(total(y^4))
x3y=double(total((x^3)*y))
xy3=double(total((y^3)*x))
x2y2=double(total((x^2)*(y^2)))
m=[	[x2,  sx, y2,  sy, xy,  1],$
	[x3y, x2y,xy3, xy2,x2y2,xy],$
	[x2y, xy ,y3,  y2, xy2, sy],$
	[x2y2,xy2,y4,  y3, xy3, y2],$
	[x3,  x2, xy2, xy, x2y, sx],$
	[x4,  x3, x2y2,x2y,x3y, x2]]
vz=[total(z),total(z*x*y),total(z*y),total(z*y^2),total(z*x),total(z*x^2)]
return,transpose(invert(m))#vz
end


function find_min_parab,x,y
;trouve un minimum parabolique
n=n_elements(x)
ymini=min(y)
sindx=where(y eq ymini)
indx=sindx(0)
if indx eq 0 then return,x(0)
if indx eq 1 then begin
	coeff=regression4([y(0),y(1),y(2)],[x(0),x(1),x(2)])
	return,-coeff(1)/(2*coeff(0))
endif
if indx eq n-1 then return,x(n-1)
if indx eq n-2 then begin
	coeff=regression4([y(n-3),y(n-2),y(n-1)],[x(n-3),x(n-2),x(n-1)])
	return,-coeff(1)/(2*coeff(0))
endif
; regression parabolique
if indx gt 2 then x1=x(0:indx-1) else x1=x(0:indx)
if indx gt 2 then y1=y(0:indx-1) else y1=y(0:indx)
coeff1=regression4(y1,x1)
a=coeff1(0)
b=coeff1(1)
return,-b/(2*a)
end


function deltazreg3,z,x,y
reg=regression3(z,x,y)
a=reg(0)
b=reg(1)
c=reg(2)
return,z-(a*x+b*y+c)
end


;calcul d interpolation dans une table a double entree ou on on a les 
;abcisses n1 et n2 et les ordonnee alpha1 , alpha2 , les valeurs de la table sont alphai_nj
; on cherche alors la valeur associe a alpha n que l on retourne 
function interabcd,alpha1,alpha2,n1,n2,alpha1_n1,alpha1_n2,alpha2_n1,alpha2_n2,alpha,n
palpha=(alpha-alpha1)/(alpha2-alpha1)
g1=palpha*(alpha2_n1-alpha1_n1)+alpha1_n1
g2=palpha*(alpha2_n2-alpha1_n2)+alpha1_n2
pn=(n-n1)/(n2-n1)
return,g1+pn*(g2-g1)
end


;calcul de la factorielle
function fact,n
r=double(1.)
for i=1,n do r=r*i
return,r
end


;calcul du coefficient en x^n2 de ((1+x)/2)^n1
function binomtwo,n1,n2
return,(0.5^n1*fact(double(n1))/fact(double(n1-n2)))/fact(double(n2))
end

;fonction qui retourne les p premiers elements (si p >0)
; ou les p dernier elements (si p <0 ) d un vecteur
function cut,a,p
if p ge 0 then return,a(p:n_elements(a)-1) else return,a(0:n_elements(a)+p-1)
end


;calcule les rendements
function logprix,x
n=n_elements(x)-1
ss=fltarr(n)
for i=0,n-1 do ss(i)=alog(x(i+1)/x(i))
sq=ss(sort(ss))
return,sq
end



;fonction qui rettourne la courbe de distribution associe au vecteur de prix x
;le resultat se presente sous la forme d une matrice return/probabilite accumule etl que
; retour(*,0)=returns et retour(*,1)=probabilite accumulee (distributuion)
function calculdistrib,x
nb=n_elements(x)
s=double(x(sort(x)))
z=fltarr(nb,2)
for i=0,nb-1 do begin
z(i,0)=s(i)
z(i,1)=float(i)/float(nb)
endfor
vm=fltarr(nb,2)
m=-1.
im=0
for i=0,nb-1 do begin
if z(i,0) gt m then begin
vm(im,0)=z(i,0)
vm(im,1)=z(i,1)
im=im+1
m=z(i,0)
endif
endfor
z1=fltarr(im,2)
z1(*,0)=vm(0:im-1,0)
z1(*,1)=vm(0:im-1,1)
return,z1
end

;fonction qui echantillonne les fonctions croissante monotone
;retourne un resultat z tel x=z(*,0),y=z(*,1) 
function echantillonne,x,y,n1
n=n1-1
z=fltarr(n+1)
zx=fltarr(n+1,2)
minx=min(x)
maxx=max(x)
ik=0
in=1
iz=intarr(n+1)
iz(0)=0
deltax=float(maxx-minx)/n
deltasamplex=float(maxx-minx)/(n_elements(x)-1)
zx(*,0)=minx+(findgen(n+1)+0.5)*deltax
debutloop:accum=0
nacum=0
debutaccum :if x(ik) le minx+in*deltax then begin
 	accum=accum+y(ik)
	nacum=nacum+1
 	ik=ik+1
	if ik ge n_elements(x) then begin
	;nous somme a la fin des points sans avoir atteint une borne (generalite)
 	iborneinf=iz(in-1)-1
 	ibornesup=ik-1
	idebutsample=ik-nacum
	ifinsample=ik-1
	xdebutperiode=minx+(in-1)*deltax
	xfinperiode=minx+in*deltax
 	if nacum gt 0 then begin
  	if idebutsample eq iborneinf then ydebutperiode=y(idebutsample) else begin
 	ydebutperiode=y(iborneinf)+(xdebutperiode-x(iborneinf))/ $
	(x(idebutsample)-x(iborneinf))*(y(idebutsample)-y(iborneinf))
	end
 	if ibornesup eq ifinsample then yfinperiode=y(ifinsample) else begin
	yfinperiode=y(ifinsample)+(xfinperiode-x(ifinsample))/ $
	(x(ibornesup)-x(ifinsample))*(y(ibornesup)-y(ifinsample))
	end
  	poidstotal=nacum*deltasamplex+xfinperiode-x(ifinsample)+$
	x(idebutsample)-xdebutperiode
	z(in-1)=(accum*deltasamplex+ydebutperiode*(x(idebutsample)- $
	xdebutperiode)+yfinperiode*(xfinperiode-x(ifinsample)))/poidstotal
   	endif
 		goto,fin
	endif
	;nous ne somme pas a la fin du sample et n'avons pas depasse une borne
 	goto,debutaccum
				endif else begin
	;nous venons de depasser une borne
 if in gt n-1 then print,'pb dep z'
	iz(in)=ik
	if nacum gt 0 then begin
	;calcul ordinaire de la moyenne entre deux bornes
	iborneinf=iz(in-1)-1
	if iborneinf  lt 0 then iborneinf=0
 	ibornesup=ik
	idebutsample=ik-nacum
	ifinsample=ik-1
	xdebutperiode=minx+(in-1)*deltax
	xfinperiode=minx+in*deltax
   	if idebutsample eq iborneinf then ydebutperiode=y(idebutsample) else begin
 	ydebutperiode=y(iborneinf)+(xdebutperiode-x(iborneinf))/ $
	(x(idebutsample)-x(iborneinf))*(y(idebutsample)-y(iborneinf))
	end
 	if ibornesup eq ifinsample then yfinperiode=y(ifinperiode) else begin
	yfinperiode=y(ifinsample)+(xfinperiode-x(ifinsample))/ $
	(x(ibornesup)-x(ifinsample))*(y(ibornesup)-y(ifinsample))
	end
 	poidstotal=nacum*deltasamplex+xfinperiode-x(ifinsample)+x(idebutsample)- $
	xdebutperiode
	z(in-1)=(accum*deltasamplex+ydebutperiode*(x(idebutsample)- $
	xdebutperiode)+yfinperiode*(xfinperiode-x(ifinsample)))/poidstotal
  		in=in+1
 			   endif else begin
		;nous venons de depasser une borne mais pas de points recense entre
	iborneinf=iz(in-1)-1
 	if iborneinf  lt 0 then iborneinf=0
 	ibornesup=ik
	xabcisse=minx+(in-0.5)*deltax
	z(in-1)=y(iborneinf)+(y(ibornesup)-y(iborneinf))/(x(ibornesup)-$
	x(iborneinf))*(xabcisse-x(iborneinf))
 		in=in+1
 			   endelse
	if ik lt n_elements(x) then goto,debutloop
		endelse
;probleme non compris : pansement
fin:if z(n) eq 0 then z(n)=z(n-1)
zx(*,1)=z
return,zx
end



;extraction de la distribution des returns et de
;la densite associe d un historique normalise
;n est le nombre de points associe au lissage la densite
;typiquement n=nombre de point/5
function show_densite,f,n,a1,b1
x1=f(4,*)
x=x1(a1:b1)
nb=n_elements(x)
window,0
plot,x
s=fltarr(nb-1)
for i=0,nb-2 do s(i)=alog(x(i+1)/x(i))
c=calculdistrib(s)
window,1
plot,c(*,0),c(*,1)
e=echantillonne(c(*,0),c(*,1),n)
window,2
plot,e(*,0),e(*,1)
dens=deriv(e(*,0),e(*,1))
plot,e(*,0),dens
mean=total(e(*,0)*dens)/total(dens)
variance=total((e(*,0)-mean)^2*dens)/total(dens)
print,'moyenne,sigma'
print,mean,sqrt(variance)
return,e
end



