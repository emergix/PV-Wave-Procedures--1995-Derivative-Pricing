function avg,a
n=n_elements(a)
return,total(a)/n
end
   
pro etude_volat,z,delai
;etude des fichier de type .histo
;delai est en jour reel pour le calcul de la volatilite
;datedebut est au format jjmmaa
;delaiwarrant est en annee
p='~olivier/'+z+'.histo'
openr,1,p
datecac=fltarr(2000)
opencac=datecac & highcac=datecac & lowcac=datecac & closecac=datecac
i=0
while not eof(1) do begin
readf,1,dd,oo,hh,bb,cc
datecac(i)=dd & opencac(i)=oo & highcac(i)=hh & lowcac(i)=bb & closecac(i)=cc
i=i+1
endwhile
print,i-1,'record'
close,1
nbhisto=i-1
jdebut=0
nbvol=0
vecvol=fltarr(10000)
j=0
while (j lt nbhisto) do begin
while ((jj2(datecac(j))-jj2(datecac(jdebut))) le delai and j lt nbhisto) do j=j+1
if j lt nbhisto then begin
jfin=j
nb=jfin-jdebut+1
ret=fltarr(nb)
for i=jdebut,jfin-1 do ret(i-jdebut)=alog(closecac(i+1)/closecac(i))
retbar=avg(ret)
u=0.
for i=0,nb-1 do u=u+(ret(i)-retbar)^2
vol=sqrt(u)
vecvol(jdebut)=vol/sqrt(delai/365.)
endif
jdebut=jdebut+1
endwhile
w=fltarr(jdebut-1)
w(0:jdebut-2)=vecvol(0:jdebut-2)
x=findgen(jdebut-1)/220.
window,0,title='volatilite a '+string(delai)+' jours reels'
plot,x,w
meanw=avg(w)
print,'volatilite moyenne a '+string(delai)+' jours (% annuel) : '+string(100*meanw)
print,'ecart-type de cette volatilite moyenne (% annuel): '+string(100*stdev(w))
oplot,x,replicate(meanw,n_elements(w)),color=112
jdebut=nbhisto-10
expvol=fltarr(jdebut+1)
x=expvol
while jdebut ge 0 do begin
ret=fltarr(nbhisto-jdebut)
for i=jdebut,nbhisto-1 do ret(i-jdebut)=alog(closecac(i+1)/closecac(i))
retbar=avg(ret)
u=0.
for i=0,nbhisto-jdebut-1 do u=u+(ret(i)-retbar)^2
vol=sqrt(u)
expvol(jdebut)=vol/sqrt((jj2(datecac(nbhisto-1))-jj2(datecac(jdebut)))/365.)
jdebut=jdebut-1
endwhile
for i=0,nbhisto-11 do x(i)=(jj2(datecac(nbhisto-1))-jj2(datecac(i)))/365.
window,1,title='volatilite arriere a duree croissante annualisee'
plot,x,expvol
end

pro etude_volat1,z,delai
;etude des fichier de type .print
;delai est en jour reel pour le calcul de la volatilite
;datedebut est au format jjmmaa
;delaiwarrant est en annee
p='/home/p6379j/'+z+'.print'
openr,1,p
datecac=fltarr(2000)
opencac=datecac & highcac=datecac & lowcac=datecac & closecac=datecac
i=0
while not eof(1) do begin
readf,1,dd,hehe,oo,hh,bb,cc,vv,nn
datecac(i)=dd & opencac(i)=oo & highcac(i)=hh & lowcac(i)=bb & closecac(i)=cc
i=i+1
endwhile
print,i-1,'record'
close,1
nbhisto=i-1
jdebut=0
nbvol=0
vecvol=fltarr(10000)
j=0
while (j lt nbhisto) do begin
while ((jj3(datecac(j))-jj3(datecac(jdebut))) le delai and j lt nbhisto) do j=j+1
if j lt nbhisto then begin
jfin=j
nb=jfin-jdebut+1
ret=fltarr(nb)
for i=jdebut,jfin-1 do ret(i-jdebut)=alog(closecac(i+1)/closecac(i))
retbar=avg(ret)
u=0.
for i=0,nb-1 do u=u+(ret(i)-retbar)^2
vol=sqrt(u)
vecvol(jdebut)=vol/sqrt(delai/365.)
endif
jdebut=jdebut+1
endwhile
w=fltarr(jdebut-1)
w(0:jdebut-2)=vecvol(0:jdebut-2)
x=findgen(jdebut-1)/220.
window,0,title='volatilite a '+string(delai)+' jours reels'
plot,x,w
meanw=avg(w)
print,'volatilite moyenne a '+string(delai)+' jours (% annuel) : '+string(100*meanw)
print,'ecart-type de cette volatilite moyenne (% annuel): '+string(100*stdev(w))
oplot,x,replicate(meanw,n_elements(w)),color=112
jdebut=nbhisto-10
expvol=fltarr(jdebut+1)
x=expvol
while jdebut ge 0 do begin
ret=fltarr(nbhisto-jdebut)
for i=jdebut,nbhisto-1 do ret(i-jdebut)=alog(closecac(i+1)/closecac(i))
retbar=avg(ret)
u=0.
for i=0,nbhisto-jdebut-1 do u=u+(ret(i)-retbar)^2
vol=sqrt(u)
expvol(jdebut)=vol/sqrt((jj3(datecac(nbhisto-1))-jj3(datecac(jdebut)))/365.)
jdebut=jdebut-1
endwhile
for i=0,nbhisto-11 do x(i)=(jj3(datecac(nbhisto-1))-jj3(datecac(i)))/365.
window,1,title='volatilite arriere a duree croissante annualisee'
plot,x,expvol
end

