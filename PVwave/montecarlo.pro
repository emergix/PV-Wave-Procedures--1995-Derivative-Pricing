  ;
;
;
;
;
;
;  fonctions de calculs optionnels d'inspiration Monte Carlo
;
;
;
;
;
;
;
;
;
;
function regresgauss,s
;suppose que x est un vecteur gaussien et retourne 
;un vecteur [mean , sigma] des coefficients gaussien associes aux rendements
n=n_elements(s)
mean=avg(s)
sigma=stdev(s)
print,'moyenne estimee des rendements : '+string(mean)
print,'ecart type standart des rendements :'+string(sigma)
return,[mean,sigma]
end


function genegauss,mean,sigma,r,nbjours,ik
;retourne une liste de prix dont les rendements sont gaussiens
common rando,table,rflag
if rflag eq 0 then begin
	rendements=randomn(seed,nbjours)*sigma
endif else begin
	rendements=table(long(ik)*nbjours:long(ik+1)*nbjours-1)*sigma
endelse
f=fltarr(nbjours)
s=1.
for i=0,nbjours-1 do begin
s=s*exp((mean+r-sigma^2/2.)+rendements(i))
f(i)=s
endfor
return,f
end



function testgaussKS,s,mean,sigma
;effectue un test de kolmogorov-smirnov gaussien sur le vecteur donnee en argument
nb=n_elements(s)+1
c=calculdistrib(s)
distribgauss=gaussint((c(*,0)-mean)/sigma)
distance=abs(c(*,1)-distribgauss)
Dn=max(distance)
np=n_elements(distance)
seuil1=0.01
seuil5=0.05
dn1=1.6276/sqrt(np)
dn5=1.3581/sqrt(np)
KS1=Dn-dn1
KS5=Dn-dn5
return,[KS1,KS5]
end

;fonction qui teste le generateur gaussien
pro testgenegauss,mean,sigma,r,ik
common rando,table,rflag
x=genegauss(mean,sigma,r,500,0)
vec=regresgauss(logprix(x))
a=testgaussKS(logprix(x),mean,sigma)
nb=n_elements(logprix(x))+1
c=calculdistrib(logprix(x))
distribgauss=gaussint((c(*,0)-mean)/sigma)
window,31,title='distribution empirique/theorique pour le cas gaussien'
plot,c(*,0),c(*,1)
oplot,c(*,0),distribgauss,color=112
end

pro display_delta,f,deltaf,cout,coutexercicefinal,investinitial,result
	set_plot,'PS'
	device,/landscape
	plot,f,title='cours du support' ,$
		xtitle='cout final total='+string(result)+' / exercice final='+string(coutexercicefinal)+$
		'  / investissement initial ='+string(investinitial)
	plot,deltaf,title='delta'
	plot,cout,title='cout du gamma'
	device,/close_file
	spawn,'lpr wave.ps'
	set_plot,'X'
end

pro display_distrib,a,nb,titre,titre2,pr
a1=fltarr(n_elements(a))
for i=0,n_elements(a)-1 do a1(i)=avg(a(0:i))
b=calculdistrib(a)
x=b(*,0)
w=b(*,1)
z=echantillonne(x,w,int(nb))
zx=z(*,0)
zw=z(*,1)
q=deriv(zx,zw)
vec2=regresgauss(a)
mean=vec2(0)
sigma=vec2(1)
distribgauss=gaussint((x-mean)/sigma)
ra=testgaussKS(x,mean,sigma)
KS1=ra(0)
KS5=ra(1)
if KS1 le 0 then rep1=' = acceptation' else rep1=' = rejet'
if KS5 le 0 then rep5=' = acceptation' else rep5=' = rejet'
if pr eq 1 then  begin
	set_plot,'PS'
	device,/landscape
	plot,a1,title='convergence du processus de MC :'+'mean='+string(mean)+'   sigma='+string(sigma) ,$
		xtitle=titre,$
		subtitle=titre2
	plot,x,w,title='distribution des cout d emission :'+'mean='+string(mean)+'   sigma='+string(sigma),$
		xtitle=titre,$
		subtitle=titre2,$
		Ytitle='niveau 1% :'+string(KS1)+rep1+'  /niveau 5% :'+string(KS5)+rep5
	oplot,[mean,mean],[0,max(w)]
	oplot,x,distribgauss,color=112
	plot,a,replicate(-0.1,n_elements(a)),psym=3,yrange=[-0.25,1], $
		title='densite des couts d emission :'+'mean='+string(mean)+'   sigma='+string(sigma)  ,$
		xtitle=titre,$
		subtitle=titre2
	oplot,zx,q,color=112
	device,/close_file
	spawn,'lpr wave.ps'
endif
set_plot,'X'
if pr lt 2 then  begin
	window,0,title='convergence du prcocessus de monte carlo'
	for i=0,n_elements(a)-1 do a1(i)=avg(a(0:i))
	plot,a1,title='convergence du processus de MC :',$
		xtitle=titre
	window,1,title='distribution des cout d emission'
	plot,x,w,title='distribution des cout d emission :',$
		xtitle=titre
	oplot,[mean,mean],[0,max(w)]
	oplot,x,distribgauss,color=112
	window,2,title='densite des cout d emission'
	plot,a,replicate(-0.1,n_elements(a)),psym=3,yrange=[-0.25,1], $
		title='densite des couts d emission :' ,$
		xtitle=titre
	oplot,zx,q,color=112
	print,'---- on regarde l hypothese gaussienne'
	print,'kolmogorov-smirnov gaussien au niveau 1% :'+string(KS1)+rep1
	print,'kolmogorov-smirnov gaussien au niveau 5% :'+string(KS5)+rep5
endif
end


function GEScall,f,df,k0,taux,sig,debut,fin,tauxfrais
nb=fin-debut+1
cf=f(debut:fin)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
 for i=0,nb-1 do deltaf(i)=call_delta(cf(i),k0,bdates(i),taux,sig)
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
end
if cf(nb-1) gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) $
 else coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)
return,result
end


function MCcall,s01,k01,t01,taux,sig,nbj,nbs,iflg,rflag1
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GEScall(f,df,k0,rj,volj,0,nbj-1,0.)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
mean=avg(vec)
titre='call simple , s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
display_distrib,vec,nbs/10,titre,' ',iflg
return,mean
end


function GESPcall,f,df,k0,taux,sig,debut,fin,tauxfrais,palier,ndel
nb=fin-debut+1
cf=f(debut:fin)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
palierflag=fltarr(nb)
if cf(0) ge (palier+k0) then palierflag(0)=1
for i=1,nb-1 do $
	if cf(i) ge (palier+k0) then palierflag(i)= max([palierflag(i-1),1]) else $
		palierflag(i)=palierflag(i-1)
if palierflag(nb-1) eq 1 then palierlevel=palier+k0
for i=0,nb-1 do if palierflag(i) eq 0 then $ 
	deltaf(i)= CREcall_delta(cf(i),k0,bdates(i),taux,sig,ndel,palier) else $
	deltaf(i)= call_delta(cf(i),palierlevel,bdates(i),taux,sig)
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
end
if cf(nb-1) gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) $
 else coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
coutpalier=0
if (palierflag(nb-1) eq 1) and (cf(nb-1) lt palier+k0) then $
	coutpalier=palierlevel*exp(taux*(df(0)-df(nb-1)))
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)+coutpalier
display_delta,f,deltaf,cout,coutexercicefinal,investinitial,result
return,result
end

function MCPcall,s01,k01,t01,taux,sig,palier1,nbj,nbs,ndel,iflg,rflag1
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
palier=float(palier1)
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESPcall(f,df,k0,rj,volj,0,nbj-1,0.,palier,ndel)
vec(j)=cout
if j mod 1 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
mean=avg(vec)
titre='call 1 palier , s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
display_distrib,vec,nbs/10,titre,'palier='+string(palier1),iflg
return,mean
end



function GESPLcall,f,df,k0,taux,sig,debut,fin,tauxfrais,palier_list,ndel
nb=fin-debut+1
cf=f(debut:fin)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
palierflag=fltarr(nb)
pw=where(palier_list+k0 le cf(0),cw)
if cw eq 0 then palierflag(0)=0 else palierflag(0)=max(pw)+1
for i=1,nb-1 do begin
	pw=where(palier_list+k0 le cf(i),cw)
	if cw eq 0 then palierfl=0 else palierfl=max(pw)+1
	palierflag(i)= max([palierflag(i-1),palierfl])  
	end
if palierflag(nb-1) ge 1 then palierlevel=palier_list(palierflag(nb-1)-1) else palierlevel=0
for i=0,nb-1 do if palierflag(i) lt n_elements(palier_list) then $ 
	deltaf(i)= CREcall_delta(cf(i),k0+palier_list(palierflag(i)-1), $
		bdates(i),taux,sig,ndel,palier_list(palierflag(i))) else $
	deltaf(i)= call_delta(cf(i),k0+palier_list(palierflag(i)-1),bdates(i),taux,sig)
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
end
if cf(nb-1) gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) $
 else coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
coutpalier=0
if (palierflag(nb-1) ge 1) and (cf(nb-1) lt palierlevel) then $
	coutpalier=palierlevel*exp(taux*(df(0)-df(nb-1)))
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)+coutpalier
return,result
end


function MCPLcall,s01,k01,t01,taux,sig,palier_list1,nbj,nbs,ndel,iflg,rflag1
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
np=n_elements(palier_list1)
palier_list2=fltarr(np)
for i=0,np-1 do palier_list2(i)=float(palier_list1(i))
palier_list=palier_list2(sort(palier_list2))
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESPLcall(f,df,k0,rj,volj,0,nbj-1,0.,palier_list,ndel)
vec(j)=cout
if j mod 1 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
mean=avg(vec)
titre='call paliers , s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2='Palier_list='+string(palierlist1)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,mean
end




function GESYLcall,f,df,k0,taux,sig,debut,fin,tauxfrais,yield_list,nbjdiv_list
;plot,df(debut:fin),f(debut:fin)
;calcule le cout d une option le long d un chemin particulier
;avec une liste de yields 
nb=fin-debut+1
cf=fltarr(nb)
for i=0,nb-1 do $
cf(i)=f(i+debut)*yield_dec(df(debut),df(i+debut),yield_list,nbjdiv_list)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
 for i=0,nb-1 do deltaf(i)=YLcall_delta(cf(i),k0,bdates(i),taux,sig, $
		yield_list,nbjdiv_list-df(i+debut))
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
dividends=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
decap=yield_dec(df(i+debut-1),df(i+debut),yield_list,nbjdiv_list)
dividends(i)=dividends(i-1)+cf(i)*deltaf(i)*(1.-decap)
end
if cf(nb-1) gt k0 then $
coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)-dividends(nb-1)
return,result
end

function MCYLcall,s01,k01,t01,taux,sig,yield_list,datediv_list,nbj,nbs,iflg,rflag1
;calcule le prix d un call avec liste de reset et liste de yield 
;par la methode monte carlo : nbs = nombre de tirage de montecarlo 
; nbj = nombre de jour de bourse , sig et taux etant la volatilite 
;et le taux par an
;datediv_list etant la liste des distance en jours separant
; l origine des distribution de dividendes
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESYLcall(f,df,k0,rj,volj,0,nbj-1,0.,yield_list,datediv_list*nbj/t0)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
BS=YLcall(s0,k0,t0,taux,sig,yield_list,datediv_list)
mean=avg(vec)
print,'le prix moyen de l option est '+string(mean)
print,'alors que Black & scholes = '+string(BS)
titre='call avec yields , s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2='Yield_list='+mapstring(yield_list)+' / date_list='+mapstring(datediv_list)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,mean
end

function GESYLPcall,f,df,k0,taux,sig,debut,fin,tauxfrais,yield_list, $
						nbjdiv_list,palier,ndel
;plot,df(debut:fin),f(debut:fin)
;calcule le cout d une option le long d un chemin particulier
;avec une liste de yields 
nb=fin-debut+1
cf=fltarr(nb)
for i=0,nb-1 do $
cf(i)=f(i+debut)*yield_dec(df(debut),df(i+debut),yield_list,nbjdiv_list)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
palierflag=fltarr(nb)
if cf(0) ge (palier+k0) then palierflag(0)=1
for i=1,nb-1 do $
	if cf(i) ge (palier+k0) then palierflag(i)= max([palierflag(i-1),1]) else $
		palierflag(i)=palierflag(i-1)
if palierflag(nb-1) eq 1 then palierlevel=cf(min(where(palierflag eq 1)))
for i=0,nb-1 do if palierflag(i) eq 0 then $ 
	deltaf(i)= CRYLEcall_delta(cf(i),k0,bdates(i),taux,sig,ndel,yield_list, $
			nbjdiv_list-df(i+debut),palier) else $
	deltaf(i)= YLcall_delta(cf(i),k0,bdates(i),taux,sig, $
		yield_list,nbjdiv_list-df(i+debut))
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
dividends=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
decap=yield_dec(df(i+debut-1),df(i+debut),yield_list,nbjdiv_list)
dividends(i)=dividends(i-1)+cf(i)*deltaf(i)*(1.-decap)
end
if cf(nb-1) gt k0 then $
coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
coutpalier=0
if (palierflag(nb-1) eq 1) and (cf(nb-1) lt palier+k0) then $
	coutpalier=palierlevel*exp(taux*(df(0)-df(nb-1)))
result=coutexercicefinal+investinitial+cout(nb-1)+ $
	fraistrans(nb-1)-dividends(nb-1)+coutpalier
return,result
end
function MCYLPcall,s01,k01,t01,taux,sig,yield_list,datediv_list,palier1,nbj,nbs,ndel,iflg,rflag1
;calcule le prix d un call avec liste de reset et liste de yield 
;par la methode monte carlo : nbs = nombre de tirage de montecarlo 
; nbj = nombre de jour de bourse , sig et taux etant la volatilite 
;et le taux par an
;datediv_list etant la liste des distance en jours separant
; l origine des distribution de dividendes
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
palier=palier1
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESYLPcall(f,df,k0,rj,volj,0,nbj-1,0.,yield_list,datediv_list*nbj/t0,palier,ndel)
vec(j)=cout
if j mod 1 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
BS=YLcall(s0,k0,t0,taux,sig,yield_list,datediv_list)
mean=avg(vec)
print,'le prix moyen de l option est '+string(mean)
print,'alors que Black & scholes = '+string(BS)
titre='call avec yield et 1 palier , s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2='Yield_list='+string(yield_list)+' / date_list='+string(datediv_list)+' Palier='+string(palier1)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,mean
end


function GESYLPLcall,f,df,k0,taux,sig,debut,fin,tauxfrais,yield_list, $
						nbjdiv_list,palier_list,ndel
;plot,df(debut:fin),f(debut:fin)
;calcule le cout d une option le long d un chemin particulier
;avec une liste de yields 
nb=fin-debut+1
cf=fltarr(nb)
for i=0,nb-1 do $
cf(i)=f(i+debut)*yield_dec(df(debut),df(i+debut),yield_list,nbjdiv_list)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
palierflag=fltarr(nb)
pw=where(palier_list+k0 le cf(0),cw)
if cw eq 0 then palierflag(0)=0 else palierflag(0)=max(pw)+1
for i=1,nb-1 do begin
	pw=where(palier_list+k0 le cf(i),cw)
	if cw eq 0 then palierfl=0 else palierfl=max(pw)+1
	palierflag(i)= max([palierflag(i-1),palierfl])  
	end
if palierflag(nb-1) ge 1 then palierlevel=palier_list(palierflag(nb-1)-1) else $
		 palierlevel=0
for i=0,nb-1 do begin
	if palierflag(i) gt 0 then begin
	if palierflag(i) lt n_elements(palier_list) then $
	deltaf(i)= CRYLEcall_delta(cf(i),k0+palier_list(palierflag(i)-1), $
			bdates(i),taux,sig,ndel,yield_list, $
			nbjdiv_list-df(i+debut),palier_list(palierflag(i))) else $
	deltaf(i)= YLcall_delta(cf(i),k0+palier_list(palierflag(i)-1),bdates(i),taux,sig, $
		yield_list,nbjdiv_list-df(i+debut))
	endif else begin
	deltaf(i)= CRYLEcall_delta(cf(i),k0,bdates(i),taux,sig,ndel,yield_list, $
			nbjdiv_list-df(i+debut),palier_list(0))
	endelse
	end
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
dividends=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
decap=yield_dec(df(i+debut-1),df(i+debut),yield_list,nbjdiv_list)
dividends(i)=dividends(i-1)+cf(i)*deltaf(i)*(1.-decap)
end
if cf(nb-1) gt k0 then $
coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
coutpalier=0
if (palierflag(nb-1) ge 1) and (cf(nb-1) lt palierlevel) then $
	coutpalier=palierlevel*exp(taux*(df(0)-df(nb-1)))
result=coutexercicefinal+investinitial+cout(nb-1)+ $
	fraistrans(nb-1)-dividends(nb-1)+coutpalier
return,result
end

function MCYLPLcall,s01,k01,t01,taux,sig,yield_list, $
		datediv_list,palier_list1,nbj,nbs,ndel,iflg,rflag1
;calcule le prix d un call avec liste de reset et liste de yield 
;par la methode monte carlo : nbs = nombre de tirage de montecarlo 
; nbj = nombre de jour de bourse , sig et taux etant la volatilite 
;et le taux par an
;datediv_list etant la liste des distance en jours separant
; l origine des distribution de dividendes
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
np=n_elements(palier_list1)
palier_list2=fltarr(np)
for i=0,np-1 do palier_list2(i)=float(palier_list1(i))
palier_list=palier_list2(sort(palier_list2))
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESYLPLcall(f,df,k0,rj,volj,0,nbj-1,0.,yield_list, $
		datediv_list*nbj/t0,palier_list,ndel)
vec(j)=cout
if j mod 1 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
BS=YLcall(s0,k0,t0,taux,sig,yield_list,datediv_list)
mean=avg(vec)
print,'le prix moyen de l option est '+string(mean)
print,'alors que Black & scholes = '+string(BS)
titre='call avec yield et paliers , s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2='Yield_list='+mapstring(yield_list)+' / date_list='+mapstring(datediv_list)+' Palier_list='+string(palier_list1)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,mean
end




function GESRScall,f,df,k0,taux,sig,debut,fin,tauxfrais,nbjreset_list1
nbjreset_list=nbjreset_list1(sort(nbjreset_list1))
nr=n_elements(nbjreset_list1)
creset_list=fltarr(nr)
for i=0,nr-1 do begin
zw=where(df le nbjreset_list(i),count)
if count gt 0 then creset_list(i)=f(max(zw)) else creset_list(i)=0.
endfor
nb=fin-debut+1
cf=f(debut:fin)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
 for i=0,nb-1 do deltaf(i)= $
RScall_delta(cf(i),k0,bdates(i),taux,sig,nbjreset_list-df(i+debut),creset_list)
deltareset=fltarr(4,nr)
;for i=0,nr-1 do begin
;zw=where(df le nbjreset_list(i),count)
;if count gt 0 then deltareset(0,i)=deltaf(max(zw)) 
;if count gt 0 then deltareset(1,i)=deltaf(max(zw)-1) 
;if count gt 0 then deltareset(2,i)=cf(max(zw)) 
;if count gt 0 then deltareset(3,i)=cf(max(zw)-1) 
;endfor
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
end
nj1=df(debut)
nj2=df(fin)
indexs=where((nbjreset_list ge nj1) and (nbjreset_list le nj2))
cmax=max(f(nbjreset_list(indexs)))
cmax2=max([cmax,cf(fin)])
if cmax2 gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else $
 coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
if cf(fin) ge cmax then begin
	coutreset=0
endif else begin
	if cmax gt k0 then begin
		coutreset=cmax*exp(taux*(df(0)-df(nb-1)))
	endif else begin
		coutreset=0
	endelse
endelse
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)+coutreset
;print,cmax,cout(nb-1),deltareset
;print,[cf(fin),coutreset,coutexercicefinal+investinitial+coutreset,result]
;window,0,title='cours'
;plot,cf
;window,1,title='delta'
;plot,deltaf
;window,3,title='cout de la couverture'
;plot,investinitial+cout
;print,'debit client='+string(coutexercicefinal+coutreset)
return,result
end

 
function MCRScall,s01,k01,t01,taux,sig,datereset_list,nbj,nbs,iflg,rflag1
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESRScall(f,df,k0,rj,volj,0,nbj-1,0.,datereset_list*nbj/t0)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+'t='+systime()
endif
endfor
BS=call(s0,k0,t0,taux,sig)
mean=avg(vec)
print,'le prix moyen de l option est '+string(mean)
nr=n_elements(datereset_list)
sumoption=0.
if nr gt 0 then begin 
for i=0,nr do begin
case 1 of 
 	(i eq 0): deltat=datereset_list(0)
	(i eq nr):deltat=t0-datereset_list(i-1)
	else:deltat=datereset_list(i)-datereset_list(i-1)
endcase
opt=call(s0,k0,deltat,taux,sig)
if (i eq 0) then deltat=0 else deltat=datereset_list(i-1)
act=exp(taux*deltat)
print,'option de rang '+string(i)+'='+string(opt)
sumoption=sumoption+opt/act
endfor
endif
if nr eq 0 then sumoption=BS
print, 'le black and sholes somme donne ='+string(sumoption)
print,'alors que black and sholes normal = '+string(BS)
titre='call avec reset, s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2=' Datereset_list='+mapstring(datereset_list)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,mean
end

 
 function GESRScall1,f,df,k0,taux,sig,debut,fin,tauxfrais,nbjreset_list1,nbjreset_list3
;calcule le resultat de gestion avec la liste nbjreset_list1 en calculant le delta avec nbjreset_list3
nbjreset_list=nbjreset_list1(sort(nbjreset_list1))
nbjreset_list2=nbjreset_list3(sort(nbjreset_list3))
nr=n_elements(nbjreset_list1)
creset_list=fltarr(nr)
for i=0,nr-1 do begin
zw=where(df le nbjreset_list(i),count)
if count gt 0 then creset_list(i)=f(max(zw)) else creset_list(i)=0.
endfor
nr2=n_elements(nbjreset_list3)
creset_list2=fltarr(nr2)
for i=0,nr2-1 do begin
zw2=where(df le nbjreset_list2(i),count2)
if count2 gt 0 then creset_list2(i)=f(max(zw2)) else creset_list2(i)=0.
endfor
nb=fin-debut+1
cf=f(debut:fin)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
 for i=0,nb-1 do deltaf(i)= $
RScall_delta(cf(i),k0,dn-df(i+debut),taux,sig,nbjreset_list2-df(i+debut),creset_list2)
deltareset=fltarr(4,nr)
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
end
nj1=df(debut)
nj2=df(fin)
indexs=where((nbjreset_list ge nj1) and (nbjreset_list le nj2))
cmax=max(f(nbjreset_list(indexs)))
cmax2=max([cmax,cf(fin)])
if cmax2 gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else $
 coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
if cf(fin) ge cmax then begin
	coutreset=0
endif else begin
	if cmax gt k0 then begin
		coutreset=cmax*exp(taux*(df(0)-df(nb-1)))
	endif else begin
		coutreset=0
	endelse
endelse
pindexs=where((nbjreset_list2 ge nj1) and (nbjreset_list2 le nj2))
pcmax=max(f(nbjreset_list2(pindexs)))
pcmax2=max([pcmax,cf(fin)])
if pcmax2 gt k0 then pcoutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else $
 pcoutexercicefinal=0.
if cf(fin) ge pcmax then begin
	pcoutreset=0
endif else begin
	if pcmax gt k0 then begin
		pcoutreset=pcmax*exp(taux*(df(0)-df(nb-1)))
	endif else begin
		pcoutreset=0
	endelse
endelse
coutoption=fltarr(nr2+1)
gainoption=fltarr(nr2+1)
for j=0,nr2-1 do begin
  if j eq 0 then investini=deltaf(0)*cf(0) else investini=deltaf(nbjreset_list2(j-1))*cf(nbjreset_list2(j-1))*exp(taux*(df(0)-df(nbjreset_list2(j-1))))
  if j eq 0 then coutgamma=cout(nbjreset_list2(0)) else coutgamma=cout(nbjreset_list2(j))-cout(nbjreset_list2(j-1))
  if deltaf(nbjreset_list2(j)-1) gt 0.5 then coutexercice=-k0*exp(taux*(df(0)-df(nbjreset_list2(j)))) else coutexercice=0
  if deltaf(nbjreset_list2(j)-1) gt 0.5 then gainoption(j)=deltaf(nbjreset_list2(j))*cf(nbjreset_list2(j))*exp(taux*(df(0)-df(nbjreset_list2(j)))) $
    else gainoption(j)=0
  coutoption(j)=investini+coutgamma+coutexercice
endfor
investini=deltaf(nbjreset_list2(nr2-1))*cf(nbjreset_list2(nr2-1))
coutgamma=cout(nb-1)-cout(nbjreset_list2(nr2-1))
if deltaf(nb-1) gt 0.99 then coutexercice=-k0*exp(taux*(df(0)-df(nb-1))) else coutexercice=0.
if deltaf(nb-1) gt 0.99 then gainoption(nr2)=deltaf(nb-1)*cf(nb-1)*exp(taux*(df(0)-df(nb-1))) else gainoption(nr2)=0
coutoption(nr2)=investini+coutgamma+coutexercice
print,'opt interm:',coutoption
print,'total cout',total(coutoption)
print,'gain option',gainoption
print,'total gain',total(gainoption)
print,'som=',total(coutoption)+total(gainoption)
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)+coutreset
presult=pcoutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)+pcoutreset
print,'rsdate=',nbjreset_list2-1
print,'rscf=',cf(nbjreset_list2-1)
print,'rsdelta fin=',deltaf(nbjreset_list2-1)
print,'rsdelta debut=',deltaf(df(0)),deltaf(nbjreset_list2)
print,'prsdate=',nbjreset_list-1
print,'prscf=',cf(nbjreset_list-1)
print,'coutreset=',coutreset,' coutexercicefinal',coutexercicefinal,' cout gamma=',cout(nb-1),' deltafinal=',deltaf(nb-1)
print,'pcoutreset=',pcoutreset,' pcoutexercicefinal',pcoutexercicefinal
print,'investinitial=',investinitial,' echeance=',-(df(0)-df(nb-1)),' cout faux delta vendu=',result,' cout vrai delta =',presult
return,result
end


function GESYLRScall,f,df,k0,taux,sig,debut,fin,tauxfrais, $
		yield_list,nbjdiv_list,nbjreset_list1
;calcule le cout d une option le long d un chemin particulier
;avec une liste de yields 
nbjreset_list=nbjreset_list1(sort(nbjreset_list1))
nr=n_elements(nbjreset_list1)
creset_list=fltarr(nr)
for i=0,nr-1 do begin
zw=where(df le nbjreset_list(i),count)
if count gt 0 then creset_list(i)=f(max(zw)) else creset_list(i)=0.
endfor
nb=fin-debut+1
cf=fltarr(nb)
for i=0,nb-1 do $
cf(i)=f(i+debut)*yield_dec(df(debut),df(i+debut),yield_list,nbjdiv_list)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
 for i=0,nb-1 do deltaf(i)=YLRScall_delta(cf(i),k0,bdates(i),taux,sig, $
yield_list,nbjdiv_list-df(i+debut),nbjreset_list-df(i+debut),creset_list)
coutmarg=fltarr(nb)
prixmoy=fltarr(nb)
prixmoy(0)=cf(0)
cout=fltarr(nb)
fraistrans=fltarr(nb)
dividends=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
if deltaf(i) ne 0 then prixmoy(i)=(deltaf(0)*cf(0)+cout(i))/deltaf(i) else prixmoy(i)=prixmoy(i-1)

fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
decap=yield_dec(df(i+debut-1),df(i+debut),yield_list,nbjdiv_list)
dividends(i)=dividends(i-1)+cf(i)*deltaf(i)*(1.-decap)
end
nj1=df(debut)
nj2=df(fin)
indexs=where((nbjreset_list ge nj1) and (nbjreset_list le nj2))
cmax=max(f(nbjreset_list(indexs)))
cmax2=max([cmax,cf(fin)])
if cmax2 gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else $
 coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
if cf(fin) ge cmax then begin
	coutreset=0
endif else begin
	if cmax gt k0 then begin
		coutreset=cmax*exp(taux*(df(0)-df(nb-1)))
	endif else begin
		coutreset=0
	endelse
endelse
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)- $
			dividends(nb-1)+coutreset
;window,1,title='cours'
;plot,df(debut:fin),cf(debut:fin)
;window,2,title='delta'
;plot,df(debut:fin),deltaf(debut:fin)
;window,3,title='gammaneg'
;plot,df(debut:fin),cout+investinitial
;window,4,title='cout moyen'
;plot,df(debut:fin),prixmoy(debut:fin-1),yrange=[0,+10*k0]
;print,[cf(fin),cmax,coutreset,dividends(nb-1),coutexercicefinal,investinitial,cout(nb-1),result]
return,result
end



function MCYLRScall,s01,k01,t01,taux,sig,yield_list, $
		datediv_list,datereset_list,nbj,nbs,iflg,rflag1
;calcule le prix d un call avec liste de reset et liste de yield 
;par la methode monte carlo : nbs = nombre de tirage de montecarlo 
; nbj = nombre de jour de bourse , sig et taux etant la volatilite 
;et le taux par an
;datediv_list etant la liste des distance en jours separant
; l origine des distribution de dividendes
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESYLRScall(f,df,k0,rj,volj,0,nbj-1,0.,yield_list, $
		datediv_list*nbj/t0,datereset_list*nbj/t0)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
mean=avg(vec)
print,'le prix moyen de l option est '+string(mean)
nr=n_elements(datereset_list)
sumoption=0.
if nr gt 0 then begin 
for i=0,nr do begin
case 1 of 
 	(i eq 0): deltat=datereset_list(0)
	(i eq nr):deltat=t0-datereset_list(i-1)
	else:deltat=datereset_list(i)-datereset_list(i-1)
endcase
opt=call(s0,k0,deltat,taux,sig)
if (i eq 0) then deltat=0 else deltat=datereset_list(i-1)
act=exp(taux*deltat)
print,'option de rang '+string(i)+'='+string(opt)
sumoption=sumoption+opt/act
endfor
endif
if iflg lt 2 then begin
	BS=YLcall(s0,k0,t0,taux,sig,yield_list,datediv_list)
	if nr eq 0 then sumoption=BS
	print, 'le black and sholes somme donne ='+string(sumoption)
	print,'alors que black and sholes normal = '+string(BS)
endif
titre='call avec yields et resets, s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2='Yield_list='+mapstring(yield_list)+' /date_list='+mapstring(datediv_list)+'Datereset_list='+mapstring(datereset_list)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,avg(vec)
end


function GESYLMLcall,f,df,k0,taux,sig,debut,fin,tauxfrais, $
			yield_list,nbjdiv_list,nbjmoyenne_list1
;plot,df(debut:fin),f(debut:fin)
;calcule le cout d une option le long d un chemin particulier
;avec une liste de yields 
nbjmoyenne_list=nbjmoyenne_list1(sort(nbjmoyenne_list1))
nr=n_elements(nbjmoyenne_list1)
cmoyenne_list=fltarr(nr)
for i=0,nr-1 do begin
zw=where(df le nbjmoyenne_list(i),count)
if count gt 0 then cmoyenne_list(i)=f(max(zw)) else cmoyenne_list(i)=0.
endfor
nb=fin-debut+1
cf=fltarr(nb)
for i=0,nb-1 do $
cf(i)=f(i+debut)*yield_dec(df(debut),df(i+debut),yield_list,nbjdiv_list)
dn=df(fin)
deltaf=fltarr(nb)
bdates=dn-df(debut:fin)
 for i=0,nb-1 do deltaf(i)=YLMLcall_delta(cf(i),k0,bdates(i),taux,sig, $
	yield_list,nbjdiv_list-df(i+debut),nbjmoyenne_list-df(i+debut),cmoyenne_list)
coutmarg=fltarr(nb)
cout=fltarr(nb)
fraistrans=fltarr(nb)
dividends=fltarr(nb)
portagetotal=fltarr(nb)
for i=1,nb-1 do  begin 
coutmarg(i)=cf(i)*(deltaf(i)-deltaf(i-1))*exp(taux*(df(0)-df(i)))
cout(i)=cout(i-1)+coutmarg(i)
fraistrans(i)=fraistrans(i-1)+tauxfrais*cf(i)*abs(deltaf(i)-deltaf(i-1))
decap=yield_dec(df(i+debut-1),df(i+debut),yield_list,nbjdiv_list)
dividends(i)=dividends(i-1)+cf(i)*deltaf(i)*(1.-decap)
end
nj1=df(debut)
nj2=df(fin)
indexs=where((nbjmoyenne_list ge nj1) and (nbjmoyenne_list le nj2))
cavg=avg(f(nbjmoyenne_list(indexs)))
cmax2=max([cavg,cf(fin)])
if cmax2 gt k0 then coutexercicefinal=-k0*exp(taux*(df(0)-df(nb-1))) else $
 coutexercicefinal=0.
investinitial=deltaf(0)*cf(0)
if cf(fin) ge cavg then begin
	coutmoyenne=0
endif else begin
	if cavg gt k0 then begin
		coutmoyenne=cavg*exp(taux*(df(0)-df(nb-1)))
	endif else begin
		coutmoyenne=0
	endelse
endelse
result=coutexercicefinal+investinitial+cout(nb-1)+fraistrans(nb-1)- $
	dividends(nb-1)+coutmoyenne
return,result
end

function MCYLMLcall,s01,k01,t01,taux,sig,yield_list,datediv_list, $
		datemoyenne_list,nbj,nbs,iflg,rflag1
;calcule le prix d un call avec liste de moyenne et liste de yield 
;par la methode monte carlo : nbs = nombre de tirage de montecarlo 
; nbj = nombre de jour de bourse , sig et taux etant la volatilite 
;et le taux par an
;datediv_list etant la liste des distance en jours separant
; l origine des distribution de dividendes
common rando,rtable,rflag
rflag=rflag1
if rflag eq 1 then rtable=lecture2R('random.normal.table',nbj*nbs)
s0=float(s01)
k0=float(k01)
t0=float(t01)
rj=(1+taux)^(t0/nbj)-1.
volj=sig*sqrt(t0/nbj)
vec=fltarr(nbs)
for j=0,nbs-1 do begin
f=s0*genegauss(0,volj,rj,nbj,j)
df=findgen(nbj)
cout=GESYLMLcall(f,df,k0,rj,volj,0,nbj-1,0.,yield_list, $
		datediv_list*nbj/t0,datemoyenne_list*nbj/t0)
vec(j)=cout
if j mod 10 eq 0. then begin
print,'step='+string(j)+'/'+string(nbs)+' option='+string(avg(vec(0:j)))+' t='+systime()
endif
endfor
mean=avg(vec)
print,'le prix moyen de l option est '+string(mean)
nr=n_elements(datemoyenne_list)
sumoption=0.
if nr gt 0 then begin 
for i=0,nr do begin
case 1 of 
 	(i eq 0): deltat=datemoyenne_list(0)
	(i eq nr):deltat=t0-datemoyenne_list(i-1)
	else:deltat=datemoyenne_list(i)-datemoyenne_list(i-1)
endcase
opt=call(s0,k0,deltat,taux,sig)
if (i eq 0) then deltat=0 else deltat=datemoyenne_list(i-1)
act=exp(taux*deltat)
print,'option de rang '+string(i)+'='+string(opt)
sumoption=sumoption+opt/act
endfor
endif
BS=YLcall(s0,k0,t0,taux,sig,yield_list,datediv_list)
if nr eq 0 then sumoption=BS
print, 'le black and sholes somme donne ='+string(sumoption)
print,'alors que black and sholes normal = '+string(BS)

titre='call avec yields et moyenne, s='+string(s01)+' k='+string(k01)+' t='+string(t01)+' r='+string(taux)+$
	' v='+string(sig)
titre2='Yield_list='+string(yield_list)+' / date_list='+mapstring(datediv_list)+' DateMoy._list='+string(datemoyenne_list)
display_distrib,vec,nbs/10,titre,titre2,iflg
return,mean
end


pro essaih,nbj,nbs
print,MCYLcall(35.,40.,7/12.,alog(1.05),0.3,[0.5/40.,0.5/40,0.5/40],[0.5/12.,3.5/12,6.5/12],nbj,nbs,1)
print,MCYLcall(40.,40.,7/12.,alog(1.05),0.3,[0.5/40.,0.5/40,0.5/40],[0.5/12.,3.5/12,6.5/12],nbj,nbs,1)
print,MCYLcall(45.,40.,7/12.,alog(1.05),0.3,[0.5/40.,0.5/40,0.5/40],[0.5/12.,3.5/12,6.5/12],nbj,nbs,1)
end


pro essaiRS,nbj,nbs
mcrs=MCRScall(100.,100.,4.,0.1,0.3,[1,2,3],nbj,nbs,0)
print,mcrs
end

pro essaiYLRS,nbj,nbs
mcrs=MCYLRScall(100.,100.,2.,0.1,0.3,[0.0217,0.0217],[0.5,1.5],[0.7,1.3],nbj,nbs,0)
print,mcrs
end

pro essaiYLML,nbj,nbs
mcrs=MCYLMLcall(100.,100.,2.,0.1,0.18,[0.0217,0.0217],[0.5,1.5],[0.7,1.3],nbj,nbs,0)
print,mcrs
end


 pro essaiE,nbj,nbs,ndel
mcrs=MCPcall(100.,100.,2.,0.1,0.3,100,nbj,nbs,ndel,0)
print,mcrs
end


pro essaip,nbj,nbs,ndel
mcp=MCYLPcall(100.,100.,2.,0.1,0.13,[0.0217,0.0217],[0.5,1.5],25.,nbj,nbs,ndel,0)
end
 ;***************************************************************************

;comparaison option europeene / monte carlo europeen
pro essaig,z
s=100.
k=100
t=1.
r=0.1
sigma=0.3
njour=200
nessai=100
bs=call(s,k,t,r,sigma)
mc=MCcall(s,k,t,r,sigma,njour,nessai)
print,[bs,mc]
end
