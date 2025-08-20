; a besoin de option2.pro et de psimplexe.pro
	

;
;
;
;
;
;
;
;            calcul du risque d un portefeuille
;        
;                  d options
;
;
;
;
;
;
 ; calcul du risque
; un portefeuille se presente sous la forme p(m,n)
;ou p(0,n) est le nombre d option (algebrique + pour achat - pour vente)
;   p(1,n) est le type de l'option  1 pour call , 2 pour put
;   p(2,n) est le prix d'exercice 
;   p(3,n) est l echeance  en jour julien : accessible par exemple via la fonction jj(jjmmaa)
;   p(4,n) est la volatilite implicite annuelle a utiliser pour valoriser le portefeuille
;   p(5,n) est le numero d ordre du support
;   p(6,n) est le multiplicateur de la classe d option (pour le cac c est 200)
;   p(7,n) est la liquidite maximale sur ce type d instrument (en nombre de contrats ou de titre)
;   p(8,n) est le taux de frais sur ce type d instrument 


;fonctions qui font passer d un format ou les quantites sont >0 pour des achats et <0 pour des ventes 
; a un format ou les quantites sont toujours positives mais il a une ligne pourt les achat et une ligne pour les ventes
; achat en 1 , vente en 2

function pexpand,p
psize=size(p)
n=psize(2)
pnew=dblarr(9,2*n)
for i=0,n-1 do if p(i) ge 0 then  pnew(0,2*i)=p(0,i)  else pnew(0,2*i+1)=-p(0,i)
for i=0,n-1 do begin
		pnew(1,2*i)=p(1,i) 
		pnew(1,2*i+1)=-p(1,i)
	endfor
for j=2,8 do  for i=0,n-1 do begin
		pnew(j,2*i)=p(j,i) 
		pnew(j,2*i+1)=p(j,i)
	endfor
return,pnew
end

function pcompress,p
psize=size(p)
n=psize(2)
pnew=dblarr(9,n/2)
for i=0,n-2,2 do  pnew(0,i/2)=p(0,i)-p(0,i+1)
for j=1,8 do for i=0,n-2,2 do  pnew(j,i/2)=p(j,i) 
return,pnew
end

function instrum_valuation,type,exercice,echeance,volatimplicite,s,deltavolannuel,dateev,r,qte
; echeance sous la forme jj(jjmmaa)
; volatimplicite sous la forme 0.39 pour 39 %
; s cours du support
;dateevaluation est exprimee en jj(jjmmaa)
;deltavolannuel vaut 0.05 pour un ecart de vol de 5%/an
nbjdatet=dateev
tauxj=(1.+r)^(1/365.)-1.
v=option(type,s,exercice,echeance-nbjdatet,tauxj,max([0,(volatimplicite+deltavolannuel)/(sqrt(365.)) ]))
return,v*qte
end

function instrum_gamma,type,exercice,echeance,volatimplicite,s,deltavolannuel,dateev,r,qte
; echeance sous la forme jj(jjmmaa)
; volatimplicite sous la forme 0.39 pour 39 %
; s cours du support
;dateevaluation est exprimee en jj(jjmmaa)
;deltavolannuel vaut 0.05 pour un ecart de vol de 5%/an
nbjdatet=dateev
tauxj=(1.+r)^(1/365.)-1.
v=option_gamma(type,s,exercice,echeance-nbjdatet,tauxj,max([0,(volatimplicite+deltavolannuel)/(sqrt(365.)) ]))
return,v*qte
end

function instrum_thetan,type,exercice,echeance,volatimplicite,s,deltavolannuel,dateev,r,qte,nbj
; echeance sous la forme jj(jjmmaa)
; volatimplicite sous la forme 0.39 pour 39 %
; s cours du support
;dateevaluation est exprimee en jj(jjmmaa)
;deltavolannuel vaut 0.05 pour un ecart de vol de 5%/an
nbjdatet=dateev
tauxj=(1.+r)^(1/365.)-1.
v1=option(type,s,exercice,echeance-nbjdatet,tauxj,max([0,(volatimplicite+deltavolannuel)/(sqrt(365.)) ]))
v2=option(type,s,exercice,echeance-(nbjdatet+nbj),tauxj,max([0,(volatimplicite+deltavolannuel)/(sqrt(365.)) ]))
return,(v2-v1)*qte
end

function instrum_traduc,type,exercice,echeance,sup
d=caldat(echeance)
;calcule le jour julien a partir du jour jjmmaa , inverse de caldat
a=d mod 100
m=((d-a) mod 10000)/100
if a lt 10 then a=a+100
if type eq 1 or type eq 2 then begin
case 1 of
m eq 1 : month='jan'
m eq 2 : month='fev'
m eq 3 : month='mar'
m eq 4 : month='avr'
m eq 5 : month='mai'
m eq 6 : month='jun'
m eq 7 : month='juil'
m eq 8 : month='aout'
m eq 9 : month='sep'
m eq 10 : month='oct'
m eq 11 : month='nov'
m eq 12 : month='dec'
endcase
endif
case 1 of
sup eq 1 : support='CAC40'
sup eq 2 : support='FTSE'
sup eq 3 : support='DAX'
sup eq 4 : support='SP500'
sup eq 5 : support='NIKKEI'
else : support='support inconnu'
endcase
case 1 of
type eq 1 : instrum='call'
type eq 2 : instrum='put '
type eq 3 : instrum='future ou cash         '
else : instrum='instrument inconnu'
endcase
case 1 of
type eq 1 or type eq 2 :res=instrum+' '+strtrim(string(exercice))+' '+month+strtrim(string(a,format='(i2)'))+' / '+support
type eq 3:res=instrum+' '+' / '+support
endcase
return,res
end


function pvalue,p,s,r,deltavolannuel,dateevaluation
;dateevaluation est exprimee en jj(jjmmaa)
;deltavolannuel vaut 0.05 pour un ecart de vol de 5%/an
nbjdatet=dateevaluation
tauxj=(1.+r)^(1/365.)-1.
psize=size(p)
n=psize(2)
v=0
for i=0,n-1 do begin
v=v+p(0,i)*p(6,i)*option(p(1,i),s,p(2,i),p(3,i)-nbjdatet,tauxj,max([0, (p(4,i)+deltavolannuel) / (sqrt(365.)) ]))
endfor
return,v
end

function pgamma,p,s,r,deltavolannuel,dateevaluation
;dateevaluation est exprimee en jj(jjmmaa)
;deltavolannuel vaut 0.05 pour un ecart de vol de 5%/an
nbjdatet=dateevaluation
tauxj=(1.+r)^(1/365.)-1.
psize=size(p)
n=psize(2)
v=0
for i=0,n-1 do begin
v=v+p(0,i)*p(6,i)*option_gamma(p(1,i),s,p(2,i),p(3,i)-nbjdatet,tauxj,max([0, (p(4,i)+deltavolannuel) / (sqrt(365.)) ]))
endfor
return,v
end

function pthetan,p,s,r,deltavolannuel,dateevaluation,nbjtheta
;dateevaluation est exprimee en jj(jjmmaa)
;deltavolannuel vaut 0.05 pour un ecart de vol de 5%/an
nbjdatet=dateevaluation
tauxj=(1.+r)^(1/365.)-1.
psize=size(p)
n=psize(2)
v1=0
v2=0
for i=0,n-1 do begin
v1=v1+p(0,i)*p(6,i)*option(p(1,i),s,p(2,i),p(3,i)-nbjdatet,tauxj,max([0, (p(4,i)+deltavolannuel) / (sqrt(365.)) ]))
v2=v2+p(0,i)*p(6,i)*option(p(1,i),s,p(2,i),p(3,i)-(nbjdatet+nbjtheta),tauxj,max([0, (p(4,i)+deltavolannuel) / (sqrt(365.)) ]))
endfor
return,v2-v1
end

function pchange,p,newq
psize=size(p)
n=psize(2)
pnew=dblarr(9,n)
for i=0,n-1 do pnew(0,i)=newq(i)
for j=1,8 do pnew(j,0:n-1)=p(j,0:n-1)
return,pnew
end

function pplus,p,newq
psize=size(p)
n=psize(2)
pnew=dblarr(9,n)
for i=0,n-1 do pnew(0,i)=newq(i)+p(0,i)
for j=1,8 do pnew(j,0:n-1)=p(j,0:n-1)
return,pnew
end

function pvbornesv,xsupport,yvol,p1,r,datecalcul,ss,vv,dateevaluation
;p : portefeuille calcule ses bornes de variation en performance
; yvol un vecteur de variation de volatilite annuel par rapport au vol de reference fournies dans p
;dateevaluation est de la forme " jj(010192) " soit en fait un nombre de jour "
nx=n_elements(xsupport)
ny=n_elements(yvol)
valo=dblarr(nx,ny)
p0=pvalue(p1,ss,r,vv,dateevaluation)
for i=0,nx-1 do for j=0,ny-1 do valo(i,j)=pvalue(p1,xsupport(i),r,yvol(j),datecalcul)-p0
;print,'valo risque=',valo
return,[min(valo),max(valo)]
end

function pgbornesv,xsupport,yvol,p1,r,datecalcul,ss,vv,dateevaluation
;p1 : portefeuille calcule ses bornes de variation en gamma
; yvol un vecteur de variation de volatilite annuel par rapport au vol de reference fournies dans p
;dateevaluation est de la forme " jj(010192) " soit en fait un nombre de jour "
nx=n_elements(xsupport)
ny=n_elements(yvol)
valo=dblarr(nx,ny)
for i=0,nx-1 do for j=0,ny-1 do valo(i,j)=pgamma(p1,xsupport(i),r,yvol(j),datecalcul)
return,[min(valo),max(valo)]
end

function ptbornesv,xsupport,yvol,p1,r,datecalcul,ss,vv,dateevaluation,n
;p : portefeuille calcule ses bornes de variation en theta
; yvol un vecteur de variation de volatilite annuel par rapport au vol de reference fournies dans p
;dateevaluation est de la forme " jj(010192) " soit en fait un nombre de jour "
nx=n_elements(xsupport)
ny=n_elements(yvol)
valo=dblarr(nx,ny)
for i=0,nx-1 do for j=0,ny-1 do valo(i,j)=pthetan(p1,xsupport(i),r,yvol(j),datecalcul,n)
return,[min(valo),max(valo)]
end

pro psurface,xsupport,yvol,p1,r,datecalcul,ss,vv,dateevaluation
;p : portefeuille
; yvol un vecteur de variation de volatilite annuel par rapport au vol de reference fournies dans p
;dateevaluation est de la forme " jj(010192) " soit en fait un nombre de jour "
nx=n_elements(xsupport)
ny=n_elements(yvol)
valo=dblarr(nx,ny)
p0=pvalue(p1,ss,r,vv,dateevaluation)
for i=0,nx-1 do for j=0,ny-1 do valo(i,j)=pvalue(p1,xsupport(i),r,yvol(j),datecalcul)-p0
surface,valo,xsupport,yvol
end


pro p_print,p
psize=size(p)
n=psize(2)
print,'-------------------------Portefeuille-----------------------------'
for i=0,n-1 do print,string(p(0,i),format='(i6)'),' ',instrum_traduc(int(p(1,i)),p(2,i),p(3,i),p(5,i))
print,'------------------------------------------------------------------'
end  



pro p_print2,p,pp
psize=size(p)
n=psize(2)
print,'------------------------------------initial----ordres-----final---'
for i=0,n-1 do print,instrum_traduc(int(p(1,i)),p(2,i),p(3,i),p(5,i))$
		,' ',string(p(0,i),format='(i10)'),string(pp(0,i)-p(0,i),format='(i10)'),string(pp(0,i),format='(i10)')
print,'------------------------------------------------------------------'
end  

pro p_print2f,p,pp
psize=size(p)
n=psize(2)
print,'--------------------------------------initial----ordres optimaux------ordres effectifs---------final---'
for i=0,n-1 do print,instrum_traduc(int(p(1,i)),p(2,i),p(3,i),p(5,i))$
  	 ,' ',string(p(0,i),format='(i10)'),string(pp(0,i)-p(0,i),format='(f20.10)'),$
	'     || ',string(pp(0,i)-p(0,i),format='(i10)'),'    || ',string(pp(0,i),format='(i10)')
print,'-------------------------------------------------------------------------------------------------------'
end  

function vector_one,x,p,n
;retourne un vecteur [0,0,..,0,x,0,..,0] de n composantes
; ou le x est a la position p (0<=p<=n-1)
a=replicate(0,n)
a(p)=x
return,a
end		


pro ecrit_equation,ai,b0,sgn,label
common cpoptimisation2,a,b,c,ni,n_eq,n_var,eq_type,eq_sgn
;n_var : pointeur courrant indiquant le numero de la prochaine variable a creer
;n_eq : pointeur courrant indiquant le numero de la prochaine equation a ecrire
;sous ls forme somme(ai * xi) >= a0 si sgn=+1
;sous ls forme somme(ai * xi) = a0 si sgn=0
;sous ls forme somme(ai * xi) <= a0 si sgn=-1
eq_type(n_eq)=sgn
eq_sgn(n_eq)=sgn(b0)
case 1 of
	sgn eq 0 and b0 ge 0:begin
			a(n_eq,0:ni-1)=ai
			b(n_eq)=b0
		 endif
	sgn eq 0 and b0 lt 0:begin
			a(n_eq,0:ni-1)=-ai
			b(n_eq)=-b0
		 endif
	sgn eq 1 and b0 ge 0 :begin
			a(n_eq,0:ni-1)=ai
			b(n_eq)=b0
			a(n_eq,n_var)=-1.
			n_var=n_var+1
		 endif
	sgn eq 1 and b0 lt 0:begin
			a(n_eq,0:ni-1)=-ai
			b(n_eq)=-b0
			a(n_eq,n_var)=1.
			n_var=n_var+1
		 endif
	sgn eq -1 and b0 ge 0:begin
			a(n_eq,0:ni-1)=ai
			b(n_eq)=b0
			a(n_eq,n_var)=1.
			n_var=n_var+1
		 endif
	sgn eq -1and b0 lt 0:begin
			a(n_eq,0:ni-1)=-ai
			b(n_eq)=-b0
			a(n_eq,n_var)=-1.
			n_var=n_var+1
		 endif
else:print,'arguments non interpretable pour ecrit_equation:',sgn,b0
endcase
;print,'n_eq=',n_eq,label
;print,'ai=',transpose(a(n_eq,*))
;print,'b=',b(n_eq)
n_eq=n_eq+1
end


pro ecrit_equation2,ai,b0,b02,sgn,label
common cpoptimisation2,a,b,a2,b2,c,ni,n_eq,n_var,eq_type,eq_sgn
;n_var : pointeur courrant indiquant le numero de la prochaine variable a creer
;n_eq : pointeur courrant indiquant le numero de la prochaine equation a ecrire
;sous ls forme somme(ai * xi) >= a0 si sgn=+1
;sous ls forme somme(ai * xi) = a0 si sgn=0
;sous ls forme somme(ai * xi) <= a0 si sgn=-1
eq_type(n_eq)=sgn
eq_sgn(n_eq)=sgn(b0)
case 1 of
	sgn eq 0 and b0 ge 0:begin
			a(n_eq,0:ni-1)=ai
			b(n_eq)=b0
		 endif
	sgn eq 0 and b0 lt 0:begin
			a(n_eq,0:ni-1)=-ai
			b(n_eq)=-b0
		 endif
	sgn eq 1 and b0 ge 0 :begin
			a(n_eq,0:ni-1)=ai
			b(n_eq)=b0
			a(n_eq,n_var)=-1.
			n_var=n_var+1
		 endif
	sgn eq 1 and b0 lt 0:begin
			a(n_eq,0:ni-1)=-ai
			b(n_eq)=-b0
			a(n_eq,n_var)=1.
			n_var=n_var+1
		 endif
	sgn eq -1 and b0 ge 0:begin
			a(n_eq,0:ni-1)=ai
			b(n_eq)=b0
			a(n_eq,n_var)=1.
			n_var=n_var+1
		 endif
	sgn eq -1and b0 lt 0:begin
			a(n_eq,0:ni-1)=-ai
			b(n_eq)=-b0
			a(n_eq,n_var)=-1.
			n_var=n_var+1
		 endif
else:print,'arguments non interpretable pour ecrit_equation:',sgn,b0
endcase
case 1 of
	sgn eq 0 and b02 ge 0:begin
			a2(n_eq,0:ni-1)=ai
			b2(n_eq)=b02
		 endif
	sgn eq 0 and b02 lt 0:begin
			a2(n_eq,0:ni-1)=-ai
			b2(n_eq)=-b02
		 endif
	sgn eq 1 and b02 ge 0 :begin
			a2(n_eq,0:ni-1)=ai
			b2(n_eq)=b02
			a2(n_eq,n_var)=-1.
			n_var=n_var+1
		 endif
	sgn eq 1 and b02 lt 0:begin
			a2(n_eq,0:ni-1)=-ai
			b2(n_eq)=-b02
			a2(n_eq,n_var)=1.
			n_var=n_var+1
		 endif
	sgn eq -1 and b02 ge 0:begin
			a2(n_eq,0:ni-1)=ai
			b2(n_eq)=b02
			a2(n_eq,n_var)=1.
			n_var=n_var+1
		 endif
	sgn eq -1and b02 lt 0:begin
			a2(n_eq,0:ni-1)=-ai
			b2(n_eq)=-b02
			a2(n_eq,n_var)=-1.
			n_var=n_var+1
		 endif
else:print,'arguments non interpretable pour ecrit_equation:',sgn,b02
endcase

;print,'n_eq=',n_eq,label
;print,'ai=',transpose(a(n_eq,*))
;print,'b=',b(n_eq)
n_eq=n_eq+1
end



function poptimisation2,pinitial,pinitial2,s0,deltav0,r,db,dr,nvolr,nsr,nvolb,nsb,deltajb,deltajr,dateevaluation,sumrisk,impflag,$
opt=opt,pond=pond,nbjtheta=nbjtheta,psup=psup,pinf=pinf,gsup=gsup,ginf=ginf,tsup=tsup,tinf=tinf,$
pmsup=pmsup,pminf=pminf,gmsup=gmsup,gminf=gminf,tmsup=tmsup,nbmax=nbmax,fraismax=fraismax
common cpoptimisation2,a,b,a2,b2,c,ni,n_eq,n_var,eq_type,eq_sgn
;fonction qui uilisant un portefeuille p avec les limite de liquidite et deux domaine DR et DB ,
; une date JJMMAA ou le profit est regarde sur le dommaine DB a la date : dateprofit=dateevaluation + deltajr
; une date JJMMAA ou le risque est calcule sur le domaine DR a la date : daterisque=dateevaluation + deltajb
; les domaines sont repere par une liste de la forme [[supmin,supmax],[volmin,volmax]]
; aussi bien pour le domaine DB que pour le domaine DR
;rends une liste d'ordres c'est a dire un autre portefeuille (difference entre deux portefeuilles)
;qui optimise le portefeuille precedent
;pinitial portefeuille initial et en meme temps liste des instruments sur lequel on realise l optimisation
;s0 prix central du support
;deltav0 ecart en point de volatilite implicite observe sur le marche par rapport au vol donne
;dans pinitial
;r taux d'interet
;db liste de la forme ([[sminb,smaxb],[volminb,volmaxb]])
;dr liste de la forme ([[sminr,smaxr],[volminr,volmaxr]])
;fixe les domaine db surlequel on mqximise le benfice et
;dr,ds domaine sur lequel on minimise le risque (b ou r)
;nvolr,nvols nombre de pas de volatilite (b ou r)
;ns nombre de pas de support
;deltajb nomde jours forward donnant la date de prise en compte du benefice
;deltajr nomde jours forward donnant la date de prise en compte du risque
;dateevaluation date a laquelle on effectue l'evaluation au format jj(jjmmaa)
;sumrisk vaut 1e6 pour un risque possible maximum sur le domaine de risque dr de 1 million
;opt='maxprofit' pour maximiser le profit (valeur par defaut)
; mais on peut mettre  minprofit,maxgamma ,mingamma ,maxtheta ou mintheta
;exemple :
; pour pinitial :p=[ 	[100	 ,1 ,1800 ,jj(310692) ,30.9, 1,200,1000,0.01],$ pour un call
;			[100	 ,2 ,2200 ,jj(310692) ,30.9 1,200,,500,0.01],$  pour un put
;			[100	 ,3 ,0    ,jj(310692) ,30.9 1,200,,500,0.0003]]  pour un futur
; pour l appel a la fonction :
;potimisation(p,2021.,1.,0.097,[[2150.,2200.],[-5.,5.]],[[1800.,2200.],
;        [-10.,10.]],3,3,20,1,jj(010192),1e6)
;nombre d'instrument : ni
;generation des ensembles de points  DB et DR de la forme [[s1,v1],..., [sn,vn]]
t0=systime(1)
if not keyword_set(nbjtheta) then begin
	print,': oubli de mentionner NBJTHETA=  il est pris =1 par defaut'
	nbjtheta=1
endif
volminb=db(0,1) & volmaxb=db(1,1) & sminb=db(0,0) & smaxb=db(1,0)
volminr=dr(0,1) & volmaxr=dr(1,1) & sminr=dr(0,0) & smaxr=dr(1,0)
nip=size(pinitial)
ni=nip(2)
liquidite=dblarr(ni)
nbtb=nvolb*nsb
nbtr=nvolr*nsr
dr=dblarr(nbtr,2)
db=dblarr(nbtb,2)
rjour=(1.+r)^(1/365.)-1.
i=0
for jv=0,nvolr-1 do for js=0,nsr-1 do begin
	dr(i,0)=(smaxr-sminr)*js/float(nsr-1)+sminr
	dr(i,1)=(volmaxr-volminr)*jv/float(nvolr-1)+volminr
	i=i+1
endfor
i=0
for jv=0,nvolb-1 do for js=0,nsb-1 do begin
	db(i,0)=(smaxb-sminb)*js/float(nsb-1)+sminb
	db(i,1)=(volmaxb-volminb)*jv/float(nvolb-1)+volminb
	i=i+1
endfor
volmoyenne=total(pinitial(4,*))/ni
;fabrication des valorisations des instruments 
vi=dblarr(ni)
for ji=0,ni-1 do begin 
		liquidite(ji)=pinitial(7,ji)
		vi(ji)=instrum_valuation(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		s0,deltav0,dateevaluation,r,pinitial(6,ji))
		endfor
;fabrication des vecteurs de pertes locales (profits sur les points de risque)
vr=dblarr(nbtr,ni)
for jp=0,nbtr-1 do for ji=0,ni-1 do vr(jp,ji)=instrum_valuation(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		dr(jp,0),dr(jp,1),dateevaluation + deltajr,r,pinitial(6,ji))- vi(ji)
if keyword_set(opt) then opt=opt else opt='maxprofit'
if keyword_set(pond) then pond=pond else pond=0
if (opt eq 'maxprofit') or (opt eq 'minprofit') or keyword_set(psup) or keyword_set(pmsup) or keyword_set(pinf) or keyword_set(pminf) then begin
;fabrication des vecteurs de  profits locaux  (profits sur les points de benefice)
vb=dblarr(nbtb,ni)
for jp=0,nbtb-1 do for ji=0,ni-1 do vb(jp,ji)=instrum_valuation(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		db(jp,0),db(jp,1),dateevaluation + deltajb,r,pinitial(6,ji))-vi(ji)
endif
if (opt eq 'maxgamma') or (opt eq 'mingamma') or keyword_set(gsup) or keyword_set(gmsup) or keyword_set(ginf) or keyword_set(gminf) then begin
;fabrication des vecteur de gamma locaux
gvb=dblarr(nbtb,ni)
for jp=0,nbtb-1 do for ji=0,ni-1 do gvb(jp,ji)=instrum_gamma(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		db(jp,0),db(jp,1),dateevaluation + deltajb,r,pinitial(6,ji))
endif
if (opt eq 'maxtheta') or (opt eq 'mintheta') or keyword_set(tsup) or keyword_set(tmsup) or keyword_set(tinf) or keyword_set(tminf) then begin
;fabrication des vecteur de theta locaux
tvb=dblarr(nbtb,ni)
for jp=0,nbtb-1 do for ji=0,ni-1 do tvb(jp,ji)=instrum_thetan(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		db(jp,0),db(jp,1),dateevaluation + deltajb,r,pinitial(6,ji),nbjtheta)
endif
; optimisation 
case opt of
	'maxprofit': hvb=vb
	'minprofit': hvb=vb
	'maxgamma': hvb=gvb
	'mingamma': hvb=gvb
	'maxtheta': hvb=tvb
	'mintheta': hvb=tvb
else:print ,'opt a une valeur interdite :',opt
endcase
	neqtotal=ni+nbtr
	if keyword_set(pinf) then neqtotal=neqtotal+nbtb
	if keyword_set(ginf) then neqtotal=neqtotal+nbtb
	if keyword_set(tinf) then neqtotal=neqtotal+nbtb
	if keyword_set(psup) then neqtotal=neqtotal+nbtb
	if keyword_set(gsup) then neqtotal=neqtotal+nbtb
	if keyword_set(tsup) then neqtotal=neqtotal+nbtb
	if keyword_set(pminf) then neqtotal=neqtotal+1
	if keyword_set(gminf) then neqtotal=neqtotal+1
	if keyword_set(tminf) then neqtotal=neqtotal+1
	if keyword_set(pmsup) then neqtotal=neqtotal+1
	if keyword_set(gmsup) then neqtotal=neqtotal+1
	if keyword_set(tmsup) then neqtotal=neqtotal+1
	if keyword_set(nbmax) then neqtotal=neqtotal+1
	if keyword_set(fraismax) then neqtotal=neqtotal+1
	m=neqtotal
	n=neqtotal+ni
	a=dblarr(m,n)
	b=dblarr(m)
	b2=dblarr(m)
	eq_type=dblarr(m)
	eq_sgn=dblarr(m)
	c=dblarr(n)
	case pond of
	0:for ji=0,ni-1 do c(ji)=-total(hvb(*,ji))/nbtb
	1:begin
		pond=dblarr(nbtb)
		for jbi =0,nbtb-1 do begin
 		pond(jbi)=exp(-(db(jbi,0)/s0-exp(-rjour*deltajb))^2/(2*((volmoyenne+db(jbi,1))^2))) 
		endfor
		totalpond=total(pond)
		for ji=0,ni-1 do c(ji)=total(vb(*,ji)*pond)/totalpond
	 end
	else:print ,'pond a une valeur interdite :',pond
	endcase
	tab_aux=dblarr(ni)
 	n_eq=0
	n_var=ni
;---------------------------------------------------------------------------------------
	for jbi=0,nbtr-1 do  begin
		label='risque(s='+string(dr(jbi,0))+',v='+string(dr(jbi,1))+')'		;profit >= -sumrisk dans domaine de risque
		b0=total(vr(jbi,0:ni-1)*pinitial(0,0:ni-1))
		b02=total(vr(jbi,0:ni-1)*pinitial2(0,0:ni-1))
		ecrit_equation,vr(jbi,0:ni-1),-(sumrisk+b0),-(sumrisk+b02),1,label
		endfor									
	for jbi=0,ni-1 do ecrit_equation2,vector_one(1,jbi,ni),liquidite(jbi),-1,$	;instruments traites < liquidite dans chaque classe
		'liquidite de i='+string(jbi)
	label=''					
	if keyword_set(pinf) then begin							;profit > pinf dans domaine d esperance
		print,'profit >=', pinf ,' dans domaine d esperance'
		for jbi=0,nbtb-1 do begin
		b0=total(vb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,vb(jbi,0:ni-1),pinf-b0,1,label
		endfor
	endif
	if keyword_set(psup)  then begin      			                        ;profit < psup dans domaine d esperance
		print,'profit <= ',psup,'  dans domaine d esperance'
		for jbi=0,nbtb-1 do begin
		b0=total(vb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,vb(jbi,0:ni-1),psup-b0,-1,label
		endfor
	endif
	if keyword_set(ginf) then begin         					;gamma > ginf dans domaine d esperance
		print,'gamma >= ',ginf,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do  begin
		b0=total(gvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
		b02=total(gvb(jbi,0:ni-1)*pinitial2(0,0:ni-1))
			ecrit_equation2,gvb(jbi,0:ni-1),ginf-b0,ginf-b02,1,label
		endfor
	endif
	if keyword_set(gsup) then begin         					;gamma < gsup dans domaine d esperance
		print,'gamma <=',gsup,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do begin
		b0=total(gvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,gvb(jbi,0:ni-1),gsup-b0,-1,label
		endfor
	endif
	if keyword_set(tinf) then begin							;theta > tinf dans domaine d esperance
		print,'theta >= ',tinf,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do begin
		b0=total(tvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,tvb(jbi,0:ni-1),tinf-b0,1,label
		endfor
	endif
	if keyword_set(tsup) then begin							;theta < tsup dans domaine d esperance
		print,'theta <= ',tsup,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do begin
		b0=total(tvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,tvb(jbi,0:ni-1),tsup-b0,-1,label
		endfor
	endif
	if keyword_set(pminf) then begin                 				;profit moyen > pminf dans domaine d esperance
		print,'profit moyen >= ',pminf ,' dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(vb(*,ji))/nbtb
		ecrit_equation,tab_aux,pminf,1,label
	endif
	if keyword_set(pmsup) then begin						;profit moyen <pmsup dans domaine d esperance
		print,'profit moyen <= ',pmsup,'  dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(vb(*,ji))/nbtb
		ecrit_equation,tab_aux,pmsup,-1,label
	endif
	if keyword_set(gminf) then begin						;gamma moyen > gminf dans domaine d esperance
		print,'gamma moyen >= ',gminf,'  dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(gvb(*,ji))/nbtb
		ecrit_equation,tab_aux,gminf,1,label
	endif
	if keyword_set(gmsup) then begin						;gamma moyen <gmsup dans domaine d esperance
		print,'gamma moyen <= ',gmsup ,' dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(gvb(*,ji))/nbtb
		ecrit_equation,tab_aux,gmsup,-1,label
	endif
	if keyword_set(tminf) then begin						;theta moyen > tminf dans domaine d esperance
		print,'theta moyen >= ',tminf,' dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(tvb(*,ji))/nbtb
		ecrit_equation,tab_aux,tminf,1,label
	endif
	if keyword_set(tmsup) then begin						;theta moyen <tmsup dans domaine d esperance
		print,'theta moyen <= ',tmsup,'  dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(tvb(*,ji))/nbtb
		ecrit_equation,tab_aux,tmsup,-1,label
	endif
	if keyword_set(nbmax) then begin                                                ;total des instruments traites < nbmax
		print,'total des instruments traites <=', nbmax
		ecrit_equation,replicate(1,ni),nbmax,-1
	endif
	if keyword_set(fraismax) then begin                                             ;total des frais < fraismax
		print,'total des frais <= ',fraismax
		ecrit_equation,abs(vi(0:ni-1))*pinitial(8,0:ni-1),fraismax,-1
	endif
;------------------------------------------------------------------------------
t1=systime(1)
print,neqtotal,' inequations / ',ni+neqtotal,' variables'
print,'optimisation...'
toexec='result=bisimplexe(a,b,b2,c)'
flagexe=execute(toexec)
t2=systime(1)
if flagexe ne 1 then begin
	print,'erreur durant l optimisation'
	return,[[0,0]]
endif else begin
	if result(0,0) eq 1 then begin
		print,'erreur dans le contraintes : solution non bornee'
		return,0
		endif
	if result(0,0) eq -1 then begin
		print,'erreur : pas de solution possible aux contraintes"
		mm=result(1,0)
		print,'les contraintes suivantes semblent incriminees:',transpose(result(1,1:mm))
		return,0
		endif
	if result(0,0) eq 0 then begin
		vector=dblarr(n)
		vector(0:n-1)=result(4,1:n)
		eb=dblarr(ni)
		for ji=0,ni-1 do eb(ji)=total(vb(*,ji))/nbtb
		pvaleurmoyb=total(vector*eb)
		pnouveau=pplus(pinitial,vector)
		xsupport=dblarr(nsr)
		yvol=dblarr(nvolr)
		for js=0,nsr-1 do xsupport(js)=(smaxr-sminr)*js/float(nsr-1)+sminr
		for jv=0,nvolr-1 do yvol(jv)=(volmaxr-volminr)*jv/float(nvolr-1)+volminr
		;window,0
		;psurface,xsupport,yvol,pnouveau,r,dateevaluation+deltajr,s0,deltav0,dateevaluation
		print,'temps prep. data =',t1-t0,' temps opt. =',t2-t1
		;print,'valeur moyenne du dommaine d esperance',pvaleurmoyb
		bornes=pvbornesv(xsupport,yvol,pnouveau,r,dateevaluation+deltajr,s0,deltav0,dateevaluation)
		;print,'bornes en profit du domaine de risque',bornes
		bxsupport=dblarr(nsb)
		byvol=dblarr(nvolb)
		for js=0,nsb-1 do bxsupport(js)=(smaxb-sminb)*js/float(nsb-1)+sminb
		for jv=0,nvolb-1 do byvol(jv)=(volmaxb-volminb)*jv/float(nvolb-1)+volminb
		bornes=pvbornesv(bxsupport,byvol,pnouveau,r,dateevaluation+deltajb,s0,deltav0,dateevaluation)
		;print,'bornes en profit du domaine d esperance',bornes
		gbornes=pgbornesv(bxsupport,byvol,pnouveau,r,dateevaluation+deltajb,s0,deltav0,dateevaluation)
		;print,'bornes en gamma du domaine d esperance',gbornes
		tbornes=ptbornesv(bxsupport,byvol,pnouveau,r,dateevaluation+deltajb,s0,deltav0,dateevaluation,nbjtheta)
		;print,'bornes en theta du domaine d esperance',tbornes
		frais0=0
		for i=0,ni-2,2 do  begin
			frais0=frais0+abs(vector(i)-vector(i+1))*abs(vi(i))*pinitial(8,i)
		endfor
		;print,'montant des frais associe au changement de position:',frais0
		return,pnouveau
	endif
endelse
end


pro testequations,pordres,pinitial,s0,deltav0,r,db,dr,nvolr,nsr,nvolb,nsb,deltajb,deltajr,dateevaluation,sumrisk,impflag,$
opt=opt,pond=pond,nbjtheta=nbjtheta,psup=psup,pinf=pinf,gsup=gsup,ginf=ginf,tsup=tsup,tinf=tinf,$
pmsup=pmsup,pminf=pminf,gmsup=gmsup,gminf=gminf,tmsup=tmsup,nbmax=nbmax,fraismax=fraismax
common cpoptimisation2,a,b,a2,b2,c,ni,n_eq,n_var,eq_type,eq_sgn
;fonction qui uilisant un portefeuille p avec les limite de liquidite et deux domaine DR et DB ,
; une date JJMMAA ou le profit est regarde sur le dommaine DB a la date : dateprofit=dateevaluation + deltajr
; une date JJMMAA ou le risque est calcule sur le domaine DR a la date : daterisque=dateevaluation + deltajb
; les domaines sont repere par une liste de la forme [[supmin,supmax],[volmin,volmax]]
; aussi bien pour le domaine DB que pour le domaine DR
;rends une liste d'ordres c'est a dire un autre portefeuille (difference entre deux portefeuilles)
;qui optimise le portefeuille precedent
;pinitial portefeuille initial et en meme temps liste des instruments sur lequel on realise l optimisation
;s0 prix central du support
;deltav0 ecart en point de volatilite implicite observe sur le marche par rapport au vol donne
;dans pinitial
;r taux d'interet
;db liste de la forme ([[sminb,smaxb],[volminb,volmaxb]])
;dr liste de la forme ([[sminr,smaxr],[volminr,volmaxr]])
;fixe les domaine db surlequel on mqximise le benfice et
;dr,ds domaine sur lequel on minimise le risque (b ou r)
;nvolr,nvols nombre de pas de volatilite (b ou r)
;ns nombre de pas de support
;deltajb nomde jours forward donnant la date de prise en compte du benefice
;deltajr nomde jours forward donnant la date de prise en compte du risque
;dateevaluation date a laquelle on effectue l'evaluation au format jj(jjmmaa)
;sumrisk vaut 1e6 pour un risque possible maximum sur le domaine de risque dr de 1 million
;opt='maxprofit' pour maximiser le profit (valeur par defaut)
; mais on peut mettre  minprofit,maxgamma ,mingamma ,maxtheta ou mintheta
;exemple :
; pour pinitial :p=[ 	[100	 ,1 ,1800 ,jj(310692) ,30.9, 1,200,1000,0.01],$ pour un call
;			[100	 ,2 ,2200 ,jj(310692) ,30.9 1,200,,500,0.01],$  pour un put
;			[100	 ,3 ,0    ,jj(310692) ,30.9 1,200,,500,0.0003]]  pour un futur
; pour l appel a la fonction :
;potimisation(p,2021.,1.,0.097,[[2150.,2200.],[-5.,5.]],[[1800.,2200.],
;        [-10.,10.]],3,3,20,1,jj(010192),1e6)
;nombre d'instrument : ni
;generation des ensembles de points  DB et DR de la forme [[s1,v1],..., [sn,vn]]
epsilon=1.e-5
t0=systime(1)
if not keyword_set(nbjtheta) then begin
	print,': oubli de mentionner NBJTHETA=  il est pris =1 par defaut'
	nbjtheta=1
endif
volminb=db(0,1) & volmaxb=db(1,1) & sminb=db(0,0) & smaxb=db(1,0)
volminr=dr(0,1) & volmaxr=dr(1,1) & sminr=dr(0,0) & smaxr=dr(1,0)
nip=size(pinitial)
ni=nip(2)
liquidite=dblarr(ni)
nbtb=nvolb*nsb
nbtr=nvolr*nsr
dr=dblarr(nbtr,2)
db=dblarr(nbtb,2)
rjour=(1.+r)^(1/365.)-1.
i=0
for jv=0,nvolr-1 do for js=0,nsr-1 do begin
	dr(i,0)=(smaxr-sminr)*js/float(nsr-1)+sminr
	dr(i,1)=(volmaxr-volminr)*jv/float(nvolr-1)+volminr
	i=i+1
endfor
i=0
for jv=0,nvolb-1 do for js=0,nsb-1 do begin
	db(i,0)=(smaxb-sminb)*js/float(nsb-1)+sminb
	db(i,1)=(volmaxb-volminb)*jv/float(nvolb-1)+volminb
	i=i+1
endfor
volmoyenne=total(pinitial(4,*))/ni
;fabrication des valorisations des instruments 
vi=dblarr(ni)
for ji=0,ni-1 do begin 
		liquidite(ji)=pinitial(7,ji)
		vi(ji)=instrum_valuation(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		s0,deltav0,dateevaluation,r,pinitial(6,ji))
		endfor
;fabrication des vecteurs de pertes locales (profits sur les points de risque)
vr=dblarr(nbtr,ni)
for jp=0,nbtr-1 do for ji=0,ni-1 do vr(jp,ji)=instrum_valuation(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		dr(jp,0),dr(jp,1),dateevaluation + deltajr,r,pinitial(6,ji))- vi(ji)
if keyword_set(opt) then opt=opt else opt='maxprofit'
if keyword_set(pond) then pond=pond else pond=0
if (opt eq 'maxprofit') or (opt eq 'minprofit') or keyword_set(psup) or keyword_set(pmsup) or keyword_set(pinf) or keyword_set(pminf) then begin
;fabrication des vecteurs de  profits locaux  (profits sur les points de benefice)
vb=dblarr(nbtb,ni)
for jp=0,nbtb-1 do for ji=0,ni-1 do vb(jp,ji)=instrum_valuation(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		db(jp,0),db(jp,1),dateevaluation + deltajb,r,pinitial(6,ji))-vi(ji)
endif
if (opt eq 'maxgamma') or (opt eq 'mingamma') or keyword_set(gsup) or keyword_set(gmsup) or keyword_set(ginf) or keyword_set(gminf) then begin
;fabrication des vecteur de gamma locaux
gvb=dblarr(nbtb,ni)
for jp=0,nbtb-1 do for ji=0,ni-1 do gvb(jp,ji)=instrum_gamma(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		db(jp,0),db(jp,1),dateevaluation + deltajb,r,pinitial(6,ji))
endif
if (opt eq 'maxtheta') or (opt eq 'mintheta') or keyword_set(tsup) or keyword_set(tmsup) or keyword_set(tinf) or keyword_set(tminf) then begin
;fabrication des vecteur de theta locaux
tvb=dblarr(nbtb,ni)
for jp=0,nbtb-1 do for ji=0,ni-1 do tvb(jp,ji)=instrum_thetan(pinitial(1,ji),pinitial(2,ji),pinitial(3,ji),pinitial(4,ji), $
		db(jp,0),db(jp,1),dateevaluation + deltajb,r,pinitial(6,ji),nbjtheta)
endif
; optimisation 
case opt of
	'maxprofit': hvb=vb
	'minprofit': hvb=vb
	'maxgamma': hvb=gvb
	'mingamma': hvb=gvb
	'maxtheta': hvb=tvb
	'mintheta': hvb=tvb
else:print ,'opt a une valeur interdite :',opt
endcase
	neqtotal=ni+nbtr
	if keyword_set(pinf) then neqtotal=neqtotal+nbtb
	if keyword_set(ginf) then neqtotal=neqtotal+nbtb
	if keyword_set(tinf) then neqtotal=neqtotal+nbtb
	if keyword_set(psup) then neqtotal=neqtotal+nbtb
	if keyword_set(gsup) then neqtotal=neqtotal+nbtb
	if keyword_set(tsup) then neqtotal=neqtotal+nbtb
	if keyword_set(pminf) then neqtotal=neqtotal+1
	if keyword_set(gminf) then neqtotal=neqtotal+1
	if keyword_set(tminf) then neqtotal=neqtotal+1
	if keyword_set(pmsup) then neqtotal=neqtotal+1
	if keyword_set(gmsup) then neqtotal=neqtotal+1
	if keyword_set(tmsup) then neqtotal=neqtotal+1
	if keyword_set(nbmax) then neqtotal=neqtotal+1
	if keyword_set(fraismax) then neqtotal=neqtotal+1
	m=neqtotal
	n=neqtotal+ni
	a=dblarr(m,n)
	b=dblarr(m)
	eq_type=dblarr(m)
	eq_sgn=dblarr(m)
	c=dblarr(n)
	case pond of
	0:for ji=0,ni-1 do c(ji)=-total(hvb(*,ji))/nbtb
	1:begin
		pond=dblarr(nbtb)
		for jbi =0,nbtb-1 do begin
 		pond(jbi)=exp(-(db(jbi,0)/s0-exp(-rjour*deltajb))^2/(2*((volmoyenne+db(jbi,1))^2))) 
		endfor
		totalpond=total(pond)
		for ji=0,ni-1 do c(ji)=total(vb(*,ji)*pond)/totalpond
	 end
	else:print ,'pond a une valeur interdite :',pond
	endcase
	tab_aux=dblarr(ni)
 	n_eq=0
	n_var=ni
	for jbi=0,nbtr-1 do  begin
		label='risque(s='+string(dr(jbi,0))+',v='+string(dr(jbi,1))+')'		;profit >= -sumrisk dans domaine de risque
		b0=total(vr(jbi,0:ni-1)*pinitial(0,0:ni-1))
		ecrit_equation,vr(jbi,0:ni-1),-(sumrisk+b0),1,label
		endfor									
	for jbi=0,ni-1 do ecrit_equation,vector_one(1,jbi,ni),liquidite(jbi),-1,$	;instruments traites < liquidite dans chaque classe
		'liquidite de i='+string(jbi)
	label=''					
	if keyword_set(pinf) then begin							;profit > pinf dans domaine d esperance
		print,'profit >=', pinf ,' dans domaine d esperance'
		for jbi=0,nbtb-1 do begin
		b0=total(vb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,vb(jbi,0:ni-1),pinf-b0,1,label
		endfor
	endif
	if keyword_set(psup)  then begin      			                        ;profit < psup dans domaine d esperance
		print,'profit <= ',psup,'  dans domaine d esperance'
		for jbi=0,nbtb-1 do begin
		b0=total(vb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,vb(jbi,0:ni-1),psup-b0,-1,label
		endfor
	endif
	if keyword_set(ginf) then begin         					;gamma > ginf dans domaine d esperance
		print,'gamma >= ',ginf,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do  begin
		b0=total(gvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,gvb(jbi,0:ni-1),ginf-b0,1,label
		endfor
	endif
	if keyword_set(gsup) then begin         					;gamma < gsup dans domaine d esperance
		print,'gamma <=',gsup,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do begin
		b0=total(gvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,gvb(jbi,0:ni-1),gsup-b0,-1,label
		endfor
	endif
	if keyword_set(tinf) then begin							;theta > tinf dans domaine d esperance
		print,'theta >= ',tinf,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do begin
		b0=total(tvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,tvb(jbi,0:ni-1),tinf-b0,1,label
		endfor
	endif
	if keyword_set(tsup) then begin							;theta < tsup dans domaine d esperance
		print,'theta <= ',tsup,'  dans domaine d esperance '
		for jbi=0,nbtb-1 do begin
		b0=total(tvb(jbi,0:ni-1)*pinitial(0,0:ni-1))
			ecrit_equation,tvb(jbi,0:ni-1),tsup-b0,-1,label
		endfor
	endif
	if keyword_set(pminf) then begin                 				;profit moyen > pminf dans domaine d esperance
		print,'profit moyen >= ',pminf ,' dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(vb(*,ji))/nbtb
		ecrit_equation,tab_aux,pminf,1,label
	endif
	if keyword_set(pmsup) then begin						;profit moyen <pmsup dans domaine d esperance
		print,'profit moyen <= ',pmsup,'  dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(vb(*,ji))/nbtb
		ecrit_equation,tab_aux,pmsup,-1,label
	endif
	if keyword_set(gminf) then begin						;gamma moyen > gminf dans domaine d esperance
		print,'gamma moyen >= ',gminf,'  dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(gvb(*,ji))/nbtb
		ecrit_equation,tab_aux,gminf,1,label
	endif
	if keyword_set(gmsup) then begin						;gamma moyen <gmsup dans domaine d esperance
		print,'gamma moyen <= ',gmsup ,' dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(gvb(*,ji))/nbtb
		ecrit_equation,tab_aux,gmsup,-1,label
	endif
	if keyword_set(tminf) then begin						;theta moyen > tminf dans domaine d esperance
		print,'theta moyen >= ',tminf,' dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(tvb(*,ji))/nbtb
		ecrit_equation,tab_aux,tminf,1,label
	endif
	if keyword_set(tmsup) then begin						;theta moyen <tmsup dans domaine d esperance
		print,'theta moyen <= ',tmsup,'  dans domaine d esperance'
		for ji=0,ni-1 do tab_aux(ji)=total(tvb(*,ji))/nbtb
		ecrit_equation,tab_aux,tmsup,-1,label
	endif
	if keyword_set(nbmax) then begin                                                ;total des instruments traites < nbmax
		print,'total des instruments traites <=', nbmax
		ecrit_equation,replicate(1,ni),nbmax,-1
	endif
	if keyword_set(fraismax) then begin                                             ;total des frais < fraismax
		print,'total des frais <= ',fraismax
		ecrit_equation,abs(vi(0:ni-1))*pinitial(8,0:ni-1),fraismax,-1
	endif
print,'verification'
;pordres doit etre un vecteur de taille ni
bp=dblarr(neqtotal)
for i=0,neqtotal-1 do begin
	bp(i)=total(a(i,0:ni-1)*pordres(0:ni-1))
	case 1 of 
		eq_type(i) eq 1:begin
					contype=' >='
					if eq_sgn(i)*bp(i) ge eq_sgn(i)*b(i) then  resutest='ok' else $
	if eq_sgn(i)*bp(i) ge eq_sgn(i)*b(i)*(1. - eq_sgn(i)*epsilon) then resutest ='ok de justesse' else resutest=' contrainte non verifiee !'
				end
		eq_type(i) eq -1 :begin
					contype=' <='
					if eq_sgn(i)*bp(i) le eq_sgn(i)*b(i) then resutest ='ok' else $
	if eq_sgn(i)*bp(i) le eq_sgn(i)*b(i)*(1. + eq_sgn(i)*epsilon) then resutest ='ok de justesse' else resutest=' contrainte non verifiee !'
				end
	endcase
	print,'equation  ',i,'  :  a*x=',eq_sgn(i)*bp(i),contype,' b=',eq_sgn(i)*b(i),' -->',resutest
endfor	
end
 
pro opticac
supportbenef=[1800.,1900.]
volbenef=[-0.075,0.075]
supportrisque=[1700,2000]
volrisque=[-0.15,0.15]
s0=1800
nbmax=500.
sumrisk=100000.
pcac=[	[-100	 ,3 ,0    ,jj(300393) ,0.  ,1,200,1000,0.0002],$
	[0	 ,1 ,1850 ,jj(300393) ,0.44,1,200,100,0.01],$
	;[0	 ,1 ,1875 ,jj(300393) ,0.42,1,200,100,0.01],$
	;[0	 ,1 ,1900 ,jj(300393) ,0.42,1,200,100,0.01],$
	;[0	 ,1 ,1950 ,jj(300393) ,0.41,1,200,100,0.01],$
	;[0	 ,1 ,1975 ,jj(300393) ,0.38,1,200,100,0.01],$
	;[0	 ,1 ,2000 ,jj(300393) ,0.3628,1,200,200,0.01],$
	;[0	 ,1 ,2025 ,jj(300393) ,0.36,1,200,400,0.01],$
	;[0	 ,2 ,1650 ,jj(300393) ,0.6081,1,200,300,0.01],$
	;[0	 ,2 ,1700 ,jj(300393) ,0.592,1,200,300,0.01],$
	[-40	 ,2 ,1800 ,jj(300393) ,0.5485,1,200,300,0.01],$
	;[0	 ,2 ,1850 ,jj(300393) ,0.55,1,200,300,0.01],$
	;[0	 ,2 ,1750 ,jj(300493) ,0.369,1,200,300,0.01],$
	;[0	 ,2 ,1800 ,jj(300493) ,0.3635,1,200,300,0.01],$
	[0	 ,1 ,1850 ,jj(300593) ,0.44,1,200,100,0.01],$
	;[0	 ,1 ,1875 ,jj(300593) ,0.42,1,200,100,0.01],$
	;[0	 ,1 ,1900 ,jj(300593) ,0.42,1,200,100,0.01],$
	;[0	 ,1 ,1950 ,jj(300593) ,0.41,1,200,100,0.01],$
	;[0	 ,1 ,1975 ,jj(300593) ,0.38,1,200,100,0.01],$
	;[0	 ,1 ,2000 ,jj(300593) ,0.3628,1,200,200,0.01],$
	;[0	 ,1 ,2025 ,jj(300593) ,0.36,1,200,400,0.01],$
	;[0	 ,2 ,1650 ,jj(300593) ,0.6081,1,200,300,0.01],$
	;[0	 ,2 ,1700 ,jj(300593) ,0.592,1,200,300,0.01],$
	;[0	 ,2 ,1800 ,jj(300593) ,0.5485,1,200,300,0.01],$
	[0	 ,2 ,1850 ,jj(300593) ,0.55,1,200,300,0.01]]

volminb=volbenef(0) & volmaxb=volbenef(1) & sminb=supportbenef(0) & smaxb=supportbenef(1)
volminr=volrisque(0) & volmaxr=volrisque(1) & sminr=supportrisque(0) & smaxr=supportrisque(1)
db = [[sminb,smaxb],[volminb,volmaxb]]
dr = [[sminr,smaxr],[volminr,volmaxr]]
deltav0=0.
r=0.097
nvolr=4 & nsr=4
nvolb=2 & nsb=2
deltajr=3
deltajb=5
;calcul de la date d'aujourd'hui
dateevaluation=jj(010170)+int(systime(1)/86400.)
pcac1=pexpand(pcac)
res=poptimisation(pcac1,s0,deltav0,r,db,dr,nvolr,nsr,nvolb,nsb,deltajb,deltajr,dateevaluation,sumrisk,0,ginf=-7.)
tr=size(res)
if tr(0) ne 0 then begin
	res1=pcompress(res)
	;p_print2f,pcac,res1
endif
pordres1=[103.5305336146,0.,100.,0.,0.,140.5263657599,0.,100.,141.3313604216,0.]
;testequations,pordres1,pcac1,s0,deltav0,r,db,dr,nvolr,nsr,nvolb,nsb,deltajb,deltajr,dateevaluation,sumrisk,0,ginf=-7.
end


pro optinikkei
supportbenef=[1.6000,1.9000]
volbenef=[-0.05,0.05]
supportrisque=[1.4000,2.2000]
volrisque=[-0.15,0.15]
s0=1.7240
nbmax=500.
sumrisk=50000.
pcac=[	[0	 ,3 ,0    ,0          ,0.,1,500,1000],$
	[0	 ,1 ,1.6000 ,jj(120493) ,0.3576,2,500,500],$
	[0	 ,1 ,1.6500 ,jj(120493) ,0.3480,2,500,500],$
	[0	 ,1 ,1.7000 ,jj(120493) ,0.3393,2,500,1000],$
	[0	 ,1 ,1.7500 ,jj(120493) ,0.3287,2,500,1000],$
	[0	 ,1 ,1.8000 ,jj(120493) ,0.3147,2,500,1000],$
	[0	 ,1 ,1.8500 ,jj(120493) ,0.3087,2,500,500],$
	[0	 ,1 ,1.9000 ,jj(120493) ,0.2986,2,500,500],$
	[0	 ,1 ,1.9500 ,jj(120493) ,0.3015,2,500,500],$
	[0	 ,1 ,1.6000 ,jj(100593) ,0.3482,2,500,500],$
	[0	 ,1 ,1.6500 ,jj(100593) ,0.3382,2,500,500],$
	[0	 ,1 ,1.7000 ,jj(100593) ,0.3257,2,500,1000],$
	[0	 ,1 ,1.7500 ,jj(100593) ,0.3166,2,500,1000],$
	[0	 ,1 ,1.8000 ,jj(100593) ,0.3027,2,500,1000],$
	[0	 ,1 ,1.8500 ,jj(100593) ,0.2951,2,500,500],$
	[0	 ,1 ,1.9000 ,jj(100593) ,0.2873,2,500,500],$
	[0	 ,1 ,1.9500 ,jj(100593) ,0.2934,2,500,500],$
	[0	 ,1 ,1.6000 ,jj(140193) ,0.3483,2,500,500],$
	[0	 ,1 ,1.6500 ,jj(140193) ,0.3503,2,500,500],$
	[0	 ,1 ,1.7000 ,jj(140193) ,0.3420,2,500,1000],$
	[0	 ,1 ,1.7500 ,jj(140193) ,0.3242,2,500,1000],$
	[0	 ,1 ,1.8000 ,jj(140193) ,0.3261,2,500,1000],$
	[0	 ,1 ,1.8500 ,jj(140193) ,0.3180,2,500,500],$
	[0	 ,1 ,1.9000 ,jj(140193) ,0.2990,2,500,500],$
	[0	 ,1 ,1.9500 ,jj(140193) ,0.3015,2,500,500],$
	[0	 ,2 ,1.6000 ,jj(120493) ,0.3547,2,500,500],$
	[0	 ,2 ,1.6500 ,jj(120493) ,0.3466,2,500,500],$
	[0	 ,2 ,1.7000 ,jj(120493) ,0.3389,2,500,1000],$
	[0	 ,2 ,1.7500 ,jj(120493) ,0.3291,2,500,1000],$
	[0	 ,2 ,1.8000 ,jj(120493) ,0.3162,2,500,1000],$
	[0	 ,2 ,1.8500 ,jj(120493) ,0.3116,2,500,500],$
	[0	 ,2 ,1.9000 ,jj(120493) ,0.3039,2,500,500],$
	[0	 ,2 ,1.9500 ,jj(120493) ,0.3106,2,500,500],$
	[0	 ,2 ,1.6000 ,jj(100593) ,0.3446,2,500,500],$
	[0	 ,2 ,1.6500 ,jj(100593) ,0.3363,2,500,500],$
	[0	 ,2 ,1.7000 ,jj(100593) ,0.3251,2,500,1000],$
	[0	 ,2 ,1.7500 ,jj(100593) ,0.3172,2,500,1000],$
	[0	 ,2 ,1.8000 ,jj(100593) ,0.3047,2,500,1000],$
	[0	 ,2 ,1.8500 ,jj(100593) ,0.2982,2,500,500],$
	[0	 ,2 ,1.9000 ,jj(100593) ,0.2882,2,500,500],$
	[0	 ,2 ,1.9500 ,jj(100593) ,0.2964,2,500,500],$
	[0	 ,2 ,1.6000 ,jj(140193) ,0.3312,2,500,500],$
	[0	 ,2 ,1.6500 ,jj(140193) ,0.3350,2,500,500],$
	[0	 ,2 ,1.7000 ,jj(140193) ,0.3279,2,500,1000],$
	[0	 ,2 ,1.7500 ,jj(140193) ,0.3104,2,500,1000],$
	[0	 ,2 ,1.8000 ,jj(140193) ,0.3122,2,500,1000],$
	[0	 ,2 ,1.8500 ,jj(140193) ,0.3036,2,500,500],$
	[0	 ,2 ,1.9000 ,jj(140193) ,0.3171,2,500,500],$
	[0	 ,2 ,1.9500 ,jj(140193) ,0.2855,2,500,500]]


	
volminb=volbenef(0) & volmaxb=volbenef(1) & sminb=supportbenef(0) & smaxb=supportbenef(1)
volminr=volrisque(0) & volmaxr=volrisque(1) & sminr=supportrisque(0) & smaxr=supportrisque(1)
db = [[sminb,smaxb],[volminb,volmaxb]]
dr = [[sminr,smaxr],[volminr,volmaxr]]
deltav0=0.
r=0.04
nvolr=8 & nsr=17
nvolb=5 & nsb=12
deltajr=5
deltajb=5
;calcul de la date d'aujourd'hui
dateevaluation=jj(010170)+int(systime(1)/86400.)
pcac1=pexpand(pcac)
res=poptimisation(pcac1,s0,deltav0,r,db,dr,nvolr,nsr,nvolb,nsb,deltajb,deltajr,$
			dateevaluation,sumrisk,1,nbmax=5000,pinf=-0.1,tinf=-0.1)
res1=pcompress(res)
print,'------------portefeuille optimal--------------'
p_print,res1
print,'----------------------------------------------'
end



