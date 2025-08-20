; resolution du simplexe

	; verification			;initialisations
;	x_verif=dblarr(m+1)
;	xx_verif=dblarr(m+1,n+1)
;	zcw_verif=dblarr(n+1)
;	zcn_verif=dblarr(n+1)
;	abar=dblarr(m+1,m+n+1)
;	abar(*,*)=0
;	for i=1,m do abar(i,i)=1
;	abar(1:m,m+1:m+n)=a
;print,b1(b(findgen(n)+1))            ; b1 o b et b o b1 doivent etre egal a findgen(m+n)+1
;print,b(b1(findgen(n)+1))
;for p=1,m do begin                   ; x_verif doit etre egal a 0
;	som=0
;	for i=1,m do som=som+xp(i)*abar(p,bp(i))
;	x_verif(p)=som-b0(p-1)
;endfor
;print,'x_verif=',x_verif
;for i=1,n do begin                   ; xx_verif doit etre egal a 0
;	for t=1,m do begin
;		som=0
;		for p=1,m do som=som+xxp(p,i)*abar(t,bp(p))
;		xx_verif(t,i)=som-abar(t,bp(i+m))
;	endfor
;endfor
;print,'xx_verif='
 
;print,transpose(xx_verif)
;for p=1,n do begin
;	som=0
;	for t=1,m do if bp(t) le m then som=som+xxp(t,p)
;	if bp(p+m) le m then som=som-1.
;	zcw_verif(p)=som-zcwp(p)
;endfor
;for p=1,n do begin
;	som=0
;	for t=1,m do if bp(t) gt m then som=som+xxp(t,p)*c0(bp(t)-(m+1))
;	if bp(p+m) gt m then som=som-c0(bp(p+m)-(m+1))
;	zcn_verif(p)=som-zcnp(p)
;endfor
;print,'zcw_verif=',zcw_verif
;print,'zcn_verif=',zcn_verif
;som1=0				;calcul du z0 initial
;som2=0
;for i=1,m do if b(i) gt m then som1=som1+x(i)*c0(b(i)-(m+1)) else som2=som2+x(i)
;z0n=som1
;z0w=som2
;print,'z0  recalcule initial :',z0w,' W + ',z0n
;
;print,'incremental prevu: ',-mintheta*zcw(k),' W + ',-mintheta*zcn(k)
;
;print,'------------------------------------------------------------------------------------------------------'
;som1=0
;som2=0
;for i=1,m do if b(i) gt m then som1=som1+x(i)*c0(b(i)-(m+1)) else som2=som2+x(i)
;z0n=som1
;z0w=som2
;print,'z0 initial recalcule :',z0w,' W + ',z0n
;som1=0				; calcul du z0 final
;som2=0
;for i=1,m do if bp(i) gt m then som1=som1+xp(i)*c0(bp(i)-(m+1)) else som2=som2+xp(i)
;z0np=som1
;z0wp=som2
;z0w=z0w-mintheta*zcw(k)
;z0n=z0n-mintheta*zcn(k)
;print,'z0  final:',z0w,' W + ',z0n




function print_erreur,str
print,str
return,1
end

function bisimplexe,a,b0,b02,c0
; a(m,n) represente la matrice des coefficients des egalites deffinissant les contraintes
; b0(m) represente le vecteurs second membre des des egalites ci-dessus definies
; b02(m) represente le vecteurs second membre des des egalites ci-dessus pour un autre cas test 
; c0(n) est la forme lineaire a minimser
n=n_elements(c0)
m=n_elements(b0)
size_a=size(a)
epsilon=1.e-7
flagerreur=0
if size_a(0) ne 2 or size_a(1) ne m or size_a(2) ne n then flagerreur=print_erreur('erreur dans les dimensions')
for i=0,m-1 do if b0(i) lt 0 then flagerreur=print_erreur('erreur b('+string(i)+') n est pas positif')
if flagerreur eq 0 then begin
	npivot=0
	x=dblarr(m+1)
	x2=dblarr(m+1)
	xp=x
	xp2=x2
	x(1:m)=b0
	x2(1:m)=b02
	xx=dblarr(m+1,n+1)
	xx2=dblarr(m+1,n+1)
	xxp=xx
	xxp2=xx2
	xx(1:m,1:n)=a
	xx2(1:m,1:n)=a
	zcn=dblarr(n+1)
	zcn2=dblarr(n+1)
	zcnp=zcn
	zcnp2=zcn2
	zcw=dblarr(n+1)
	zcw2=dblarr(n+1)
	zcwp=zcw
	zcwp2=zcw2
	zcn(1:n)=-c0
	zcn2(1:n)=-c0
	for i=1,n do begin
		som=0
		for t=1,m do som=som+xx(t,i)
		zcw(i)=som
	endfor
	for i=1,n do begin
		som=0
		for t=1,m do som=som+xx2(t,i)
		zcw2(i)=som
	endfor
	b=findgen(m+n+1)
	b2=findgen(m+n+1)
	b1=b
	b12=b2
	bp=b
	bp2=b2
	b1p=b
	n_verif=0
	n_verif_max=100
	boucle:
	zct=zcw*(max(abs(zcn))>1.)+zcn
	zct2=zcw2*(max(abs(zcn2))>1.)+zcn2
;
;regle classique
;
 	kindex=where(zct(1:n) gt epsilon ,ctk)
 	kindex2=where(zct2(1:n) gt epsilon ,ctk2)
	if ctk eq 0 then goto,fin1 else begin		; fin1 = contraintes inpossibles ou fin d algorithme
	kindex=kindex+1
	kindex2=kindex2+1
print,'etape :',i,' norml=',kindex,' devie=',kindex2
		zcnm=max(zct(kindex),ordermax)
;print,where(zct(kindex) eq zcnm)
		k=kindex(ordermax)
	endelse
 	j=-1
	for i=1,m do if xx(i,k) gt 0  then if j eq -1 then begin
					j=i
					mintheta=x(i)/xx(i,k)
				endif else begin
					theta=x(i)/xx(i,k)
					if theta le mintheta then begin
						j=i
						mintheta=theta
					endif
				endelse
	if j eq -1 then goto,fin2			; fin2 = direction infinie

;print,'pivot',j,k
mintheta2=x2(j)/xx2(j,k)
	;upgrade
	bp=b		;
	bp2=b2		;
	bp(j)=b(k+m)	;
	bp2(j)=b2(k+m)	;
	bp(k+m)=b(j)	;
	bp2(k+m)=b2(j)	;
	b1p=b1
	b1p2=b12
	for p=1,m+n do begin
		if b1(p) eq k+m then b1p(p)=j
		if b1(p) eq j then b1p(p)=k+m
	endfor
	for p=1,m+n do begin
		if b12(p) eq k+m then b1p2(p)=j
		if b12(p) eq j then b1p2(p)=k+m
	endfor
	xjk=xx(j,k)
	xjk2=xx2(j,k)
	for i=1,m do xp(i)=x(i)-mintheta*xx(i,k)
	for i=1,m do xp2(i)=x2(i)-mintheta2*xx2(i,k)
	xp(j)=mintheta
	xp2(j)=mintheta2
	for i=1,m do for p=1,n do xxp(i,p)=xx(i,p)-xx(j,p)*xx(i,k)/xjk
	for i=1,m do for p=1,n do xxp2(i,p)=xx2(i,p)-xx2(j,p)*xx2(i,k)/xjk2
	for p=1,n do xxp(j,p)=xx(j,p)/xjk
	for p=1,n do xxp2(j,p)=xx2(j,p)/xjk2
	for i=1,m do xxp(i,k)=-xx(i,k)/xjk
	for i=1,m do xxp2(i,k)=-xx2(i,k)/xjk2
	xxp(j,k)=1./xjk
	xxp2(j,k)=1./xjk2
	for i=1,n do begin
		zcnp(i)=zcn(i)-xx(j,i)/xjk*zcn(k)
		zcwp(i)=zcw(i)-xx(j,i)/xjk*zcw(k)
	endfor	
	for i=1,n do begin
		zcnp2(i)=zcn2(i)-xx2(j,i)/xjk2*zcn2(k)
		zcwp2(i)=zcw2(i)-xx2(j,i)/xjk2*zcw2(k)
	endfor
	zcnp(k)=-zcn(k)/xjk
	zcnp2(k)=-zcn2(k)/xjk2
	zcwp(k)=-zcw(k)/xjk
	zcwp2(k)=-zcw2(k)/xjk2
	npivot=npivot+1
n_verif=n_verif+1
if n_verif mod 100 eq 0 then print,'pivot:',n_verif
	b=bp
	b2=bp2
	b1=b1p
	b12=b1p2
	x=xp
	x2=xp2
	xx=xxp
	xx2=xxp2
	zcn=zcnp
	zcn2=zcnp2
	zcw=zcwp
	zcw2=zcwp2
	goto,boucle
fin1:aindexes=where(b(1:m) le m,act)
restea=b(m+1:m+n)
resti=sort(restea)
;print,'reste=',restea(resti)
zea=zct(1:n)
;print,'zct=',zea(resti)
;print,'reste ex=',b(m+1:m+n)
;print,'zct=',zct
	resultat=0
	if act gt 0 then begin
		flagerreur=print_erreur(' impossible de realiser les contraintes')
		res=dblarr(5,n+1)
		lignes=b(aindexes+1)
		mm=n_elements(lignes)
		res(1,0)=mm
		res(1,1:mm)=lignes
		res(0,0)=-1
		print,'nombre de pivots:',n_verif
		return,res
	endif else begin
	resultat=dblarr(n)
	for i=1,m do resultat(b(i)-(m+1))=x(i)
	;print,'resultat=',resultat
	res=dblarr(5,n+1)
	res(4,1:n)=resultat
	res(0,0)=0
	print,'nombre de pivots:',n_verif
	return,res
	endelse
fin2:flagerreur=print_erreur('direction infinie dans la colonne '+string(b1(k)-m))
	res=dblarr(5,n+1)
	res(0,0)=1
	print,'nombre de pivots:',n_verif
	return,res
endif
end

function pre_simplexe,matrice,nsup0,ninf0
;prepare un probleme de simplexe en introduisant les "slack variables"
;nsup est le nombre d inegalite de type b+...>= 0  ou b>=0 presente en tete de la matrice
;ninf est le nombre d inegalite de type b+ ...<= 0 ou b>=0 presente a la suite
; on supose que les equations suivantes sont de la forme b + ... =0 ou b>=0 
s_mat=size(matrice)
m=s_mat(1)-1
n1=s_mat(2)-1
pe=0
if nsup0+ninf0 gt m+1 then pe=print_erreur('la somme des dimensions nsup+ninf est trop grande / matrice fournie')
if pe eq 0 then begin
	b=dblarr(m)
	for i=0,m-1 do b(i)=matrice(i+1,0)
	n2=n1+nsup0+ninf0
	a=dblarr(m,n2)
	a(*,*)=0
	a(0:m-1,0:n1-1)=-matrice(1:m,1:n1)
	if nsup0 gt 0 then for i=0,nsup0-1 do a(i,n1+i)=1.
	if ninf0 gt 0 then for i=0,ninf0-1 do a(nsup0+i,n1+nsup0+i)=-1.
	d=dblarr(n2)
	d(0:n1-1)=matrice(0,1:n1)
	resultat=simplexe(a,b,d)
	return,resultat
endif
end


;test du simplexe exemple p. 16 de Discrete optimisation in Pascal
; le resultat doit etre
;w=[0,1,2,3]
;x=[0.,2.5,2.5,2.5]
pro tests1
a=[	[0.,-1.,-2.,-3.,1.],$
	[15.,-1.,-2.,-3.,0.],$
	[20.,-2.,-1.,-5.,0.],$
	[10.,-1.,-2.,-1.,-1.]]
res=pre_simplexe(transpose(a),0,0)
print,transpose(res(4,*))
end

;test du simplexe exemple p. 313 de Numerical Recipe en Fortran
; le resultat doit etre 0 , 3.25 , 4.725 , 0.95
pro tests2
a=[[  0, -1., -1., -3.,0.5],$ ;  maximiser 	      1*x1 +1*x2 +3*x3 -.5*x4		 
   [740,-1, 0,-2,   0],$ ; 		740 - 1*x1 +0*x2 -2*x3 +0*x4  >= 0	1     I
   [  0, 0,-2, 0,  +7],$ ; 		0   + 0*x1 -2*x2 +0*x3 +7*x*4 >= 0	2     Im1=2
   [0.5, 0,-1,+1,  -2],$ ; 		0.5 - 0*x1 -1*x2 +1*x3 -2*x4  <= 0	3        Im2=1
   [  9,-1,-1,-1,  -1]]  ; 		9   - 1*x1 -1*x2 -1*x3 -1*x4   = 0	4           Im3=1
res=pre_simplexe(transpose(a),2,1)
print,transpose(res(4,*))
end
 
;test de chvatal cycling p31
;cyclage
pro tests3
a=[	[0.,-10.,57.,9.,24.],$
	[0, -0.5, 5.5, 2.5, -9.],$
	[0, -0.5, 1.5, 0.5, -1.],$
	[1., -1.,0.,0.,0.]]
 res=pre_simplexe(transpose(a),3,0)
 print,transpose(res(4,*))
end

;test de chvatal p 19
;la solution est 32/29=1.10345,8/29=0.275862,30/29=1.03448
pro tests4
a=[	[0,-5.,-5.,-3.],$
	[3.,-1.,-3.,-1.],$
	[2.,1.,0.,-3.],$
	[4.,-2.,1.,-2.],$
	[2.,-2.,-3.,1.]]
res=pre_simplexe(transpose(a),4,0)
 print,transpose(res(4,*))
end

; pb de chvatal p 69
;la solution est 0,0,5/2,7/2,0,1/2
 pro tests5
a=[	[0,4.,5.,1.,3.,-5.,8.],$
	[-1.,1.,0.,-4.,3.,1.,1.],$
	[-4.,5.,3.,1.,0.,-5.,3.],$
	[-4.,4.,5.,-3.,3.,-4.,1.],$
	[-5.,0.,-1.,0.,2.,1.,-5.],$
	[-7.,-2.,1.,1.,1.,2.,2.],$
	[-5.,2.,-3.,2.,-1.,4.,5.]]
res=pre_simplexe(transpose(-a),6,0)
 print,transpose(res(4,*))
end
