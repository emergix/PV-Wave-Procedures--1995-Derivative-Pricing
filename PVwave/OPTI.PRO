 ;**********************************************************************************

function funk,a
;la solution est probablement 0.5,1.5,2.5,3.5,....
n=n_elements(a)
s=0.
for i=1,n-1 do s=s+((a(i)-1)^2 + (a(i)-(i+1))^2)
;print,'arg=',a,'  return=',s,n
return,s
end

; implementation de amoeba p 411 de numerical recipes
; funk  est la function a optimiser qui prend en argument un vecteur de dimensionndim
function amoeba,p,y,ndim,ftol
psum=dblarr(ndim+1)
ptry=dblarr(ndim+1)
nfunk=0
mpts=ndim+1
nmax=5000
 ;get_psum
for j=1,ndim do begin
	sum=0.
	for i=1,mpts do begin
		sum=sum+p(i,j)
		psum(j)=sum
	endfor
endfor
;

while (1 eq 1) do begin
	ilo=1
	if y(1) gt y(2) then begin
		inhi=2
		ihi=1
	endif else begin
		inhi=1
		ihi=2
	endelse
	for i=1,mpts do begin
		if y(i) le y(ilo) then ilo=i
	 	if y(i) gt y(ihi) then begin
			inhi=ihi
			ihi=i
		endif else if ((y(i) gt y(inhi)) and (i ne ihi)) then inhi=i
	endfor
	rtol= 2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
print,'but=',ftol,'actual=',rtol,'sol=',y(ihi),y(ilo)
	if rtol lt ftol then  begin
		result=fltarr(ndim)
		result(0:ndim-1)=p(ilo,1:ndim)
		return,result
	endif
	if nfunk gt nmax then begin
		print,'nmax exceeded'
		return,y
	endif
	nfunk=nfunk+2
	 
	; amotry
	fac=-1.
	fac1=(1.-fac)/float(ndim)
	fac2=fac1-fac
	for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
	ytry=funk(ptry)								;appel a funk
	if  ytry lt y(ihi) then begin
		y(ihi)=ytry
		for j=1,ndim do begin
			psum(j)=psum(j)+ptry(j)-p(ihi,j)
			p(ihi,j)=ptry(j)
		endfor
	endif
	;
	if ytry lt y(ilo) then begin

		; amotry
		fac=2.
		ptry=fltarr(ndim+1)
		fac1=(1.-fac)/float(ndim)
		fac2=fac1-fac
		for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
		ytry=funk(ptry)							;appel a funk
		if  ytry lt y(ihi) then begin
			y(ihi)=ytry
			for j=1,ndim do begin
				psum(j)=psum(j)+ptry(j)-p(ihi,j)
				p(ihi,j)=ptry(j)
			endfor
		endif
		;
	endif else begin
		if ytry ge y(inhi) then begin
			ysave=y(ihi)
			; amotry
			fac=0.5
			fac1=(1.-fac)/float(ndim)
			fac2=fac1-fac
			for j=1,ndim do ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
			ytry=funk(ptry)						;appel a funk
			if  ytry lt y(ihi) then begin
				y(ihi)=ytry
				for j=1,ndim do begin
					psum(j)=psum(j)+ptry(j)-p(ihi,j)
					p(ihi,j)=ptry(j)
				endfor
			endif
			;
			if ytry ge ysave then begin
				for i=1,mpts do if i ne ilo then begin
							for j=1,ndim do begin
								psum(j)=0.5*(p(i,j)+p(ilo,j)) 
								p(i,j)=psum(j)
							endfor
							y(i)=funk(psum)		;appel a funk
						endif
				;
				nfunk=nfunk+ndim
				 ;get_psum
				for j=1,ndim do begin
					sum=0.
					for i=1,mpts do begin
						sum=sum+p(i,j)
						psum(j)=sum
					endfor
				endfor
				;
			endif
		endif else nfunk=nfunk-2
	endelse
endwhile
print,'nfunk=',nfunk
end
 	 
pro test_funk
ndim=3
y=dblarr(ndim+2)
p=dblarr(ndim+2,ndim+1)
for i=1,ndim do for j=1,ndim do  if i eq j then p(i+1,j)=p(0,j)+1. else p(i+1,j)=0.
ptry1=dblarr(ndim+1)
ptry1(1:ndim)=p(i,1:ndim)
for i=1,ndim+1 do begin
	ptry1(1:ndim)=p(i,1:ndim)
	y(i)=funk(ptry1) 
endfor
print,'y initial=',y
ftol=0.0000000001
print,amoeba(p,y,ndim,ftol)
end
