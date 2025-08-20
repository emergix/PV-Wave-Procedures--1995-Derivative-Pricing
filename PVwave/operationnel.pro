; fichier operationnel
; necessite option.pro

; description de la structure des marches

function marche_devise,marche
case 1 of
marche eq 'cac40' : devise='frf'
marche eq 'ftse': devise='gbp'
marche eq 'sp500' : devise='usd'
marche eq 'ibex' : devise='esp'
marche eq 'dax30' : devise='dem'
marche eq 'nikkei225' : devise='jpy'
else :print,'marche inconnu'
endcase
return,devise
end

function marche_spot,marche
case 1 of
marche eq 'cac40' : spot=2035.
marche eq 'ftse': spot=2922.
marche eq 'sp500' : spot=452.
marche eq 'ibex' : spot=2640.
marche eq 'dax30' : spot=1707.
marche eq 'nikkei225' : spot=18065.
else :print,'marche inconnu'
endcase
return,devise
end



; fournit la correlation annuelle entre les marche
; c est a dire que la covariance entre les logarithme des marche annuellement est donne
; par correlation * vol1 * vol2
function devise_correlation,marche1,devise2
case 1 of
marche1 eq 'cac40' and devise2 eq 'frf': correlation=0.2
marche1 eq 'ftse' and devise2 eq 'frf': correlation=0.1
marche1 eq 'sp500' and devise2 eq 'frf': correlation=0.1
marche1 eq 'ibex' and devise2 eq 'frf': correlation=0.1
marche1 eq 'dax30' and devise2 eq 'frf': correlation=0.1
marche1 eq 'nikkei225' and devise2 eq 'frf': correlation=0.1
marche1 eq 'cac40' and devise2 eq 'gbp': correlation=0.1
marche1 eq 'ftse' and devise2 eq 'gbp': correlation=0.2
marche1 eq 'sp500' and devise2 eq 'gbp': correlation=0.1
marche1 eq 'ibex' and devise2 eq 'gbp': correlation=0.1
marche1 eq 'dax30' and devise2 eq 'gbp': correlation=0.1
marche1 eq 'nikkei225' and devise2 eq 'gbp': correlation=0.1
marche1 eq 'cac40' and devise2 eq 'usd': correlation=0.1
marche1 eq 'ftse' and devise2 eq 'usd': correlation=0.1
marche1 eq 'sp500' and devise2 eq 'usd': correlation=0.2
marche1 eq 'ibex' and devise2 eq 'usd': correlation=0.1
marche1 eq 'dax30' and devise2 eq 'usd': correlation=0.1
marche1 eq 'nikkei225' and devise2 eq 'usd': correlation=0.1
marche1 eq 'cac40' and devise2 eq 'esp': correlation=0.1
marche1 eq 'ftse' and devise2 eq 'esp': correlation=0.1
marche1 eq 'sp500' and devise2 eq 'esp': correlation=0.2
marche1 eq 'ibex' and devise2 eq 'esp': correlation=0.1
marche1 eq 'dax30' and devise2 eq 'esp': correlation=0.1
marche1 eq 'nikkei225' and devise2 eq 'esp': correlation=0.1
marche1 eq 'cac40' and devise2 eq 'dem': correlation=0.1
marche1 eq 'ftse' and devise2 eq 'dem': correlation=0.1
marche1 eq 'sp500' and devise2 eq 'dem': correlation=0.1
marche1 eq 'ibex' and devise2 eq 'dem': correlation=0.1
marche1 eq 'dax30' and devise2 eq 'dem': correlation=0.2
marche1 eq 'nikkei225' and devise2 eq 'dem': correlation=0.1
marche1 eq 'cac40' and devise2 eq 'jpy': correlation=0.1
marche1 eq 'ftse' and devise2 eq 'jpy': correlation=0.1
marche1 eq 'sp500' and devise2 eq 'jpy': correlation=0.1
marche1 eq 'ibex' and devise2 eq 'jpy': correlation=0.1
marche1 eq 'dax30' and devise2 eq 'jpy': correlation=0.1
marche1 eq 'nikkei225' and devise2 eq 'jpy': correlation=0.2

else : begin
	print,'correlation entre devises non connue ,prise egale a 0.'
	correlation=0.
	end
endcase
return,correlation
end

;volatilite historique des devises
function devise_volatilite,devise
case 1 of
devise eq 'frf' : vol=10.
devise eq 'gbp' : vol=10.
devise eq 'usd' : vol=10.
devise eq 'esp' : vol=10.
devise eq 'jpy' : vol=10.
devise eq 'dem' : vol=10.
else : print,'devise inconnue pour la volatilite historique:',devise
endcase
return,vol
end

function courbetaux,devise
case 1 of
devise eq 'frf' : cctaux=transpose([$
				[1/365.,	8.375],$
				[1/12.,		8.0625],$
				[2/12.,		7.8125],$
				[3/12.,		7.75],$
				[6/12.,		7.375],$
				[9/12.,		7.125],$
				[1.,		7.125],$
				[2.,		7.21],$
				[3.,		7.23],$
				[4.,		7.30],$
				[5.,		7.37],$
				[7.,		7.60],$
				[8.,		7.68],$
				[9.,		7.73],$
				[10.,		7.80]])
devise eq 'gbp' : cctaux=transpose([$
				[1/365.,	5.98],$
				[1/12.,		5.93],$
				[2/12.,		5.93],$
				[3/12.,		5.93],$
				[6/12.,		5.87],$
				[9/12.,		5.86],$
				[1.,		5.86],$
				[2.,		5.8],$
				[3.,		5.8],$
				[4.,		5.8],$
				[5.,		5.8],$
				[7.,		5.8],$
				[8.,		5.8],$
				[9.,		5.8],$
				[10.,		5.8]])
devise eq 'usd' : cctaux=transpose([$
				[1/365.,	3.12],$
				[1/12.,		3.02],$
				[2/12.,		3.12],$
				[3/12.,		3.18],$
				[6/12.,		3.25],$
				[9/12.,		3.45],$
				[1.,		3.56],$
				[2.,		4.10],$
				[3.,		4.5],$
				[4.,		5.05],$
				[5.,		5.30],$
				[7.,		5.65],$
				[8.,		5.78],$
				[9.,		5.85],$
				[10.,		5.95]])
devise eq 'esp' : cctaux=transpose([$
				[1/365.,	10.3],$
				[1/12.,		10.3],$
				[2/12.,		10.3],$
				[3/12.,		10.3],$
				[6/12.,		10.3],$
				[9/12.,		10.3],$
				[1.,		10.3],$
				[2.,		10.3],$
				[3.,		10.3],$
				[4.,		11.3],$
				[5.,		11.5],$
				[7.,		11.5],$
				[8.,		11.5],$
				[9.,		11.5],$
				[10.,		11.5]])
devise eq 'dem' : cctaux=transpose([$
				[1/365.,	8.75],$
				[1/12.,		8.37],$
				[2/12.,		8.07],$
				[3/12.,		7.93],$
				[6/12.,		7.56],$
				[9/12.,		7.12],$
				[1.,		6.93],$
				[2.,		6.9],$
				[3.,		6.9],$
				[4.,		6.5],$
				[5.,		6.5],$
				[7.,		6.5],$
				[8.,		6.5],$
				[9.,		6.5],$
				[10.,		6.5]])
devise eq 'jpy' : cctaux=transpose([$
				[1/365.,	3.],$
				[1/12.,		3.43],$
				[2/12.,		3.43],$
				[3/12.,		3.43],$
				[6/12.,		3.43],$
				[1.,		3.48],$
				[2.,		3.81],$
				[3.,		3.23],$
				[4.,		4.53],$
				[5.,		4.73],$
				[7.,		4.93],$
				[10.,		5.03]])
else :print,'devise inconnue:',devise
endcase
return,cctaux
end

function courbevol,marche
case 1 of
marche eq 'cac40' : ccvol=transpose([$
				[1/12.,		19.],$
				[3/12.,		18.2],$
				[6/12.,		17.7],$
				[1.,		18.5],$
				[3.,		20.5],$
				[5.,		20.]])
marche eq 'ftse' : ccvol=transpose([$
				[1/12.,		15.75],$
				[3/12.,		15.75],$
				[6/12.,		15.75],$
				[1.,		16.25],$
				[3.,		16.25],$
				[5.,		16.25]])
marche eq 'sp500' : ccvol=transpose([$
				[1/12.,		12.75],$
				[3/12.,		12.75],$
				[6/12.,		12.75],$
				[1.,		13.],$
				[3.,		14.],$
				[5.,		14.5]])
marche eq 'ibex' : ccvol=transpose([$
				[1/12.,		19.7],$
				[3/12.,		19.7],$
				[6/12.,		18.7],$
				[1.,		18.7],$
				[3.,		18.7],$
				[5.,		19.7]])
marche eq 'nikkei225' : ccvol=transpose([$
				[1/12.,		21.],$
				[3/12.,		21.],$
				[6/12.,		20.],$
				[1.,		19.5],$
				[3.,		17.75],$
				[5.,		17.25]])
marche eq 'dax30' : ccvol=transpose([$
				[1/12.,		16.],$
				[3/12.,		16.],$
				[6/12.,		15.5],$
				[1.,		15.],$
				[3.,		15.25],$
				[5.,		15.25]])

marche eq 'smi' : ccvol=transpose([$
				[1/12.,		15.5],$
				[3/12.,		14.5],$
				[6/12.,		14.5],$
				[1.,		14.5],$
				[3.,		16.],$
				[5.,		16.]])
else: print,'marche inconnu:',marche
endcase
return,ccvol
end

; taux d emprunt des titres en % par an
function taux_emprunt_titre,marche
case 1 of
marche eq 'cac40' :taux=0.5
marche eq 'ftse': taux=0.5
marche eq 'sp500' : taux=0.5
marche eq 'ibex' : taux=0.
marche eq 'dax30' : taux=0.5
marche eq 'nikkei225' : taux=1.
else :print,'marche inconnu:',marche
endcase
return,taux
end


function prevision_yield,marche,s0
case 1 of
marche eq 'cac40' :begin
cdiv1=transpose([$
[jj(140193),	.32],$
[jj(290193),	.16],$
[jj(010293),	.11],$
[jj(310393),	.4],$
[jj(280493),	.76],$
[jj(190593),	.37],$
[jj(250593),	1.54],$
[jj(010693),	3.97],$
[jj(050693),	4.08],$
[jj(090693),	2.1],$
[jj(100693),	.5],$
[jj(110693),	1.19],$
[jj(120693),	.32],$
[jj(150693),	.11],$
[jj(170693),	2.24],$
[jj(190693),	2.61],$
[jj(260693),	.89],$
[jj(290693),	2.83],$
[jj(300693),	.51],$
[jj(010793),	4.2],$
[jj(020793),	5.99],$
[jj(070793),	.57],$
[jj(030793),	1.13],$
[jj(060793),	3.45],$
[jj(100793),	1.27],$
[jj(200793),	.65],$
[jj(310793),	.37],$
[jj(291193),	.98]])
cdiv2=div_time_shift(cdiv1,1.,0.08) ;       on suppose un progression des dividendes de 8% par an 
cdiv3=div_time_shift(cdiv2,1.,0.08) ;
cdiv4=div_time_shift(cdiv3,1.,0.08) ;
cdiv5= transpose([$
[jj(100797),	42.*(1.08)^4],$
[jj(100798),	42.*(1.08)^5],$
[jj(100799),	42.*(1.08)^6],$
[jj(100700),	42.*(1.08)^7],$
[jj(100701),	42.*(1.08)^8],$
[jj(100702),	42.*(1.08)^9],$
[jj(100703),	42.*(1.08)^10],$
[jj(100704),	42.*(1.08)^11]])
ns1=size(cdiv1)
n1=ns1(1)
ns2=size(cdiv2)
n2=ns2(1)
ns3=size(cdiv3)
n3=ns3(1)
ns4=size(cdiv4)
n4=ns4(1)
ns5=size(cdiv5)
n5=ns5(1)
ccdiv=fltarr(n1+n2+n3+n4+n5,2)
ccdiv(0:n1-1,*)=cdiv1
ccdiv(n1:n1+n2-1,*)=cdiv2
ccdiv(n1+n2:n1+n2+n3-1,*)=cdiv3
ccdiv(n1+n2+n3:n1+n2+n3+n4-1,*)=cdiv4
ccdiv(n1+n2+n3+n4:n1+n2+n3+n4+n5-1,*)=cdiv5
end
marche eq 'ftse' :ccdiv=mult_div(s0,transpose([$
				[jj(060693)	,0.0],$
				[jj(060694)	,0.023],$
				[jj(060695)	,0.023],$
				[jj(060696)	,0.023],$
				[jj(060697)	,0.023],$
				[jj(060698)	,0.023],$
				[jj(060699)	,0.023],$
				[jj(060600)	,0.023],$
				[jj(060601)	,0.023],$
				[jj(060602)	,0.023]]))
marche eq 'dax30' :ccdiv=mult_div(s0,transpose([$
				[jj(060693)	,0.0],$
				[jj(060694)	,0.023],$
				[jj(060695)	,0.023],$
				[jj(060696)	,0.023],$
				[jj(060697)	,0.023],$
				[jj(060698)	,0.023],$
				[jj(060699)	,0.023],$
				[jj(060600)	,0.023],$
				[jj(060601)	,0.023],$
				[jj(060602)	,0.023]]))

marche eq 'sp500' :ccdiv=mult_div(s0,transpose([$
				[jj(060693)	,0.0],$
				[jj(060694)	,0.023],$
				[jj(060695)	,0.023],$
				[jj(060696)	,0.023],$
				[jj(060697)	,0.023],$
				[jj(060698)	,0.023],$
				[jj(060699)	,0.023],$
				[jj(060600)	,0.023],$
				[jj(060601)	,0.023],$
				[jj(060602)	,0.023]]))
marche eq 'ibex' :ccdiv=mult_div(s0,transpose([$
				[jj(060693)	,0.0],$
				[jj(060694)	,0.04],$
				[jj(060695)	,0.04],$
				[jj(060696)	,0.04],$
				[jj(060697)	,0.04],$
				[jj(060698)	,0.04],$
				[jj(060699)	,0.04],$
				[jj(060600)	,0.04],$
				[jj(060601)	,0.04],$
				[jj(060602)	,0.04]]))
marche eq 'nikkei225' :ccdiv=mult_div(s0,transpose([$
				[jj(060693)	,0.0],$
				[jj(060694)	,0.015],$
				[jj(060695)	,0.015],$
				[jj(060696)	,0.015],$
				[jj(060697)	,0.015],$
				[jj(060698)	,0.015],$
				[jj(060699)	,0.015],$
				[jj(060600)	,0.015],$
				[jj(060601)	,0.015],$
				[jj(060602)	,0.015]]))


endcase
sd=size(ccdiv)
n=sd(1)
for i=0,n-1 do ccdiv(i,1)=ccdiv(i,1)/s0
return,ccdiv
end
;


;
;			options simples
;
; les dates s ecrivent soit jj(230393) pour le 23 janvier 93 soit aa(0.25) pour dans 3 mois 


 ;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui()
 ;  calcul_put,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui()


;            americaine

 ;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),americain='oui'

 ;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),debut_americain=jj(310394)


;         asiatique

 ;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),asiatique=20./365.25   (moyenne sur 20 jours)
; ou asiatique=aa(2) pour deux ans 
; on peut preciser la valeur de lamoyenne au demarrage de l option sinon par defaut elle prise egale a s0
; cette moyenne si asiatique > top-t0 represente ne fait la moyenne des cours deja connu qui entreront dans la composition
; de la moyenne future  
;       ex:                moyenne=1800.

;               garantie de change

; on peut exiger que la valeur d une option soient garantie dans un autre monnaie
; que celle dans laquelle se gere la position du support
; on emploie alors le mot cle    garantie_change='dem' pour une option a support en frf
;              ou garantie_change='frf' pour une option a support etranger.

;;;;;;;; les option suivante sont exclusives l une de l autre

; 		avec resets

  ;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310195),reset=[jj(010294),jj(010594),jj(010794)]


;  		explosive

;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),plafond_explosif=1800
;  calcul_put,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),plancher_explosif=1600
;
; avec evidemment plancher explosif < k < plafond explosif

; 		 avec des paliers


;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310195),paliers=[1800.,1900.],debut_paliers=jj(310694)

; 		 spectrale (on fournit des paliers et des coefficients relatif aux differentes tranches)


;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310195),t0=aujoudhui(),paliers=[1800.,1900.],debut_paliers=jj(310694),spectre=[0.5,2.,1.]

;		rappelable par l emetteur (formule encore discutable :employer plutot un spread d option)

;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),niveau_rappel=1800



;                lookback

;  calcul_call,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),lookback='min'    droit d acheter au plus bas
;  calcul_put,'cac40',s0=1789,k=1700,top=jj(310194),t0=aujoudhui(),lookback='max'     droit de vendre au plus haut


;;;;;;;;;;; option sur deux supports

;  normale ou americaine  ces options se decline en quatre type :

; calcul_bioption,'call','min','cac40','SP500,s01=1900.,s02=3500.,k=1.1,top=aa(3.)
; calcul_bioption,'call','max','cac40','SP500,s01=1900.,s02=3500.,k=1.1,top=aa(3.) ,debut_americain=aa(1.5)
; calcul_bioption,'call','moyenne','cac40','SP500,s01=1900.,s02=3500.,k=1.1,top=aa(3.) ,debut_americain=aa(1.5),correlation=0.55
; calcul_bioption,'call','spread','cac40','SP500,s01=1900.,s02=3500.,k=1.1,top=aa(3.) ,debut_americain=aa(1.5),correlation=0.55

 
