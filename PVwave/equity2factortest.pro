
pro testpl
rz=transpose([[1.,0.1],[2.,0.105],[3.,0.11],[4.,0.1125],[5.,0.115],[6.,0.1175],[7.,0.118]])
a=0.1
sigma=0.014
nb=8
vol=0.2
corr=+0.
s=100.
nbinstrument=5
portef=fltarr(80,nbinstrument)
;portef(1:17,0)=[201,0,2.5,0,0,0,0,0,0,0,0,0,0,0,0,3,1]
;portef(20:22,0)=[1,1.5,portef(1:17,0)=[600,0,5.,0,0,0,0,200,3,0,100.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,1)=[600,0,5.,0,0,0,0,200,-1,0,110.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,2)=[600,0,5.,0,0,0,0,200,-1,0,120.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,3)=[650,0,5.,0,0,0,0,200,5,0,80.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,4)=[650,0,5.,0,0,0,0,200,5,0,90.,0,0,0,0,0,0]  ;un call 5 an
portef(1:17,4)=[650,0,5.,0,0,0,0,200,0,0,110.,0,0,0,0,0,0]  ;un call 5 an
yieldlist=[0.00,0.00]
t_divlist=[1.5,2.5]
t=4.

;#######construction des arbres
crbtree=HWbuild_crbtree(rz,a,sigma,nb)			; arbre de taux
hwlnprob=HWLNtree(crbtree,vol,sigma,corr,s)			; arbre unifie index-taux
valorisation=HWLNOption(crbtree,hwlnprob,portef,nbinstrument,yieldlist,t_divlist)  ;valorisation d un portefeuille
print,'arbres construit'

;#######construction du risque probabiliste

;  @@@future dependance du portefeuill aux facteurs
portef_value=HWLN_portefvalue(valorisation,crbtree,t)
;  @@@implicit forward probability
probability=HWLN_probability(crbtree,hwlnprob,valorisation,rz,t)  
;  @@@distribution du P & L   
pnl=profit_loss_probability(portef_value,probability)	
;  @@@density of the P & L		
density=distribution_to_density(pnl)	
;  @@@density of the P & L  other definition		
density1=distribution_to_density1(pnl)	
;  @@@risk of the P & L	
risk=distribution_to_risk(pnl)	
		
; #### display
window,0,title='P&L on date = '+string(t)
surface,portef_value,ytitle='index',xtitle='short rate'
window,1,title='probability'
surface,probability,ytitle='index',xtitle='short rate'
window,2,title='P&L probability distribution'
plot,pnl.profit,pnl.probability
window,3,title='P&L probability density'
plot,density.profit,density.density
window,4,title='P&L probability density 1'
plot,density1.profit,density1.density
window,5,title='P&L probability density 1'
plot,risk.risk,risk.probability
end
