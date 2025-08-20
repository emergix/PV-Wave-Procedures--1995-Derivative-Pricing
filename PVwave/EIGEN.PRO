pro testcov
openr,unit,'~olivier/ING/Portfolios/Covariance_Matrix/m1',/get_lun
rmf,unit,covariances,15,15
close,unit

values=PRINC_COMP(covariances)
print,'principal components:',values
eigenval=eig(covariances,/symmetric)
print,'eigen values:',eigenval
normp=max(covariances)
print,'estimated norm:',15*normp
end
