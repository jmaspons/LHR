pars<- sampleLH() 
obj<- LH(pars)
obj<- LH()
obj<- LH(obj[1:20,])
head(obj)
obj[1:10,]
obj$a[2]
obj[[2]]
obj[1]
