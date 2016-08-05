# Prob events = b + da + dj + g = 1
b<- da<- dj<- g<- 1:10
params<- expand.grid(list(b=b, da=da, dj=dj, g=g))
parS<- params/rowSums(params)
parSR<- unique(round(parS, digits=1))
params<- merge(parSR, data.frame(clutch=1:10))
K<- 1000
params<- data.frame(params, K)
params<- params[seq(1,4500, by=100),] #Sample for tests

x0 <- c(Na=1, Nj=1)
x0L<- lapply(round(K/4), function(x) x0 * x)
tf<- 5 # Final time

i<- 40
trans<- transitionMat.LH(params[i,])
rateFunc.LH(x0L[[1]], params[i,])

sim0<- ssa.adaptivetau(x0L[[1]], trans, rateFunc.LH, params[i,], tf)
sim0<- ssa.exact(x0L[[1]], trans, rateFunc.LH, params[i,], tf)
matplot(x = sim0[,1], y=sim0[,2:3])



replicates<- 1000
burnin=-1; dtDiscretize=NULL
cores=4; discretePop=FALSE; finalPop=TRUE; mc.preschedule=TRUE

res<- exploreSSA(x0L=x0L, params=params, transitionMat=transitionMat.LH, rateFunc=rateFunc.LH, tf=tf, replicates=replicates,
                 discretePop=discretePop, finalPop=finalPop, burnin=burnin, cores=cores, mc.preschedule=mc.preschedule)
x<- merge(res$params, res$stats, by.x="row.names", by.y="params")

plot(x[,c(2:6,9,12,16:17)])
