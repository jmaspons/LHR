set.seed(1)
tf <- 20 # Final time
x0 <- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
x0L<- lapply(1:6, function(x) x0 * x)
params<- getParams.LH_Beh()
x0A<- c(N1s=1, N1b=1, N1bF=1, N1j=1, N2s=1, N2b=1, N2bF=1, N2j=1)
x0AL<- lapply(50, function(x) x0 * x)
i<- 1
trans<- transitionMat.LH_Beh(params[i,])
rateFunc.LH_Beh(x0L[[5]],params = params[i,])
for (i in 1:nrow(params)){
  tmp<- rateFunc.LH_Beh(x0AL[[1]],params = params[i,])
  tmp<- tmp[tmp != 0]
  repro<- sum(tmp[grep("repro[0-9]{1}b\\.", names(tmp))])
  dead<- sum(tmp[grep("dead", names(tmp))])
  tmp<- c(repro=repro, dead=dead)
  cat("\t", rownames(params)[i], "\n")
  print(rbind(tmp, tmp / sum(tmp)))
}

lapply()

sim0<- ssa.adaptivetau(x0L[[5]], trans, rateFunc.LH_Beh, params[i,], tf)
plotSim(sim0)
plotSim(sim0, groups = "habitat*age")

## exploreSSA
transitionMat=transitionMat.LH_Beh; rateFunc=rateFunc.LH_Beh
dtDiscretize<- 1
replicates<- 1000
discretePop=FALSE; finalPop=TRUE
burnin<- -1
cores<- 4
mc.preschedule<- TRUE

# For x0L<- lapply(1:6, function(x) x0 * x) and dtDiscretize<- 1, replicates<- 10000

params<- getParams.LH_Beh("slow", behavior="preferHab2", scenario="mortalHab2") # check extinctions
# user     system  elapsed 
# 2586.897   16.457  204.755  :  mc.preschedule<- TRUE
#1 731.257 3779.975 1436.439  :  mc.preschedule<- FALSE
params<- getParams.LH_Beh("fast") # check long simulation (many events)
# user         system   elapsed 
# 60349.784    64.692  6524.418  :  mc.preschedule<- TRUE
# 71526.617  6402.786  7880.173  :  mc.preschedule<- FALSE
params<- getParamsCombination()

system.time(res<- exploreSSA(x0L=x0L, params=params, transitionMat=transitionMat, rateFunc=rateFunc,
                             tf=tf, replicates=replicates,
                             discretePop=discretePop, finalPop=finalPop, burnin=burnin, dtDiscretize=dtDiscretize, cores=cores, mc.preschedule=mc.preschedule))


# Try to find deterministic character of the strategies in a simple demographic model.
params<- getParams()
demo1<- params[,c(1,3,4,6:8)]
demoFreq1<- demo1[,c(2,3:5)]
demoFreq1 / rowSums(demoFreq1)

# Prob events = b + da + dj + g = 1
b<- da<- dj<- g<- 1:10
params<- expand.grid(list(b=b, da=da, dj=dj, g=g))
parS<- params/rowSums(params)
parSR<- unique(round(parS, digits=1))
K<- 1000

params<- data.frame(parSR, K)


