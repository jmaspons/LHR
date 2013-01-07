

library(LHR)
## Test cohortModel.c
######################
x<- cohortModel(5, .4, .1, .05, 3, 14, .3, .02, c(1,.8,.5))

## Small example
#################
{
n0<-30
maxPomited<- 0.01

survA<- 0.7	# supervivència dels adults
var.survA<- 0.05	# provabilitat de que un ou arribi a adult [Surv Egg->Adult]
broods<- 2		# nombre de postes per any
clutch<- 5		# mida de la posta
survJ<- 0.3
var.survJ<- .01
max.years<- NA
meanSeason<- .8
amplSeason<- .4
breedInterval<- 2

s<- survdist(n0=n0, survA, maxPomited=maxPomited)
clutch2<- rep(clutch %/% broods, broods)

if (clutch %% broods > 0) clutch2[1:(clutch %% broods)]<- clutch2[1:(clutch %% broods)] + 1

f<- fertdist(s, broods, clutch2, survJ, meanSeason=meanSeason, amplSeason=amplSeason, breedInterval=breedInterval)
R0<- fitnessdist(s, f, n0[1])
R02<- fitnessdistAIO(n0[1], survA, var.survA, broods=broods, clutch=clutch, survJ=survJ, var.survJ=var.survJ, maxPomited=maxPomited)

s<- survdist(n0=n0, survA, var.survA=var.survA, maxPomited=maxPomited)
clutch<- c(3,2)
f<- fertdist(s, broods, clutch, survJ, var.survJ=var.survJ, meanSeason=meanSeason, amplSeason=amplSeason, breedInterval=breedInterval)
R0<- fitnessdist(s, f, n0[1])
R02<- fitnessdistAIO(n0[1], survA, var.survA, broods=broods, clutch=clutch, survJ=survJ, var.survJ=var.survJ, maxPomited=maxPomited)

s<- survdist(n0=n0, survA, maxPomited=maxPomited)
f<- fertdist(s, broods, clutch, survJ)
R0<- fitnessdist(s, f, n0)

library(lattice)

setwd("../../../../doctorat/projectes/LHT/R")
png(filename="../../../congressos/neobiota2012-Pontevedra/imatges/modelDemo%d.png",width=900, height=900, pointsize=24)
xyplot(probS ~ years, data=s, ylab="Probability", type="s", main="Survival", col=1,
	par.settings=list(par.main.text=list(cex=4), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3)) , scales=list(cex=3), lwd=2)
xyplot(probF ~ fert, groups=years, data=f[f$years%%2 !=0,], ylab="Probability", xlab="offsprings", type="s", xlim=c(-1,50), main="Fertility",
	auto.key=list(space="right", title="Years\nlived"),
	par.settings=list(par.main.text=list(cex=4), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3)) , scales=list(cex=3), lwd=2)
xyplot(probR0 ~ R0, data=R0, type="s", ylab="Probability", xlab="offsprings", xlim=c(-1,50), main="Fitness",
	par.settings=list(par.main.text=list(cex=4), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3)) , scales=list(cex=3), lwd=3)
dev.off()
}

## Explore parameter space
###########################n
n0<- 3^(0:3) 	# mida de la població inicial
maxPomited<- 0.01
setwd("../../../../doctorat/projectes/LHT/R")
source("./libCreaEstrategiesLHT.r")

## Environmental stochasticity
envir<- data.frame(var=seq(0, 0.2, length=4))
envirJA<- expand.grid(varJ=envir$var, varA=envir$var)

## Seasonality scenarios
season<- expand.grid(amplSeason=seq(0, 1, length=3), breedInterval=1:3, broods=(1:6)[seq(1, 6, length=3)]) ## breedInterval=2^(0:2)
season<- season[season$broods * season$breedInterval < 12,] # remove combination which don't fit in one year
season<- season[-which(season$broods == 1 | season$amplSeason == 0)[-1],] # remove combinations whith 1 brood (allways at the optimum)
season$meanSeason<- 1 - season$amplSeason / 2

alignCriterion<- data.frame(alignCriterion=c("bestFirst", "bestMean"), stringsAsFactors=FALSE)

## regular sampling
strategy<- CreaEstrategiesLambdaLliure(dimLength=15, constrains=TRUE, separa_B_sJ=TRUE, rangLambda=rangeLambda, roundB=TRUE) ##PERFER: arreglar constrains a libCreaEstrategiesLHT.r
strategy<- strategy[strategy$lambda > min(ranges$lambda) & strategy$lambda < max(ranges$lambda),]

simulation<- merge(strategy, data.frame(n0=n0))
simulation<- merge(simulation, envirJA)
simulation<- merge(simulation, season)
simulation<- merge(simulation, alignCriterion)

simulation<- with(simulation, data.frame(n0, survA=sA, var.survA=varJ, broods=broods, B=B, survJ=sJ, var.survJ=varJ, meanSeason=meanSeason, amplSeason=amplSeason, breedInterval=breedInterval, alignCriterion=alignCriterion, stringsAsFactors=FALSE))

save.image("initialConditions.lambdaLliure.RData")

result<- runCohortModel(simulation, maxPomited=0.01)
x<- cbind(simulation, result)

## Latin hypercub sampling
nStrategies<- 5000


strategy<- CreaEstrategiesHiperCub(nStrategies, rangLambda=c(1,1.2), roundB=TRUE) # Creates a set of strategies with constraints
strategy$lifespan<- - 1 / log(strategy$sA)
strategy<- strategy[-which(strategy$lifespan < strategy$alpha),]
# str(strategy)
# 'data.frame':	x obs. of  6 variables:
#  $ lambda: num  0.9 0.9 0.901 0.902 0.905 ...		LAMBDA
#  $ sA    : num  0.794 0.726 0.625 0.785 0.435 ... 	ADULT SURVIVAL ** mandatory
#  $ F     : num  0.272 0.667 1.647 0.307 1.957 ... 	NET FERTILITY
#  $ sJ    : num  0.449 0.725 0.313 0.314 0.3 ... 	JUVENILE SURVIVAL **mandatory
#  $ b     : num  0.604 0.919 5.262 0.978 6.53 ... 	FERTILITY **mandatory
#  $ alpha : int  3 4 4 3 2 3 2 2 1 4 ... 		YEAR OF FIRST BREEDING
 


strategyName<- function(strategy, sep="\t"){
  name<- character()
  for (i in 1:ncol(strategy)){
    name<- paste(name, colnames(strategy)[i], "=", round(strategy[[i]], digits=2), sep, sep="")
  }
  return (name)
}

## Result's variables
Gdemo<- array(dim=c(nrow(strategy), 5, length(n0)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0))
GenvJ<- array(dim=c(nrow(strategy), 5, length(n0), nrow(envir)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, var.survJ=envir$var))
GenvA<- array(dim=c(nrow(strategy), 5, length(n0), nrow(envir)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, var.survA=envir$var))
GenvJA<- array(dim=c(nrow(strategy), 5, length(n0), nrow(envirJA)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, var.survJ_survA=paste(envirJA$varJ, envirJA$varA, sep="_"))
GseasonBestFirst<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason))
GseasonBestMean<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason))
GenvJ_seasonBestFirst<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season), nrow(envir)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason, var.survJ=envir$var))
GenvJ_seasonBestMean<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season), nrow(envir)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason, var.survJ=envir$var))
GenvA_seasonBestFirst<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season), nrow(envir)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason, var.survA=envir$var))
GenvA_seasonBestMean<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season), nrow(envir)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason, var.survA=envir$var))
GenvJA_seasonBestFirst<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season), nrow(envirJA)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason), var.survJ_survA=paste(envirJA$varJ, envirJA$varA, sep="_"))
GenvJA_seasonBestMean<- array(dim=c(nrow(strategy), 5, length(n0), nrow(season), nrow(envirJA)), dimnames=list(strategy=character(nrow(strategy)), descriptors=c("mean", "var", "Gmean", "replacement", "error"), n0=n0, amplSeason=season$amplSeason, var.survJ_survA=paste(envirJA$varJ, envirJA$varA, sep="_"))

## Simulation
# library(doMC) #no desa els resultats correctament
##TODO correct clutch rounding
for (n in 1:length(n0)){
  for (s in 1:nrow(strategy)){
    cat("n0:", n0[n], "\tLH:", s, "\t/", nrow(strategy), strategyName(strategy[s,]))
# Demographic stochasticity
    surv<- survdist(n0=n0[n], survA=strategy$sA[s], max.years=strategy$lifespan[e] * n0[n])
    ##FIXME clutch must be an integer or TODO a discrete probability distribution
    fert<- fertdist(surv, broods=1, clutch=strategy$b[s], survJ=strategy$sJ[s])
    R0<- fitnessdist(surv, fert, n0[n])

    Gdemo[s,1:3,n]<- t(sdistri(R0))
    Gdemo[s,4,n]<- 1 - cumsum(R0$probR0)[which(R0$R0 > 2)[1]]
    if (any(is.na(R0))) Gdemo[s,5,n]<- 1 else Gdemo[s,5,n]<- 0
    if (n == 1) dimnames(Gdemo)$strategy[s]<- strategyName(strategy[s,c(1,2,4,5,6)], sep=" ")
    cat("\tG=", Gdemo[s,3,n], "\n")

# Seasonal environment best first
    cat("\t seasonal fertility (best first)... amplitude/interval/broods: ")
    for (a in 1:nrow(season)){
#       cat(season$amplSeason[a], "/", season$breedInterval[a], "/", season$broods[a], ", ", sep="")
      ##FIXME clutch = fertility / broods must be an integer or a discrete probability distribution
      ##ERROR: plot((round(strategy$b/season$broods.x)*season$broods.x) ~ (round(strategy$b/season$broods.y)*season$broods.y))

      clutch<- rep(strategy$b[s] %/% season$broods[a], season$broods[a])
      if (strategy$b[s] %% season$broods[a] > 0){
	clutch[1:(strategy$b[s] %% season$broods[a])]<- clutch[1:(strategy$b[s] %% season$broods[a])] + 1
      }

      fertSeasonBestFirst<- fertdist(surv, broods=season$broods[a], clutch=clutch, survJ=strategy$sJ[s], meanSeason=season$meanSeason[a], amplSeason=season$amplSeason[a], breedInterval=season$breedInterval[a], alignCriterion="bestFirst")
      R0<- fitnessdist(surv, fertSeasonBestFirst, n0[n])

      GseasonBestFirst[s,1:3,n,a]<- t(sdistri(R0))
      GseasonBestFirst[s,4,n,a]<- 1 - cumsum(R0$probR0)[which(R0$R0 > 2)[1]]
      if (any(is.na(R0))) GseasonBestFirst[s,5,n]<- 1 else GseasonBestFirst[s,5,n]<- 0
      if (n == 1) dimnames(GseasonBestFirst)$strategy[s]<- strategyName(strategy[s,c(1,2,4,5,6)], sep=" ")
    }

# Seasonal environment best mean
    cat("DONE!\n\t seasonal fertility (best mean)... amplitude/interval/broods: ")
    for (a in 1:nrow(season)){
#       cat(season$amplSeason[a], "/", season$breedInterval[a], "/", season$broods[a], ", ", sep="")
      ##FIXME clutch must be an integer or a discrete probability distribution
      clutch<- rep(strategy$b[s] %/% season$broods[a], season$broods[a])
      if (strategy$b[s] %% season$broods[a] > 0){
	clutch[1:(strategy$b[s] %% season$broods[a])]<- clutch[1:(strategy$b[s] %% season$broods[a])] + 1
      }

      fertSeasonBestMean<- fertdist(surv, broods=season$broods[a], clutch=clutch, survJ=strategy$sJ[s], meanSeason=season$meanSeason[a], amplSeason=season$amplSeason[a], breedInterval=season$breedInterval[a], alignCriterion="bestMean")
      R0<- fitnessdist(surv, fertSeasonBestMean, n0[n])

      GseasonBestMean[s,1:3,n,a]<- t(sdistri(R0))
      GseasonBestMean[s,4,n,a]<- 1 - cumsum(R0$probR0)[which(R0$R0 > 2)[1]]
      if (any(is.na(R0))) GseasonBestMean[s,5,n]<- 1 else GseasonBestMean[s,5,n]<- 0
      if (n == 1) dimnames(GseasonBestMean)$strategy[s]<- strategyName(strategy[s,c(1,2,4,5,6)], sep=" ")
    }

# Environmental stochasticity
    cat("DONE!\n\t environmental stochastiticy on fertility... var: ")
    for (e in 1:nrow(envir)){
      cat(envir$var[e], ", ", sep="")
      if (sum(!is.na(fbeta(strategy$sJ[s], envir$var[e]))) == 2){
	##FIXME clutch must be an integer or a discrete probability distribution
	fertEnvJ<- fertdist(surv, broods=1, clutch=round(strategy$b[s]), survJ=strategy$sJ[s], var.survJ=envir$var[e])
	R0<- fitnessdist(surv, fertEnvJ, n0[n])

	GenvJ[s,1:3,n,e]<- t(sdistri(R0))
	GenvJ[s,4,n,e]<- 1 - cumsum(R0$probR0)[which(R0$R0 > 2)[1]]
	if (any(is.na(R0))) GenvJ[s,5,n]<- 1 else GenvJ[s,5,n]<- 0
      }else GenvJ[s,,n,e]<- rep(NaN, 5)
      if (n==1) dimnames(GenvJ)$strategy[s]<- strategyName(strategy[s,c(1,2,4,5,6)], sep=" ")
    }
    cat("DONE!\n")
  }
}

save.image(file=paste("cohortSimu", Sys.time(),".RData", sep=""))


## Postproces
################
## Arrays to data frames

# setwd("../../../../doctorat/projectes/LHT/R")
# load("cohortSimu2012-08-18 15:26:56.RData")
{
Gdemo<- as.data.frame(Gdemo)
Gdemo<- reshape(Gdemo, direction="long", varying=1:ncol(Gdemo))
names(Gdemo)[1]<- "n0"
GdemoLH<- merge(cbind(id=1:nrow(strategy), strategy), Gdemo)

dims<- dim(GenvJ)
##PERFER comprovar que no es desordena res
GenvJ<- data.frame(id=rep(1:dims[1], prod(dims[-c(1,2)])),
			n0=rep( rep(as.numeric(dimnames(GenvJ)$n0), each=dims[1]), dims[4]), 
			mean=as.vector(GenvJ[,1,,]),
			var=as.vector(GenvJ[,2,,]),
			G=as.vector(GenvJ[,3,,]),
			env.var=rep(envir$var, each=prod(dims[c(1,3)])))
GenvJLH<- merge(cbind(id=1:nrow(strategy), strategy), GenvJ)

dims<- dim(GseasonBestFirst)
GseasonBestFirst<- data.frame(id=rep(1:dims[1], prod(dims[-c(1,2)])),
			n0=rep( rep(as.numeric(dimnames(GseasonBestFirst)$n0), each=dims[1]), dims[4]), 
			mean=as.vector(GseasonBestFirst[,1,,]),
			var=as.vector(GseasonBestFirst[,2,,]),
			G=as.vector(GseasonBestFirst[,3,,]),
			amplSeason=rep(season$amplSeason, each=prod(dims[c(1,3)])),
			meanSeason=rep(season$meanSeason, each=prod(dims[c(1,3)])),
			breedInterval=rep(season$breedInterval, each=prod(dims[c(1,3)])),
			broods=rep(season$broods, each=prod(dims[c(1,3)])))
GseasonBestFirstLH<- merge(cbind(id=1:nrow(strategy), strategy), GseasonBestFirst)

dims<- dim(GseasonBestMean)
GseasonBestMean<- data.frame(id=rep(1:dims[1], prod(dims[-c(1,2)])),
			n0=rep( rep(as.numeric(dimnames(GseasonBestMean)$n0), each=dims[1]), dims[4]), 
			mean=as.vector(GseasonBestMean[,1,,]),
			var=as.vector(GseasonBestMean[,2,,]),
			G=as.vector(GseasonBestMean[,3,,]),
			amplSeason=rep(season$amplSeason, each=prod(dims[c(1,3)])),
			meanSeason=rep(season$meanSeason, each=prod(dims[c(1,3)])),
			breedInterval=rep(season$breedInterval, each=prod(dims[c(1,3)])),
			broods=rep(season$broods, each=prod(dims[c(1,3)])))
GseasonBestMeanLH<- merge(cbind(id=1:nrow(strategy), strategy), GseasonBestMean)


save(GdemoLH, file="Gdemo.R")
save(GseasonBestFirstLH, file="GseasonBestFirst.R")
save(GseasonBestMeanLH, file="GseasonBestMean.R")
save(GenvJLH, file="GenvJ.R")
}

#####################
## Explore results ##
#####################
setwd("../../../../doctorat/projectes/LHT/R")
load(file="Gdemo.R")
load(file="GenvJ.R")
load(file="GseasonBestMean.R")
load(file="GseasonBestFirst.R")
library(lattice)

1/ln(.9)
# postscript(file="cohortModel.ps")
pdf(file="cohortModel.pdf")
png(filename="cohortModel%02d.png")

## Gdemo
broodValue<- log10(1/(GdemoLH$broods - 1 / log(GdemoLH$sA))) ## SRC: http://www.countrysideinfo.co.uk/bird_lifespan.htm
xyplot(Gmean ~ sA | factor(n0), group=round(b), data=GdemoLH, main="Demographic stochasticity")
# xyplot(Gmean ~ sA | factor(round(b)), group=n0, data=GdemoLH)#, auto.key=TRUE)
# xyplot(Gmean ~ lambda | factor(round(b)), group=n0, data=GdemoLH)#, auto.key=TRUE)

## GenvJ
xyplot(G ~ sA | factor(n0), group=round(b), data=GenvJLH, main="Demographic + juvenile survival environmental stochasticity")
xyplot(G ~ sA | factor(n0) * env.var, group=round(b), data=GenvJLH, main="Demographic + juvenile survival environmental stochasticity")
# xyplot(G ~ sA | factor(round(b)), group=n0, data=GenvJLH)#, auto.key=TRUE)

## GseasonBestMean
broodValue<- log10(1/(GseasonBestMeanLH$broods - 1 / log(GseasonBestMeanLH$sA)))
# xyplot(G ~ sA | factor(n0), group=round(b), data=GseasonBestMeanLH)
xyplot(G ~ sA | factor(n0) * amplSeason * broods * breedInterval, group=round(b), data=GseasonBestMeanLH, main="Demographic stochasticity + seasonality (best mean)")
# xyplot(G ~ sA | factor(round(b)), group=n0, data=GseasonBestMeanLH)#, auto.key=TRUE)

## GseasonBestFirst
# xyplot(G ~ sA | factor(n0), group=round(b), data=GseasonBestFirstLH)
xyplot(G ~ sA | factor(n0) * amplSeason * broods * breedInterval, group=round(b), data=GseasonBestFirstLH, main="Demographic stochasticity + seasonality (best first)")
# xyplot(G ~ sA | factor(round(b)), group=n0, data=GseasonBestFirstLH)#, auto.key=TRUE)

dev.off()

## lambda Selection
################
rangLambda<- c(1,1.2)

x<- factor(1:5)
levels(x)<- paste("n0 =", 1:5)

x<- factor(paste("n0 =", GdemoLH$n0), levels=paste("n0 =", 2^(0:5)))

# alpha<- ?? ## TODO: es pot afegir al model?

## Gdemo
eggValue<- log10(1/(GdemoLH$b - 1 / log(GdemoLH$sA)))
sel<- GdemoLH$n0 > 0 & GdemoLH$lambda > min(rangLambda) & GdemoLH$lambda < max(rangLambda) & !is.na(GdemoLH$Gmean)
png(filename="../../../congressos/neobiota2012-Pontevedra/imatges/Gdemo%d.png",width=900, height=900, pointsize=24)
xyplot(Gmean ~ eggValue[sel] | factor(paste("n0 =", 2 * n0), levels=paste("n0 =", sort(unique(2 * n0)))), group=round(b), data=GdemoLH[sel,],
	main="Demographic stochasticity", sub="1 < lambda < 1.2", xlab="Egg value", ylab="G",
	par.settings=list(par.main.text=list(cex=4), par.sub.text=list(cex=2), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3), add.text=list(cex=3)) , scales=list(cex=3), lwd=2)
dev.off()
# xyplot(Gmean ~ sA | factor(round(b)), group=n0, data=GdemoLH)#, auto.key=TRUE)
# xyplot(Gmean ~ lambda | factor(round(b)), group=n0, data=GdemoLH)#, auto.key=TRUE)

## GenvJ
eggValue<- log10(1/(GenvJLH$b - 1 / log(GenvJLH$sA)))
sel<- GenvJLH$n0 > 2 & GenvJLH$lambda > min(rangLambda) & GenvJLH$lambda < max(rangLambda) & !is.na(GenvJLH$G)
png(filename="../../../congressos/neobiota2012-Pontevedra/imatges/GenvJLH%d.png",width=1300, height=1300, pointsize=24)
xyplot(G ~ eggValue[sel] | factor(paste("env. var =", env.var), levels=paste("env. var =", sort(unique(env.var)))) * factor(paste("n0 =", 2 * n0), levels=paste("n0 =", sort(unique(2 * n0)))), group=round(b), data=GenvJLH[sel,],
	main="Demographic + environmental stochasticity", xlab="Egg value", ylab="G",
	par.settings=list(par.main.text=list(cex=4), par.sub.text=list(cex=2), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3), add.text=list(cex=2)) , scales=list(cex=3), lwd=2)
dev.off()
# xyplot(G ~ sA | factor(round(b)), group=n0, data=GenvJLH)#, auto.key=TRUE)

## GseasonBestFirst
broodValue<- log10(1/(GseasonBestFirstLH$broods - 1 / log(GseasonBestFirstLH$sA)))
n0<- sort(unique(GseasonBestFirstLH$n0))
amplitude<- sort(unique(GseasonBestFirstLH$amplSeason)) ##
broods<- sort(unique(GseasonBestFirstLH$broods))
breedInterval<- sort(unique(GseasonBestFirstLH$breedInterval))

i<-6
sel<- GseasonBestFirstLH$n0 == n0[i] & GseasonBestFirstLH$lambda > min(rangLambda) & GseasonBestFirstLH$lambda < max(rangLambda) & !is.na(GseasonBestFirstLH$G)
# xyplot(G ~ sA | amplSeason * broods * breedInterval, group=round(b), data=GseasonBestFirstLH[sel,], main="Demographic stochasticity + seasonality (best first)")

# sel<- sel & GseasonBestFirstLH$breedInterval == breedInterval[2]
# xyplot(G ~ sA | amplSeason * breedInterval, group=broods, data=GseasonBestFirstLH[sel,], main="Demographic stochasticity + seasonality (best first)")#, auto.key=TRUE)
png(filename="../../../congressos/neobiota2012-Pontevedra/imatges/GseasonBestFirst%d.png",width=1300, height=1300, pointsize=24)
xyplot(G ~ broodValue[sel] | factor(paste("amplitude =", amplSeason), levels=paste("amplitude =", sort(unique(amplSeason)))) * factor(paste("broods =", broods), levels=paste("broods =", sort(unique(broods)))), group=breedInterval, data=GseasonBestFirstLH[sel,],
	main="Demographic stochasticity + seasonality", sub=paste("n0 =", 2 * n0[i]), xlab="Brood value", ylab="G",
	par.settings=list(par.main.text=list(cex=4), par.sub.text=list(cex=2), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3), add.text=list(cex=2)) , scales=list(cex=3), lwd=2)#, auto.key=TRUE)

dev.off()

# G ~ G
## TODO
plot((round(z$b/z$broods.x)*z$broods.x) ~ (round(z$b/z$broods.y)*z$broods.y))

x<- GseasonBestFirstLH[sel, c(1:3,5:6,8,11:15)]
x<- x[!(x$broods==1 & x$breedInterval >1),]
s<- unique(x[,8:11])
s<- s[order(s$amplSeason, s$broods, s$breedInterval),]

y<- list()
for (i in 1:nrow(s)){
  y[[i]]<- x[x$amplSeason == s$amplSeason[i] & x$broods == s$broods[i] & x$breedInterval == s$breedInterval[i],]
  names(y)[i]<- paste("amplitude=", s$amplSeason[i], " broods=", s$broods[i], " breedInterval=", s$breedInterval[i])
}

plot(1:55, 55:1, type="n")
for (i in 1:length(y)){
  for (j in 1:(i-1)){
  cat(i,j,"\n")
  points(i,j, col="red", pch=2)
#   z<- merge(y[[1]], y[[13]], by=c("id","lambda","sA","sJ","b","n0"))
  }
}

z<- merge(y[[1]], y[[11]], by=c("id","lambda","sA","sJ","b","n0")) # amplitude= 0.02 broods=
## plot((round(z$b/z$broods.x)*z$broods.x) ~ (round(z$b/z$broods.y)*z$broods.y)) #FIXME error when rounding clutch size
xyplot(G.x ~ G.y, group=factor(round(b)), data=z, xlab="G (broods=1)", ylab="G (broods=6)",
  panel= function(x,y, ...){
    panel.abline(0,1, lty=2, col="red", lwd=2)
    panel.xyplot(x, y, group=factor(round(z$b)), subscripts=TRUE)
  }, auto.key=list(space = "right", title="Clutch")
)
xyplot((z$G.x - z$G.y) ~ broodValue, group=(round(z$b/z$broods.x)*z$broods.x) - (round(z$b/z$broods.y)*z$broods.y), data=z, xlab="G (broods=1)", ylab="G (broods=6)", ylim=c(-4,1),
  auto.key=list(space = "right", title="Clutch")
)

z<- merge(y[[43]], y[[45]], by=c("id","lambda","sA","sJ","b","n0")) # 40-52, 53-65
xyplot(G.x ~ G.y, group=round(b), data=z, auto.key=TRUE)
abline(0, 1, col="red", lty=2, cex=2)

z<- merge(y[[34]], y[[44]], by=c("id","lambda","sA","sJ","b","n0")) # 40-52, 53-65
plot(z$G.x ~ z$G.y)
abline(0, 1, col="red", lty=2, cex=2)

z<- merge(y[[53]], y[[57]], by=c("id","lambda","sA","sJ","b","n0")) # 40-52, 53-65
plot(z$G.x ~ z$G.y)
abline(0, 1, col="red", lty=2, cex=2)

z<- merge(y[[53]], y[[57]], by=c("id","lambda","sA","sJ","b","n0")) # 40-52, 53-65
plot(z$G.x ~ z$G.y)
abline(0, 1, col="red", lty=2, cex=2)
# y<- reshape(x, v.names = "G", timevar = "id", idvar = c(8:11), direction="wide", new.row.names = NULL)

## GseasonBestMean
sel<- GseasonBestMeanLH$n0 == n0 & GseasonBestMeanLH$lambda > min(rangLambda) & GseasonBestMeanLH$lambda < max(rangLambda) & !is.na(GseasonBestMeanLH$G)

xyplot(G ~ sA | amplSeason * breedInterval, group=broods, data=GseasonBestMeanLH[sel,], main="Demographic stochasticity + seasonality (best first)")#, auto.key=TRUE)
xyplot(G ~ sA | amplSeason * broods, group=breedInterval, data=GseasonBestMeanLH[sel,], main="Demographic stochasticity + seasonality (best first)")#, auto.key=TRUE)

broodValue<- log10(1/(GseasonBestMeanLH$broods - 1 / log(GseasonBestMeanLH$sA)))
n0<- sort(unique(GseasonBestMeanLH$n0))
amplitude<- sort(unique(GseasonBestMeanLH$amplSeason)) ##
broods<- sort(unique(GseasonBestMeanLH$broods))
breedInterval<- sort(unique(GseasonBestMeanLH$breedInterval))

i<- 6
sel<- GseasonBestMeanLH$n0 == n0[i] & GseasonBestMeanLH$lambda > min(rangLambda) & GseasonBestMeanLH$lambda < max(rangLambda) & !is.na(GseasonBestMeanLH$G)
# xyplot(G ~ sA | amplSeason * broods * breedInterval, group=round(b), data=GseasonBestMeanLH[sel,], main="Demographic stochasticity + seasonality (best first)")

# sel<- sel & GseasonBestMeanLH$breedInterval == breedInterval[2]
# xyplot(G ~ sA | amplSeason * breedInterval, group=broods, data=GseasonBestMeanLH[sel,], main="Demographic stochasticity + seasonality (best first)")#, auto.key=TRUE)
png(filename="../../../congressos/neobiota2012-Pontevedra/imatges/GseasonBestMean%d.png",width=1300, height=1300, pointsize=24)
xyplot(G ~ broodValue[sel] | factor(paste("amplitude =", amplSeason), levels=paste("amplitude =", sort(unique(amplSeason)))) * factor(paste("broods =", broods), levels=paste("broods =", sort(unique(broods)))), group=breedInterval, data=GseasonBestMeanLH[sel,],
	main="Demographic stochasticity + seasonality", sub=paste("n0 =", 2 * n0[i]), xlab="Brood value", ylab="G",
	par.settings=list(par.main.text=list(cex=4), par.sub.text=list(cex=2), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3), add.text=list(cex=2)) , scales=list(cex=3), lwd=2)#, auto.key=TRUE)

dev.off()

#######################################################
x<- which(strategy$lambda > 1 & strategy$lambda < 1.1)

par(mfrow=c(2,2))
for (i in 1:4){
  plot(propietatsDistri[x,3,i] ~ strategy$sA[x], main=paste("n0=", dimnames(propietatsDistri)$n0[i]), ylab="G = mu - 2 var/mu", xlab="Sa")
}

library(rgl)
plot3d(propietatsDistri[,3,1], strategy$sA, strategy$lambda)

colSums(R0poblacio[["probabilitat"]], na.rm=TRUE) # la probabilitat total per cada n0 ha de ser ~ 1
R0poblacioAcumulada<- Pacumulada(R0poblacio[["probabilitat"]])

xlim<- c(0, R0poblacio[["R0"]][which(is.na(R0poblacio[["probabilitat"]][,1]))[1]-1,1])
matplot(R0poblacio[["R0"]], R0poblacio[["probabilitat"]], type="s", ylab="Probability", xlab="R0", xlim=xlim)
matplot(R0poblacio[["R0"]], R0poblacioAcumulada, type="s", ylab="Probability", xlab="R0", xlim=xlim)
abline(v=2, col="red", lwd=2)
Pestabliment(R0poblacio, minR0establiment=2)

# [http://www.statmethods.net/advgraphs/probability.html]


