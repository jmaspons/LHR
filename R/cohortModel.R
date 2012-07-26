## COHORT DEMOGRAPHIC MODEL
#############################################################
## P(years survived) = anys viscuts pel conjunt de la població adulta (de moment assumim sa constant) -> matriu de lefkovitch?
# This represents the number of failures which occur in a sequence of Bernoulli trials before a target number of successes is reached. [?rnbinom]
# anys de vida fertil * N0 = suma(rnbinom() 1:N0) => size=N0 Convolució de binomials negatives [http://en.wikipedia.org/wiki/List_of_convolutions_of_probability_distributions] [http://www.usna.edu/MathDept/.courses/pre97/sm230/sums.htm : Sums of negative binomials with the same p have a negative binomial distribution with number of successes equal to the total of each number of successes.]
# Si size=1 => rnbinom = rexp

survdist<- function(n0, survA, var.survA, maxPomited=0.01, max.years){
  if (missing(var.survA)){
    if (missing(max.years))
      max.years<- qnbinom(maxPomited, size=n0, prob=1-survA, lower.tail=FALSE)
    prob<- dnbinom(0:max.years, size=n0, prob=1-survA)
  }else{
    if (missing(max.years))
      max.years<- qbetanbinom(maxPomited, size=n0, mean=1-survA, variance=var.survA, lower.tail=FALSE)
    prob<- dbetanbinom(0:max.years, size=n0, mean=1-survA, variance=var.surbA)
  }
  prob<- prob / sum(prob)
  ans<- data.frame(years=0:max.years, probS=prob)

  return (ans)
}
# s<- survdist(5, 0.6)

## P(R0) / P(years survived) = number of offspring that reach independence for each number of years lived
##TODO: add temporal autocorrelated environment
fertdist<- function(years, broods, clutch, survJ, var.survJ, meanSeason, amplSeason, breedInterval){
  ans<- list()

  for (i in years$year){
    size<- i * broods * clutch
    if (missing(var.survJ)){
      if(missing(meanSeason) & missing(amplSeason) & missing(breedInterval))
	prob<- dbinom(0:size, size=size, prob=survJ)
      else{
	seasons<- seasonality(i, meanSeason, amplSeason)
	seasons<- par.seasonality(seasons, broods, breedInterval, criterion="bestMean")
	prob<- dbinom(rep(0:size, each=broods), size=size, prob=survJ * seasons)
	prob<- matrix(prob, ncol=broods, byrow=TRUE)
	prob<- rowMeans(prob)
      }
    }else{
      if(missing(meanSeason) & missing(amplSeason) & missing(breedInterval)){
	parBeta<- TrobaParametresBeta(survJ, var.survJ)
	prob<- dbetabinom(0:size, size=size, parBeta[[1]], parBeta[[2]])
      }else{
	seasons<- seasonality(i, meanSeason, amplSeason)
	seasons<- par.seasonality(seasons, broods, breedInterval, criterion="bestMean")
	parBeta<- TrobaParametresBeta(survJ * seasons, var.survJ)
	prob<- dbinom(rep(0:size, each=broods), size=size, parBeta[[1]], parBeta[[2]])
	prob<- matrix(prob, ncol=broods, byrow=TRUE)
	prob<- rowMeans(prob)
      }
    }
    tmp<-  data.frame(fert=0:size, probF=prob, years=i)
    ans[[i+1]]<- tmp
  }
  ans<- as.data.frame(do.call("rbind", ans)) #passar list[[n=data.frame[1]]] a data.frame(mean, variance, espG)

  return (ans)
}
# f<- fertdist(s, 2, 5, .25)

## P(R0|survA)
# mirar paquet distr: Differences between random variables and distribution algebra?
# file:///usr/local/lib/R/site-library/distrDoc/doc/distr.pdf
fitnessdist<- function(surv, fert, n0){
  ans<- merge(surv, fert, by="years")
  ans<- data.frame(ans, probR0=numeric(nrow(ans)))
  ans$probR0<- ans$probS * ans$probF
  ans<- by(ans$probR0, ans$fert, sum)
  ans<- data.frame(R0=as.numeric(names(ans)), probR0=as.numeric(ans))
  if (!missing(n0)) ans$R0<- ans$R0 / n0

  return (ans)
}
# R0<- fitnessdist(s,f, 5)

## All in one
# Wrapper to run the complete model in on call
fitnessdistAIO<- function(n0, survA, var.survA, broods, clutch, survJ, var.survJ, maxPomited, max.years){
  parameters<- list(n0=n0, survA=survA)
  n<-2
  if (!missing(var.survA)){
    parameters[[n + 1]]<- var.survA
    names(parameters)[n+1]<- "var.survA"
    n<- n+1
  }
  if (!missing(maxPomited)){
    parameters[[n + 1]]<- maxPomited
    names(parameters)[n+1]<- "maxPomited"
    n<- n+1
  }
  if (!missing(max.years)){
    parameters[[n + 1]]<- max.years
    names(parameters)[n+1]<- "max.years"
    n<- n+1
  }

  surv<- do.call(survdist, args=parameters)

  parameters<- list(years=surv, broods=broods, clutch=clutch, survJ=survJ)
  n<-4
  if (!missing(var.survJ)){
    parameters[[n + 1]]<- var.survJ
    names(parameters)[n+1]<- "var.survJ"
    n<- n+1
  }

  fert<- do.call(fertdist, args=parameters)

  return (fitnessdist(surv, fert, n0))
}
# R02<- fitnessdistAIO(5, 0.6, broods=2, clutch=5, survJ=0.25)
# identical(R0, R02)

# plot(cumsum(R0$probR0) ~ R0$R0)

Pacumulada<- function(probabilitat){ # usar cumsum(R0$probR0)
  probabilitatAcumulada<- probabilitat # Densitat de la distribució de probabilitat (acumulada)
  for (i in 2:nrow(probabilitat)){
    probabilitatAcumulada[i,]<- colSums(probabilitat[1:i,], na.rm=TRUE)
  }
  return (probabilitatAcumulada)
}

Pestabliment<- function(R0poblacio, minR0establiment=2){
  resultat<- numeric(ncol(R0poblacio[["R0"]])) # P(establiment) ==> R0 => 2 (remplaç de la població)
  names(resultat)<- attributes(R0poblacio[["R0"]])$dimnames$n0
  for (i in 1:length(resultat)){
    resultat[i]<- which.min(abs(R0poblacio[["R0"]][,i] - minR0establiment)) # posició de l'R0 més proper a 2
    if (R0poblacio[["R0"]][resultat[i],i] < minR0establiment){ # si el més proper està per sota incrementa en una posició
      resultat[i]<- resultat[i] + 1
    }
    resultat[i]<- 1 - sum(R0poblacio[["probabilitat"]][1:resultat[i],i])
  }
  return (resultat)
}

PropietatsDistribucio<- function(R0poblacio){
  resultat<- data.frame(mu=numeric(ncol(R0poblacio[["R0"]])), var=numeric(ncol(R0poblacio[["R0"]])), G=numeric(ncol(R0poblacio[["R0"]])), row.names=attributes(R0poblacio[["R0"]])$dimnames$n0) # P(establiment) ==> R0 => 2 (remplaç de la població)
  for (i in 1:nrow(resultat)){
    x<- which(!is.na(R0poblacio[["probabilitat"]][,i]))
    resultat$mu[i]<- weighted.mean(R0poblacio[["R0"]][x,i], R0poblacio[["probabilitat"]][x,i], na.rm=TRUE)
    resultat$var[i]<- sum(R0poblacio[["probabilitat"]][x,i] * (R0poblacio[["R0"]][x,i] - resultat$mu[i])^2)
    resultat$G[i]<- G(resultat$mu[i], resultat$var[i])
  }
  return (resultat)
}

MitjanaDistribucio<- function(R0poblacio){ # Mitjana ponderada per la probabilitat
  mitjana<- numeric(ncol(R0poblacio[["R0"]])) # P(establiment) ==> R0 => 2 (remplaç de la població)
  names(mitjana)<- attributes(R0poblacio[["R0"]])$dimnames$n0
  for (i in 1:length(mitjana)){
    x<- which(!is.na(R0poblacio[["probabilitat"]][,i]))
    mitjana[i]<- weighted.mean(R0poblacio[["R0"]][x,i], R0poblacio[["probabilitat"]][x,i], na.rm=TRUE)
  }
  return (mitjana)
}

VariancaDistribucio<- function(R0poblacio){ # http://en.wikipedia.org/wiki/Variance#Discrete_random_variable
  mu<- MitjanaDistribucio(R0poblacio)
  varianca<- numeric(ncol(R0poblacio[["R0"]]))
  names(varianca)<- attributes(R0poblacio[["R0"]])$dimnames$n0
  for (i in 1:length(varianca)){
    x<- which(!is.na(R0poblacio[["probabilitat"]][,i]))
    varianca[i]<- sum(R0poblacio[["probabilitat"]][x,i] * (R0poblacio[["R0"]][x,i] - mu[i])^2)
    varianca[i]<- varianca[i] / sum(R0poblacio[["probabilitat"]][,i], na.rm=TRUE) # Correcció de la probabilitat omesa
  }
  return (varianca)
}

G<- function(mu, v){
  return (mu - 2 * v / mu)
}

espG<- function(x){  
  x.mean<- mean(x)
  x.variance<- var(x)
  G.esp<- x.mean - 2 * x.variance / x.mean
  return (data.frame(mean=x.mean, variance=x.variance, G.esp=G.esp))
}

# 
# # [http://www.statmethods.net/advgraphs/probability.html]
# 
# 
## GRÀFIC
##########
# library(lattice)
## n(sa, n0) = anys viscuts pel conjunt de la població adulta (de moment assumim sa constant) -> matriu de lefkovitch?
# x<- PAnysViscuts(n0=n0, sa=.7)
# xyplot(PnAnys ~ nAnys, group=n0, data=x, pch=16, type="s")

## Fitness(n, exit de cada posta) = n * b * sea -> nombre d'adults produits al llarg de la vida per n0
# x<-20 #perquè el gràfic quedi més bonic
# seqExplor<- trunc(seq(1, max(resulN$n), length=7)) # només grafica per uns quants n
# n<-0 #perquè el gràfic quedi més bonic
# i<-1 #perquè el gràfic quedi més bonic
# curve(dbinom(x, size=n[i]*b, prob=sea), from=0, to=x, n=x+1, type="n", add=FALSE, col=0, ylab=)
# for (i in 1:length(seqExplor)){
#   curve(dbinom(x, size=seqExplor[i] * nB, prob=sea), from=0, to=x, n=x+1, type="b", pch=1, add=TRUE, col=topo.colors(length(seqExplor))[i])
# }

## NUMÈRIC
###########
## n(sa, n0) = anys viscuts pel conjunt de la població adulta (de moment assumim sa constant) -> matriu de lefkovitch?
# n<-numeric() 
# for (i in 1:length(n0)){
#   n<- c(n, rnbinom(nRep, size=n0[i], prob=1-sa))
# }
# 
# #Preparació de la taula pels resultats numèrics
# resulN<- data.frame(n0=rep(n0, each=nRep), n=n, R0=numeric(length(n0)*nRep))
# 
# ## Fitness(n, exit de cada posta) = n * b * sea -> nombre d'adults produits al llarg de la vida per n0
# for (i in 1:length(n0)){
#   R0<- rbinom(nRep, size=resulN$n[which(resulN$n0 == n0[i])] * b, prob=sea)
#   resulN$R0[((i-1) * nRep +1):(i*nRep)]<- R0 / resulN$n0[((i-1) * nRep +1):(i*nRep)] 
# }
# plot(resulN)
# 
# ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
# 
# #################################################
# ## Prova efecte de mostreig / disperció del risc
# 
# ## Funcions
# G<- function(mu, v){
#   return (mu - 2 * v / mu)
# }
# 
