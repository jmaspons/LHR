## COHORT DEMOGRAPHIC MODEL
############################
cohortModel<- function(n0, survA, varSurvA, limit, broods, B, survJ, varSurvJ, season){
#   TODO check inputs
  if (length(season) != broods) stop()
  .External("cohortModel", as.integer(n0), as.numeric(survA), as.numeric(var.survA), as.numeric(limit), as.integer(broods), as.integer(B), as.numeric(survJ), as.numeric(var.survJ), as.numeric(seasonAmpl), as.numeric(breedInterval))
}

exploreCohortModel<- function(simulation){
#   TODO check inputs
  result<- with(simulation, .External("exploreCohortModel", as.integer(n0), as.numeric(survA), as.numeric(var.survA), as.numeric(limit), as.integer(broods), as.integer(B), as.numeric(survJ), as.numeric(var.survJ), as.numeric(amplSeason), as.numeric(breedInterval)))
  result<- as.data.frame(result)
  colnames(result)<- c("mean", "var", "G", "Preplace", "error")

  return (result)
}

## P(years survived) = anys viscuts pel conjunt de la població adulta (de moment assumim sa constant) -> matriu de lefkovitch?
# This represents the number of failures which occur in a sequence of Bernoulli trials before a target number of successes is reached. [?rnbinom]
# anys de vida fertil * N0 = suma(rnbinom() 1:N0) => size=N0 Convolució de binomials negatives [http://en.wikipedia.org/wiki/List_of_convolutions_of_probability_distributions] [http://www.usna.edu/MathDept/.courses/pre97/sm230/sums.htm : Sums of negative binomials with the same p have a negative binomial distribution with number of successes equal to the total of each number of successes.]
# Si size=1 => rnbinom = rexp

survdist<- function(n0, survA, var.survA, maxPomitted=0.01, max.years){
  if (missing(var.survA)){
    if (missing(max.years))
      max.years<- qnbinom(maxPomitted, size=n0, prob=survA, lower.tail=FALSE)
    prob<- dnbinom(0:max.years, size=n0, prob=1-survA)
  }else{ # var.survA
    parBeta<- fbeta(survA, var.survA)
    if (missing(max.years)){
      max.years<- qbetanbinom(maxPomitted, size=n0, parBeta[[1]], parBeta[[2]], lower.tail=FALSE)
    }
    prob<- exp(dbetanbinom(0:max.years, size=n0, parBeta[[1]], parBeta[[2]], log=TRUE))
  }
  prob<- prob / sum(prob, rm.na=TRUE) ## Correct for P omitted or lifespan. Sum(prob) = 1
  ans<- data.frame(years=0:max.years, probS=prob)

  return (ans)
}
# s<- survdist(5, 0.6)

## P(R0) / P(years survived) = number of offspring that reach independence for each number of years lived
##TODO: add temporal autocorrelated environment
fertdist<- function(years, broods=1, clutch, survJ, B, var.survJ, seasonalPattern, meanSeason, amplSeason, breedInterval, alignCriterion="bestFirst") ##Clutch size is a int vector of length broods
{
  if (!missing(clutch) & !missing(B)) warning("Redundant input: especify clutch or B")
  if (missing(B)){
    B<- clutch * broods
  }
  ans<- list()
  
  for (i in years$year){
    size<- i * B
     
    if (meanSeason == 1 & amplSeason == 0 | missing(seasonalPattern) & missing(meanSeason) & missing(amplSeason) & missing(breedInterval)){
      if (missing(var.survJ)){
	      prob<- dbinom(0:size, size=size, prob=survJ)
      }else{
      	parBeta<- fbeta(survJ, var.survJ)
      	prob<- exp(dbetabinom(0:size, size=size, parBeta[[1]], parBeta[[2]], log=TRUE))
      }

    }else{ # seasonality
      if (missing(clutch)){
	      clutch<- rep(B %/% broods, broods)
      	mod<- B %% broods
      	
      	if (mod > 0){
      	  clutch[1:mod]<- clutch[1:mod] + 1
      	}
      }

      size<- i * clutch

      if (missing(seasonalPattern)){
	      seasonalPattern<- par.seasonality(broods=broods, breedInterval=breedInterval, mean=meanSeason, ampl=amplSeason, criterion=alignCriterion)
      }
      
      prob<- matrix(0, nrow=max(size) + 1, ncol=broods, dimnames=list(0:max(size), NULL))
      
      if (missing(var.survJ)){
      	for (j in 1:broods){
      	  prob[1:(size[j] + 1),j]<- dbinom(0:size[j], size=size[j], prob=survJ * seasonalPattern[j])
      	}
      }else{
      	parBeta<- fbeta(survJ * seasonalPattern, var.survJ)
      
      	for (j in 1:broods){
      	  prob[1:(size[j] + 1),j]<- exp(dbetabinom(0:size[j], size=size[j], parBeta[[1]][j], parBeta[[2]][j], log=TRUE))
      	}
      }
      
      prob<- rowMeans(prob)
    }
    tmp<-  data.frame(fert=0:max(size), probF=prob, years=i)
    ans[[i+1]]<- tmp
  }
  ans<- as.data.frame(do.call("rbind", ans)) #passar list[[n=data.frame[1]]] a data.frame(mean, variance, espG)

  return (ans)
}
# f<- fertdist(s, 2, 5, .25)

## P(R0|survA)
# mirar paquet distr: Differences between random variables and distribution algebra?
# file:///usr/local/lib/R/site-library/distrDoc/doc/distr.pdf
fitnessdist<- function(surv, fert, n0)
{
  ans<- merge(surv, fert, by="years")
  ans<- data.frame(ans, probR0=numeric(nrow(ans)))
  ans$probR0<- ans$probS * ans$probF

  if (sum(is.na(ans$probR0))){
    ans<- ans[!is.na(ans$probS),]
    error<- TRUE
    warning("omitted some years from 'survdist()' because probabilities are NaN. omitted probability will be larger than specified.")
  }else error<- FALSE

  ans<- by(ans$probR0, ans$fert, sum)
  ans<- data.frame(R0=as.numeric(names(ans)), probR0=as.numeric(ans))
  if (!missing(n0)) ans$R0<- ans$R0 / n0
  if (error) ans[nrow(ans) + 1,]<- c(NA, NA)

  return (ans)
}
# R0<- fitnessdist(s,f, 5)

## All in one
# Wrapper to run the complete model in on call
fitnessdistAIO<- function(n0, survA, var.survA, broods, clutch, survJ, var.survJ, maxPomitted, max.years)
{
  parametersSurv<- list(n0=n0, survA=survA)
  #Check for missing parameters
  parametersFert<- list(broods=broods, clutch=clutch, survJ=survJ) 

  n<-2
  if (!missing(var.survA)){
    parametersSurv[[n + 1]]<- var.survA
    names(parametersSurv)[n+1]<- "var.survA"
    n<- n+1
  }
  if (!missing(maxPomitted)){
    parametersSurv[[n + 1]]<- maxPomitted
    names(parametersSurv)[n+1]<- "maxPomitted"
    n<- n+1
  }
  if (!missing(max.years)){
    parametersSurv[[n + 1]]<- max.years
    names(parametersSurv)[n+1]<- "max.years"
    n<- n+1
  }

  surv<- do.call(survdist, args=parametersSurv)
  
  parametersFert<- list(years=surv, broods=broods, clutch=clutch, survJ=survJ)

  n<-4
  if (!missing(var.survJ)){
    parametersFert[[n + 1]]<- var.survJ
    names(parametersFert)[n+1]<- "var.survJ"
    n<- n+1
  }

  fert<- do.call(fertdist, args=parametersFert)

  return (fitnessdist(surv, fert, n0))
}
# R02<- fitnessdistAIO(5, 0.6, broods=2, clutch=5, survJ=0.25)
# identical(R0, R02)

# plot(cumsum(R0$probR0) ~ R0$R0)

# Pacumulada<- function(probabilitat){ # usar cumsum(R0$probR0)
#   probabilitatAcumulada<- probabilitat # Densitat de la distribució de probabilitat (acumulada)
#   for (i in 2:nrow(probabilitat)){
#     probabilitatAcumulada[i,]<- colSums(probabilitat[1:i,], na.rm=TRUE)
#   }
#   return (probabilitatAcumulada)
# }

# Pestabliment<- function(R0poblacio, minR0establiment=2){
#   resultat<- numeric(ncol(R0poblacio[["R0"]])) # P(establiment) ==> R0 => 2 (remplaç de la població)
#   names(resultat)<- attributes(R0poblacio[["R0"]])$dimnames$n0
#   for (i in 1:length(resultat)){
#     resultat[i]<- which.min(abs(R0poblacio[["R0"]][,i] - minR0establiment)) # posició de l'R0 més proper a 2
#     if (R0poblacio[["R0"]][resultat[i],i] < minR0establiment){ # si el més proper està per sota incrementa en una posició
#       resultat[i]<- resultat[i] + 1
#     }
#     resultat[i]<- 1 - sum(R0poblacio[["probabilitat"]][1:resultat[i],i])
#   }
#   return (resultat)
# }

sdistri<- function(distri){
  mean<- weighted.mean(distri[[1]], distri[[2]], na.rm=TRUE)
  var<- sum(distri[[2]] * (distri[[1]] - mean)^2)
  Gmean<- G(mean, var)

  return (data.frame(mean, var, Gmean))
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

