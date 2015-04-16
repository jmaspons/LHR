## IMPLEMENTED IN C: USE exploreCohortModel() INSTEAD
# simulation<- data.frame(n0, survA, var.survA, broods, B, survJ, var.survJ, meanSeason, amplSeason, breedInterval, alignCriterion=c("bestFirst", "bestMean"))
runCohortModel<- function(simulation, maxPomited=0.05, verbose=FALSE){
  ans<- data.frame(mean=numeric(nrow(simulation)), var=numeric(nrow(simulation)), G=numeric(nrow(simulation)), Preplacement=numeric(nrow(simulation)), error=logical(nrow(simulation)))
  surv<- unique(simulation[,1:3]) # reuse survival distribution for all the simulations with the same parameters  (n0, survA, var.survA)
  
  reg<- which(simulation$n0 == surv$n0[1] & simulation$survA == surv$survA[1] & simulation$var.survA == surv$var.survA[1])
  cat(1, "/", nrow(surv), "set of survival parameters. n0=", surv$n0[1], "\tsA=", surv$survA[1], "\tvar.sA", surv$var.survA[1], "\t|", length(reg), "sets of fertility parameters\n")
  for (i in 1:nrow(surv)){
    reg<- which(simulation$n0 == surv$n0[i] & simulation$survA == surv$survA[i] & simulation$var.survA == surv$var.survA[i])
#     if (i %% 25 == 0 || verbose) 
{cat(i, "/", nrow(surv), "surv parameters. n0=", surv$n0[i], "\tsA=", surv$survA[i], "\tvar.sA", surv$var.survA[i], "\t|", length(reg), "sets of fertility parameters\n")}
    
    if (surv$var.survA[i] == 0){
      years<- survdist(n0=surv$n0[i], survA=surv$survA[i], maxPomited=maxPomited)
    }else{
      if (sum(is.na(fbeta(surv$survA[i], surv$var.survA[i]))) == 2){
	ans[reg,]<- NA
	next
      }
      years<- survdist(n0=surv$n0[i], survA=surv$survA[i], var.survA=surv$var.survA[i], maxPomited=maxPomited)
    }
    
    jj<- 0
    for (j in reg){
      if (verbose) {jj<- jj + 1; cat("\t", jj, "/", length(reg), "fert parameters.\tbroods=", simulation$broods[j], "\tB=", simulation$B[j], "\tsurvJ=", simulation$survJ[j], "\tvar.survJ=", simulation$var.survJ[j], "\tamplSeason=", simulation$amplSeason[j], "\tbreedInterval=", simulation$breedInterval[j], "\talignCriterion =", simulation$alignCriterion[j], "\n")} #& jj %% 10 == 0
      if (simulation$var.survJ[j] == 0){
	fert<- fertdist(years, broods=simulation$broods[j], B=simulation$B[j], survJ=simulation$survJ[j], meanSeason=simulation$meanSeason[j], amplSeason=simulation$amplSeason[j], breedInterval=simulation$breedInterval[j], alignCriterion=simulation$alignCriterion[j])
      }else{
	if (sum(is.na(fbeta(simulation$survJ[j], simulation$var.survJ[j]))) == 2){
	  ans[reg,]<- NA
	  next
	}
	fert<- fertdist(years, broods=simulation$broods[j], B=simulation$B[j], survJ=simulation$survJ[j], var.survJ=simulation$var.survJ[j], meanSeason=simulation$meanSeason[j], amplSeason=simulation$amplSeason[j], breedInterval=simulation$breedInterval[j], alignCriterion=simulation$alignCriterion[j])
      }
      
      R0<- fitnessdist(years, fert, simulation$n0[j])
      ans[j,1:3]<- sdistri(R0)
      ans$Preplacement[j]<- 1 - cumsum(R0$probR0)[which(R0$R0 > 2)[1]]
      if (any(is.na(R0))) ans$error[j]<- TRUE else ans$error[j]<- FALSE
    }
  }
  
  return (ans)
}