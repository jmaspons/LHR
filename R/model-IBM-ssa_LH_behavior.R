## Propensity functions for each transition ----
rateFunc.LH_Beh<- function(x, params, t){
  with(params, {
    return(c(repro1s=b * Pb1 * (1-PbF1) * x["N1s"], repro1sF=b * Pb1 * PbF1 * x["N1s"], 
             repro1b=b * Pb1 * (1-PbF1) * x["N1b"], repro1bF=b * Pb1 * PbF1 * x["N1b"],
             repro2s=b * Pb2 * (1-PbF2) * x["N2s"], repro2sF=b * Pb2 * PbF2 * x["N2s"], 
             repro2b=b * Pb2 * (1-PbF2) * x["N2b"], repro2bF=b * Pb2 * PbF2 * x["N2b"],
             skipR1=b * (1-Pb1) * x["N1b"], skipR1F=b * (1-Pb1) * x["N1bF"], skipR2=b * (1-Pb2) * x["N2b"], skipR2F=b * (1-Pb2) * x["N2bF"],
             dead1s=d1 * (1 + (x["N1s"] + x["N1b"] + x["N1bF"]) / K) * x["N1s"], dead1b=db1 * (1  + (x["N1s"] + x["N1b"] + x["N1bF"]) / K) * x["N1b"],
             dead1bF=db1 * (1 + (x["N1s"] + x["N1b"] + x["N1bF"]) / K) * x["N1bF"], dead1j=dj1 * (1  + (x["N1s"] + x["N1b"] + x["N1bF"] + x["N1j"]) / K) * x["N1j"],
             dead2s=d2 * (1 + (x["N2s"] + x["N2b"] + x["N2bF"]) / K) * x["N2s"], dead2b=db2 * (1 + (x["N2s"] + x["N2b"] + x["N2bF"]) / K) * x["N2b"],
             dead2bF=db2 * (1 + (x["N2s"] + x["N2b"] + x["N2bF"]) / K) * x["N2bF"], dead2j=dj2 * (1  + (x["N2s"] + x["N2b"] + x["N2bF"] + x["N2j"]) / K) * x["N2j"],
             move12s=c1 * (1-P1s) * x["N1s"], move21s=c2 * P1s * x["N2s"], move12b=c1 * (1-P1b) * x["N1b"], move21b=c2 * P1b * x["N2b"],
             move12bF=cF * (1-P1b) * x["N1bF"], move21bF=cF * P1b * x["N2bF"], move12j=c1 * (1-P1j) * x["N1j"], move21j=c2 * P1j * x["N2j"],
             grow1=g1 * x["N1j"], grow2=g2 * x["N2j"])
    )}
  )
}

rateFunc.LH_Beh<- compiler::cmpfun(rateFunc.LH_Beh) # byte-compile the function


# State-change matrix for each transition ----
## Notes: individuals from NxbF move to N-xb. Individuals which failed on the last reproductive event have
# higher probability to move to another habitat (c < cF) but they relax after an habitat change to evaluate the reproduction outcome in the new habitat.
transitionMat.LH_Beh<- function(params=data.frame(clutch1=1, clutch2=1)){
  transMat<- with(params, expr={
  #             repro1      repro2                    skip          dead1       dead2         move                       grow
    matrix(c(-1,     -1, 0, 0, 0, 0, 0, 0,   1, 1, 0, 0,  -1, 0, 0, 0, 0, 0, 0, 0,  -1, 1, 0, 0, 0, 0, 0, 0,    1, 0, # skip                #habitat 1
              1,       0,       0,-1, 0, 0, 0, 0,  -1, 0, 0, 0,   0,-1, 0, 0, 0, 0, 0, 0,   0, 0,-1, 1, 0, 1, 0, 0,    0, 0, # breeding
              0,       1,       0, 1, 0, 0, 0, 0,   0,-1, 0, 0,   0, 0,-1, 0, 0, 0, 0, 0,   0, 0, 0, 0,-1, 0, 0, 0,    0, 0, # Breeding Fail
              clutch1, 0, clutch1, 0, 0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0,-1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,-1, 1,   -1, 0, # juveniles
              #habitat 2                
              0, 0, 0, 0,      -1,-1,       0, 0,   0, 0, 1, 1,   0, 0, 0, 0,-1, 0, 0, 0,   1,-1, 0, 0, 0, 0, 0, 0,    0, 1, # skip
              0, 0, 0, 0,       1, 0,       0,-1,   0, 0,-1, 0,   0, 0, 0, 0, 0,-1, 0, 0,   0, 0, 1,-1, 1, 0, 0, 0,    0, 0, # breeding
              0, 0, 0, 0,       0, 1,       0, 1,   0, 0, 0,-1,   0, 0, 0, 0, 0, 0,-1, 0,   0, 0, 0, 0, 0,-1, 0, 0,    0, 0, # Breeding Fail
              0, 0, 0, 0, clutch2, 0, clutch2, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,-1,   0, 0, 0, 0, 0, 0, 1,-1,    0,-1), # juveniles
            ncol=30, byrow=TRUE, 
            dimnames=list(state=c("N1s", "N1b", "N1bF", "N1j", "N2s", "N2b", "N2bF", "N2j"), 
                          event=c("repro1s", "repro1sF", "repro1b","repro1bF",  "repro2s","repro2sF", "repro2b","repro2bF",
                                  "skipR1","skipR1F","skipR2","skipR2F",   "dead1s","dead1b","dead1bF","dead1j",  "dead2s","dead2b","dead2bF","dead2j",
                                  "move12s","move21s","move12b","move21b","move12bF","move21bF","move12j","move21j",  "grow1", "grow2")))
  })
  return (transMat)
}

# Parameters ----
# returns a different strategies and scenarios
# diffHab2: named vector with the differences in the parameters at habitat 2 respect habitat 1
# Warning: clutch have no effect on the simulation. It's necessary to modify the reaction channels (getReactionChannels(clutch1, clutch2))
getParams.LH_Beh.ssa<- function(strategy=c("slow", "fast", "freqRepro"), diffHab2, habDiffScenario="identicalHab", behavior="neutral"){
  params<- lapply(strategy, function(x){
    tmp<- switch(x,
                 slow=c(clutch1=1, clutch2=1,   b=1, PbF1=.4, PbF2=.4,  d1=.1,db1=.25,dj1=.25,  d2=.1,db2=.25,dj2=.25, g1=1, g2=1, K=500),
                 fast=c(clutch1=8, clutch2=8,   b=1, PbF1=.4, PbF2=.4,  d1=.5,db1=.8,dj1=1.5, d2=.5,db2=.8,dj2=1.5, g1=1, g2=1, K=500),
                 freqRepro=c(clutch1=2, clutch2=2,   b=4, PbF1=.4, PbF2=.4,  d1=.5,db1=.8,dj1=.8,  d2=.5,db2=.8,dj2=.8, g1=1, g2=1, K=500))
    tmp<- c(tmp, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5) # add neutral behavior
  })
  names(params)<- strategy
  
  if (!missing(diffHab2)){
    params<- lapply(params, function(x) setParams2diff1(x, diffHab2, type="multiplicative"))
  }else{
    params<- lapply(params, function(x) setScenario.ssa(x, habDiffScenario))
    names(params)<- paste(names(params), habDiffScenario, sep="-")
    params<- lapply(params, function(x) setBehavior(x, behavior))
    names(params)<- paste(names(params), behavior, sep="_")
  }
  params<- data.frame(do.call("rbind", params), stringsAsFactors=FALSE)
  params<- data.frame(idScenario=rownames(params), params, stringsAsFactors=FALSE)
  
  return (params)
}

# # Return new parames where parameters on habitat 2 are modified according to parameters from habitat 1 and a difference.
# # diff: named vector with the proportion of difference respect habitat 1.
# setParams2diff1<- function(params, diff=c(clutchDiff=0, PbDiff=0, PbFDiff=0, dDiff=0, dbDiff=0, djDiff=0, gDiff=0)){
#   selHab1<- sort(sapply(gsub("Diff$", "", names(diff)), function (x){
#     grep(paste0(x, "[1]{1}$"), names(params), value=TRUE)
#   }))
#   selHab2<- sort(sapply(gsub("Diff$", "", names(diff)), function (x){
#     grep(paste0(x, "[2]{1}$"), names(params), value=TRUE)
#   }))
#   diff<- diff[order(names(diff))]
#   for (i in 1:length(diff)){
#     params[selHab2[i]]<- params[selHab1[i]] * (1+diff[i])
#   }
#   
#   return (params)
# }


# params<- getParams()
# setParams2diff1(params, diff=c(PbFDiff=.5, dDiff=-.3, gDiff=-.2))

getScenario.ssa<- function(habDiffScenario){
  diff<- switch(habDiffScenario,
                `identicalHab`=c(clutchDiff=0, PbFDiff=0, dDiff=0, dbDiff=0, djDiff=0, gDiff=0),
                `mortalHab2`=c(clutchDiff=0, PbFDiff=0, dDiff=1, dbDiff=1, djDiff=1, gDiff=0),
                `nestPredHab2`=c(clutchDiff=0, PbFDiff=1, dDiff=0, dbDiff=0, djDiff=0, gDiff=0)
  )
  return (diff)
}

setScenario.ssa<- function(params, habDiffScenario="identicalHabs"){
  params<- setParams2diff1(params=params, diff=getScenario.ssa(habDiffScenario))

  return (params)
}

# setBehavior<- function(params, behavior){
#   if (any(grepl("neutral", behavior))){
#     params[c("Pb1","Pb2",  "c1", "c2", "cF",  "P1s","P1b","P1j")]<- c(1,1, 1,1,1, .5,.5,.5)
#   }
#   
#   if (any(grepl("skip", behavior))){ ## Increase breeding costs or better habitat selection
#     params[c("Pb1", "Pb2")]<- c(1, .2)
#   }
#   
#   if (any(grepl("learnBreed", behavior))){
#     params[c("c1", "c2", "cF")]<- c(0, 0, 3)
#   }
#   
#   if (any(grepl("learnExploreBreed", behavior))){ # Avoid habitat 2 after exploring or breeding fail
#     params[c("c1", "c2", "cF")]<- c(1, 3, 3)
#   }
#   
#   if (any(grepl("preferHab1", behavior))){
#     params[c("P1s","P1b","P1j")]<- c(.8,.8,.8)
#   }
#   
#   if (any(grepl("preferHab2", behavior))){
#     params[c("P1s","P1b","P1j")]<- c(.2,.2,.2)
#   }
#   
#   if (any(grepl("static", behavior))){
#     params[c("c1", "c2", "cF")]<- c(0,0,0)
#   }
#   
#   return (params)
# }

getParamsCombination.LH_Beh.ssa<- function(strategies=c("slow", "fast", "freqRepro"), 
                                           habDiffScenario=c("identicalHab", "mortalHab2", "nestPredHab2"), 
                                           behavior=c("neutral", "skip", "learnBreed", "learnExploreBreed", "preferHab1", "preferHab2")){
  comb<- expand.grid(strategy=strategies, habDiffScenario=habDiffScenario, behavior=behavior, stringsAsFactors=FALSE)
  params<- data.frame()
  for (i in 1:nrow(comb)){
    params<- rbind(params, getParams.LH_Beh.ssa(strategy=comb$strategy[i], habDiffScenario=comb$habDiffScenario[i], behavior=comb$behavior[i]))
  }
  
  return (params)
}


## Graphics ----
plotSim<- function(sim, groups=c("habitat", "age", "habitat*age", "all")[1], ...){
  if (is.integer(groups)) groups<- c("habitat", "age", "habitat*age", "all")[groups]
  tmp<- sim
  switch(groups,
         all={
           color<- grDevices::rainbow(length(grep("N1", colnames(sim))), start=0, end=1/7)
           color<- c(color, rev(grDevices::rainbow(length(grep("N2", colnames(sim))), start=2/6, end=4/7)))
           ylab<- "N by state"
         },
         habitat={
            nRow<- nrow(sim)
            sim<- matrix(nrow=nRow, ncol=3, dimnames=list(timeSeries=rep("", nRow), c("", "Hab1", "Hab2")))
            sim[,1]<-  tmp[,1]
            sim[,2]<-  apply(tmp[,2:5], 1, sum)
            sim[,3]<-  apply(tmp[,6:9], 1, sum)
            
            color<- grDevices::rainbow(length(grep("1", colnames(sim))), start=0, end=1/7)
            color<- c(color, rev(grDevices::rainbow(length(grep("2", colnames(sim))), start=2/6, end=4/7)))
            ylab<- "N by Habitat"
         },
         age={
           nRow<- nrow(sim)
            sim<- matrix(nrow=nRow, ncol=3, dimnames=list(timeSeries=rep("", nRow), c("", "Ad", "Juv")))
            sim[,1]<- tmp[,1]
            sim[,2]<- apply(tmp[,c(2:4, 6:8)], 1, sum)
            sim[,3]<- apply(tmp[,c(5,9)], 1, sum)
            
            color<- grDevices::rainbow(2, start=0, end=4/6)
            ylab<- "N by Age"
         },
         `habitat*age`={
            nRow<- nrow(sim)
            sim<- matrix(nrow=nRow, ncol=5, dimnames=list(timeSeries=rep("", nRow), c("", "AdHab1", "JuvHab1", "AdHab2", "JuvHab2")))
            sim[,1]<-  tmp[,1]
            sim[,2]<-  apply(tmp[,2:4], 1, sum)
            sim[,3]<-  tmp[,5]
            sim[,4]<-  apply(tmp[,6:8], 1, sum)
            sim[,5]<-  tmp[,9]
            
            color<- grDevices::rainbow(length(grep("1", colnames(sim))), start=0, end=1/8)
            color<- c(color, rev(grDevices::rainbow(length(grep("2", colnames(sim))), start=2/6, end=4/7)))
            ylab<- "N by Habitat and Age"
         }
  )
  
  graphics::matplot(sim[, 1], sim[,2:ncol(sim)], pch=19, cex=0.1, col=color, bty="n", xlab="Time", ylab=ylab, ...)
  graphics::legend("topright", legend = colnames(sim)[-1], bty = "y", pch = 19, col=color)
}

