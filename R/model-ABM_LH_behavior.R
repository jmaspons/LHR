
## Discrete time transitions with a probability model ----
# params<- slow=c(clutch1=1, clutch2=1,   b=1, PbF1=.4, PbF2=.4,  d1=.1,db1=.25,dj1=.25,  d2=.1,db2=.25,dj2=.25, g1=1, g2=1, K=500)
# params<- c(params, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5) # add neutral behavior
#' @importFrom stats rbinom 
transitionABM.LH_Beh<- function(N=matrix(rep(5, 8), nrow=4, ncol=8, dimnames=list(replicates=NULL, state=c("N1s", "N1b", "N1bF", "N2s", "N2b", "N2bF"))),
                         params=list(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5),
                         t){
  N0<- N
  nRep<- nrow(N)
  
  ## Growth TODO: non adult survival for AFR > 1

  N1j<- N2j<- numeric(nRep) # Juveniles
  
  ## Breeding
  for (i in 1:params$broods){
    ## Recruitment including juvenile survival (AFR = 1) and reproductive state change
    adultN1<- rowSums(N[,c("N1b", "N1bF", "N1s")])
    adultN2<- rowSums(N[,c("N2b", "N2bF", "N2s")])
    
    ## Breeding attemps
    breedingN1<-  with(params, rbinom(nRep, size=adultN1, prob=Pb1))
    breedingN2<-  with(params, rbinom(nRep, size=adultN2, prob=Pb2))
    
    # Skip reproduction
    N[,"N1s"]<- adultN1 - breedingN1
    N[,"N2s"]<- adultN2 - breedingN2
    
    ## Breeding success
    N[,"N1b"]<-  with(params, rbinom(nRep, size=breedingN1, prob=(1 - PbF1)))
    N[,"N2b"]<-  with(params, rbinom(nRep, size=breedingN2, prob=(1 - PbF2)))
    ## Breeding Fail
    N[,"N1bF"]<- breedingN1 - N[,"N1b"]
    N[,"N2bF"]<- breedingN2 - N[,"N2b"]
    
    ## Juvenile recruitment and mortality
    N1j<- N1j + with(params, rbinom(nRep, size=N[,"N1b"] * b1, prob=j1))
    N2j<- N2j + with(params, rbinom(nRep, size=N[,"N2b"] * b2, prob=j2))
    
    ## interbreed interval mortality
    
    ## Movements  (juveniles don't change habitat)
    # habOLD_NEW
    hab2_1Nb<-  with(params, rbinom(nRep, size=N[,"N2b"], prob=c2 * P1b))
    hab1_2Nb<-  with(params, rbinom(nRep, size=N[,"N1b"], prob=c1 * (1 - P1b)))
    hab2_1Ns<-  with(params, rbinom(nRep, size=N[,"N2s"], prob=c2 * P1b))
    hab1_2Ns<-  with(params, rbinom(nRep, size=N[,"N1s"], prob=c1 * (1 - P1b)))
    hab2_1NbF<-  with(params, rbinom(nRep, size=N[,"N2bF"], prob=cF * P1b))
    hab1_2NbF<-  with(params, rbinom(nRep, size=N[,"N1bF"], prob=cF * (1 - P1b)))
#     hab2_1Nj<-  with(params, rbinom(nRep, size=N2j, prob=c2 * P1j))
#     hab1_2Nj<-  with(params, rbinom(nRep, size=N1j, prob=c1 * (1 - P1j)))
    
    ## Apply movements
    N[,"N1b"]<-  N[,"N1b"]  + hab2_1Nb  - hab1_2Nb
    N[,"N2b"]<-  N[,"N2b"]  + hab1_2Nb  - hab2_1Nb
    N[,"N1s"]<-  N[,"N1s"]  + hab2_1Ns  - hab1_2Ns
    N[,"N2s"]<-  N[,"N2s"]  + hab1_2Ns  - hab2_1Ns
    N[,"N1bF"]<- N[,"N1bF"] + hab2_1NbF - hab1_2NbF
    N[,"N2bF"]<- N[,"N2bF"] + hab1_2NbF - hab2_1NbF
#     N1j<-  N1j + hab2_1Nj  - hab1_2Nj
#     N2j<-  N2j + hab1_2Nj  - hab2_1Nj
  }
  
  ## Survival
  # N[,"N1j"]<- # juvenile survival already calculated during recruitment.
  # N[,"N2j"]<- 
  
  N[,"N1b"]<-  with(params, rbinom(nRep, size=N[,"N1b"], prob=ab1))
  N[,"N2b"]<-  with(params, rbinom(nRep, size=N[,"N2b"], prob=ab2))

  N[,"N1bF"]<- with(params, rbinom(nRep, size=N[,"N1bF"], prob=ab1))
  N[,"N2bF"]<- with(params, rbinom(nRep, size=N[,"N2bF"], prob=ab2))

  N[,"N1s"]<-  with(params, rbinom(nRep, size=N[,"N1s"], prob=a1))
  N[,"N2s"]<-  with(params, rbinom(nRep, size=N[,"N2s"], prob=a2))
  
  ## Juveniles grow to Nxb classes
  N[,"N1b"]<- N[,"N1b"] + N1j
  N[,"N2b"]<- N[,"N2b"] + N2j
  
  return(N)
}

transitionABM.LH_Beh<- compiler::cmpfun(transitionABM.LH_Beh) # byte-compile the function


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

