
## Discrete time transitions with a probability model ----
# Stage based subadult class. 
# TODO: Add age based subadults class? Important for time lags
# params<- slow=c(clutch1=1, clutch2=1,   b=1, PbF1=.4, PbF2=.4,  d1=.1,db1=.25,dj1=.25,  d2=.1,db2=.25,dj2=.25, g1=1, g2=1)
# params<- c(params, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5) # add neutral behavior
#' @importFrom stats rbinom 
transitionABM.LH_Beh<- function(N=matrix(rep(5, 6 * 4), nrow=4, ncol=6, dimnames=list(replicates=NULL, state=c("N1s", "N1b", "N1bF", "N2s", "N2b", "N2bF"))),
                         params=list(b1=1, b2=1,  broods=1, PbF1=.4, PbF2=.4,  a1=.25,ab1=.1,sa1=.25,j1=.1,  a2=.25,ab2=.1,sa2=.25,j2=.1, AFR=1, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1sa=.5, P1j=.5)){
  if (is.null(dim(N)))
    N<- as.matrix(t(N))
  
  nRep<- nrow(N)
  
  N1j<- N2j<- numeric(nRep) # Juveniles
  
  saAges<- grep("sa", colnames(N), value=TRUE)
  if (length(saAges) / 2 < params$AFR - 1){
    warning("Subadult classes missing in the population matrix N. Updating N") 
    if (sum(N[, saAges]) > 0) warning("Removing individuals from subadult classes ", saAges)
    
    saColsOri<- grep("sa", colnames(N))
    if (length(saColsOri) > 0)
      N<- N[, -saColsOri] # remove original columns
    
    saAges<- paste0("sa", 1:(params$AFR - 1))
    saAges<- c(paste0("N1", saAges), paste0("N2", saAges))
    Nsa<- matrix(0, nrow=nRep, ncol=length(saAges), dimnames=list(replicates=NULL, state=saAges))
    N<- cbind(N, Nsa)
    rm(Nsa)
  }
  
  saAges1<- grep("N1", saAges, value=TRUE)
  saAges2<- grep("N2", saAges, value=TRUE)
  
  
  ## BREEDING
  adultN1<- rowSums(N[,c("N1b", "N1bF", "N1s"), drop=FALSE])
  adultN2<- rowSums(N[,c("N2b", "N2bF", "N2s"), drop=FALSE])
  
  ## Breeders
  breedingN1<-  with(params, rbinom(nRep, size=adultN1, prob=Pb1))
  breedingN2<-  with(params, rbinom(nRep, size=adultN2, prob=Pb2))
  
  # Skip reproduction
  N[,"N1s"]<- adultN1 - breedingN1
  N[,"N2s"]<- adultN2 - breedingN2
  
  for (i in 1:params$broods){
    ## Breeding success
    N[,"N1b"]<-  with(params, rbinom(nRep, size=breedingN1, prob=(1 - PbF1)))
    N[,"N2b"]<-  with(params, rbinom(nRep, size=breedingN2, prob=(1 - PbF2)))
    ## Breeding Fail
    N[,"N1bF"]<- breedingN1 - N[,"N1b"]
    N[,"N2bF"]<- breedingN2 - N[,"N2b"]
    
    ## Juvenile recruitment and mortality
    N1j<- N1j + with(params, rbinom(nRep, size=N[,"N1b"] * b1, prob=j1))
    N2j<- N2j + with(params, rbinom(nRep, size=N[,"N2b"] * b2, prob=j2))
    
    ## interbreed interval mortality? otherwise mortality only affected by the last brood event
    ## Skip reproduction moved outside the broods loop
    
    ## MOVEMENTS
    # habOLD_NEW
    hab2_1Nb<-  with(params, rbinom(nRep, size=N[,"N2b"], prob=c2 * P1b))
    hab1_2Nb<-  with(params, rbinom(nRep, size=N[,"N1b"], prob=c1 * (1 - P1b)))
    hab2_1NbF<-  with(params, rbinom(nRep, size=N[,"N2bF"], prob=cF * P1b))
    hab1_2NbF<-  with(params, rbinom(nRep, size=N[,"N1bF"], prob=cF * (1 - P1b)))

    ## Apply movements
    N[,"N1b"]<-  N[,"N1b"]  + hab2_1Nb  - hab1_2Nb
    N[,"N2b"]<-  N[,"N2b"]  + hab1_2Nb  - hab2_1Nb
    N[,"N1bF"]<- N[,"N1bF"] + hab2_1NbF - hab1_2NbF
    N[,"N2bF"]<- N[,"N2bF"] + hab1_2NbF - hab2_1NbF
    
    breedingN1<-  N[,"N1b"] +  N[,"N1bF"]
    breedingN2<-  N[,"N2b"] +  N[,"N2bF"]
  }
  
  ## MOVEMENTS (not based on breeding experience, only once per year)
  # non reproductive adults (juveniles don't change habitat)
  hab2_1Ns<-  with(params, rbinom(nRep, size=N[,"N2s"], prob=c2 * P1s))
  hab1_2Ns<-  with(params, rbinom(nRep, size=N[,"N1s"], prob=c1 * (1 - P1s)))
  # hab2_1Nj<-  with(params, rbinom(nRep, size=N2j, prob=c2 * P1j))
  # hab1_2Nj<-  with(params, rbinom(nRep, size=N1j, prob=c1 * (1 - P1j)))
  
  N[,"N1s"]<-  N[,"N1s"]  + hab2_1Ns  - hab1_2Ns
  N[,"N2s"]<-  N[,"N2s"]  + hab1_2Ns  - hab2_1Ns
  #     N1j<-  N1j + hab2_1Nj  - hab1_2Nj
  #     N2j<-  N2j + hab1_2Nj  - hab2_1Nj
  
  # subadults
  if (params$AFR > 1){
    hab2_1Nsa<-  apply(N[, saAges2, drop=FALSE], 2, function(x) with(params, rbinom(nRep, size=x, prob=c2 * P1sa)))
    hab1_2Nsa<-  apply(N[, saAges1, drop=FALSE], 2, function(x) with(params, rbinom(nRep, size=x, prob=c1 * (1 - P1sa))))
    N[, saAges1]<-  N[, saAges1, drop=FALSE] + hab2_1Nsa - hab1_2Nsa
    N[, saAges2]<-  N[, saAges2, drop=FALSE] + hab1_2Nsa - hab2_1Nsa
  }
  
  
  ## SURVIVAL
  # juvenile survival already calculated during recruitment (saved in N*j)
  N[,"N1b"]<-  with(params, rbinom(nRep, size=N[,"N1b"], prob=ab1))
  N[,"N2b"]<-  with(params, rbinom(nRep, size=N[,"N2b"], prob=ab2))

  N[,"N1bF"]<- with(params, rbinom(nRep, size=N[,"N1bF"], prob=ab1))
  N[,"N2bF"]<- with(params, rbinom(nRep, size=N[,"N2bF"], prob=ab2))

  N[,"N1s"]<-  with(params, rbinom(nRep, size=N[,"N1s"], prob=a1))
  N[,"N2s"]<-  with(params, rbinom(nRep, size=N[,"N2s"], prob=a2))
  
  ## Subadults
  if (params$AFR > 1){
    N1sa<-  apply(N[, saAges1, drop=FALSE], 2, function(x) with(params, rbinom(nRep, size=x, prob=sa1)))
    N2sa<-  apply(N[, saAges2, drop=FALSE], 2, function(x) with(params, rbinom(nRep, size=x, prob=sa2)))
    
    if (is.null(dim(N1sa))){
      N1sa<- as.matrix(t(N1sa))
      N2sa<- as.matrix(t(N2sa))
    }
  }
  
  
  ## GROWTH
  if (params$AFR > 1){
    # WARNING: assumes saAges* is sorted by age
    N[, saAges1[-1]]<- N1sa[, -ncol(N1sa)] # subadult age transitions
    N[, saAges2[-1]]<- N2sa[, -ncol(N2sa)] # subadult age transitions
    N[, "N1b"]<- N[, "N1b"] + N1sa[, ncol(N1sa)] # subadult -> adult
    N[, "N2b"]<- N[, "N2b"] + N2sa[, ncol(N2sa)] # subadult -> adult
    N[, saAges1[1]]<- N1j # juvenil -> subadult
    N[, saAges2[1]]<- N2j # juvenil -> subadult
  } else { # Juveniles grow to N*b
    N[,"N1b"]<- N[,"N1b"] + N1j
    N[,"N2b"]<- N[,"N2b"] + N2j
  }
  
  return(N)
}

transitionABM.LH_Beh<- compiler::cmpfun(transitionABM.LH_Beh) # byte-compile the function


transitionABM.LH_Beh_DET<- function(N=matrix(rep(5, 6), nrow=1, ncol=6, dimnames=list(NULL, state=c("N1s", "N1b", "N1bF", "N2s", "N2b", "N2bF"))),
                                params=list(b1=1, b2=1,  broods=1, PbF1=.4, PbF2=.4,  a1=.25,ab1=.1,sa1=.25,j1=.1,  a2=.25,ab2=.1,sa2=.25,j2=.1, AFR=1, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1sa=.5, P1j=.5)){
  if (is.null(dim(N)))
    N<- as.matrix(t(N))
  
  N1j<- N2j<- numeric(1) # Juveniles
  
  saAges<- grep("sa", colnames(N), value=TRUE)
  if (length(saAges) / 2 < params$AFR - 1){
    warning("Subadult classes missing in the population matrix N. Updating N")
    if (sum(N[, saAges]) > 0) warning("Removing individuals from subadult classes ", saAges)
    
    saColsOri<- grep("sa", colnames(N))
    if (length(saColsOri) > 0)
      N<- N[, -saColsOri] # remove original columns
    
    saAges<- paste0("sa", 1:(params$AFR - 1))
    saAges<- c(paste0("N1", saAges), paste0("N2", saAges))
    Nsa<- matrix(0, nrow=1, ncol=length(saAges), dimnames=list(replicates=1, state=saAges))
    N<- cbind(N, Nsa)
    rm(Nsa)
  }
  
  saAges1<- grep("N1", saAges, value=TRUE)
  saAges2<- grep("N2", saAges, value=TRUE)
  
  ## BREEDING
  adultN1<- sum(N[,c("N1b", "N1bF", "N1s")])
  adultN2<- sum(N[,c("N2b", "N2bF", "N2s")])
  
  ## Breeders
  breedingN1<-  with(params, adultN1 * Pb1)
  breedingN2<-  with(params, adultN2 * Pb2)
  
  # Skip reproduction
  N[,"N1s"]<- adultN1 - breedingN1
  N[,"N2s"]<- adultN2 - breedingN2
  
  for (i in 1:params$broods){
    ## Breeding success
    N[,"N1b"]<-  with(params, breedingN1 * (1 - PbF1))
    N[,"N2b"]<-  with(params, breedingN2 * (1 - PbF2))
    ## Breeding Fail
    N[,"N1bF"]<- breedingN1 - N[,"N1b"]
    N[,"N2bF"]<- breedingN2 - N[,"N2b"]
    
    ## Juvenile recruitment and mortality
    N1j<- N1j + with(params, N[,"N1b"] * b1 * j1)
    N2j<- N2j + with(params, N[,"N2b"] * b2 * j2)
    
    ## interbreed interval mortality? otherwise mortality only affected by the last brood event
    ## Skip reproduction moved outside the broods loop
    
    ## MOVEMENTS
    # habOLD_NEW
    hab2_1Nb<-  with(params, N[,"N2b"] * c2 * P1b)
    hab1_2Nb<-  with(params, N[,"N1b"] * c1 * (1 - P1b))
    hab2_1NbF<-  with(params, N[,"N2bF"] * cF * P1b)
    hab1_2NbF<-  with(params, N[,"N1bF"] * cF * (1 - P1b))
    
    ## Apply movements
    N[,"N1b"]<-  N[,"N1b"]  + hab2_1Nb  - hab1_2Nb
    N[,"N2b"]<-  N[,"N2b"]  + hab1_2Nb  - hab2_1Nb
    N[,"N1bF"]<- N[,"N1bF"] + hab2_1NbF - hab1_2NbF
    N[,"N2bF"]<- N[,"N2bF"] + hab1_2NbF - hab2_1NbF
    
    breedingN1<-  N[,"N1b"] +  N[,"N1bF"]
    breedingN2<-  N[,"N2b"] +  N[,"N2bF"]
  }
  
  ## MOVEMENTS (not based on breeding experience, only once per year)
  # non reproductive adults (juveniles don't change habitat)
  hab2_1Ns<-  with(params, N[,"N2s"] * c2 * P1s)
  hab1_2Ns<-  with(params, N[,"N1s"] * c1 * (1 - P1s))
  # hab2_1Nj<-  with(params, N2j * c2 * P1j)
  # hab1_2Nj<-  with(params, N1j * c1 * (1 - P1j))
  
  N[,"N1s"]<-  N[,"N1s"]  + hab2_1Ns  - hab1_2Ns
  N[,"N2s"]<-  N[,"N2s"]  + hab1_2Ns  - hab2_1Ns
  #     N1j<-  N1j + hab2_1Nj  - hab1_2Nj
  #     N2j<-  N2j + hab1_2Nj  - hab2_1Nj
  
  # subadults
  if (params$AFR > 1){
    hab2_1Nsa<-  apply(N[, saAges2, drop=FALSE], 2, function(x) with(params, x * c2 * P1sa))
    hab1_2Nsa<-  apply(N[, saAges1, drop=FALSE], 2, function(x) with(params, x * c1 * (1 - P1sa)))
    N[, saAges1]<-  N[, saAges1, drop=FALSE] + hab2_1Nsa - hab1_2Nsa
    N[, saAges2]<-  N[, saAges2, drop=FALSE] + hab1_2Nsa - hab2_1Nsa
  }
  
  
  ## SURVIVAL
  # juvenile survival already calculated during recruitment (saved in N*j)
  N[,"N1b"]<-  with(params, N[,"N1b"] * ab1)
  N[,"N2b"]<-  with(params, N[,"N2b"] * ab2)
  
  N[,"N1bF"]<- with(params, N[,"N1bF"] * ab1)
  N[,"N2bF"]<- with(params, N[,"N2bF"] * ab2)
  
  N[,"N1s"]<-  with(params, N[,"N1s"] * a1)
  N[,"N2s"]<-  with(params, N[,"N2s"] * a2)
  
  ## Subadults
  if (params$AFR > 1){
    N1sa<-  apply(N[, saAges1, drop=FALSE], 2, function(x) with(params, x * sa1))
    N2sa<-  apply(N[, saAges2, drop=FALSE], 2, function(x) with(params, x * sa2))
    
    if (is.null(dim(N1sa))){
      N1sa<- as.matrix(t(N1sa))
      N2sa<- as.matrix(t(N2sa))
    }
  }
  
  
  ## GROWTH
  if (params$AFR > 1){
    # WARNING: assumes saAges* is sorted by age
    N[, saAges1[-1]]<- N1sa[, -ncol(N1sa)] # subadult age transitions
    N[, saAges2[-1]]<- N2sa[, -ncol(N2sa)] # subadult age transitions
    N[, "N1b"]<- N[, "N1b"] + N1sa[, ncol(N1sa)] # subadult -> adult
    N[, "N2b"]<- N[, "N2b"] + N2sa[, ncol(N2sa)] # subadult -> adult
    N[, saAges1[1]]<- N1j # juvenil -> subadult
    N[, saAges2[1]]<- N2j # juvenil -> subadult
  } else { # Juveniles grow to N*b
    N[,"N1b"]<- N[,"N1b"] + N1j
    N[,"N2b"]<- N[,"N2b"] + N2j
  }
  
  return(N)
}

transitionABM.LH_Beh_DET<- compiler::cmpfun(transitionABM.LH_Beh_DET) # byte-compile the function


## Graphics ----

#' Plot LH_behavior simulations
#'
#' @param x a \code{\link{discreteABMSim}} or a \code{\link{Model.ABM}}
#' @param ... 
#'
#' @export
plotLH_behavior<- function(x, ...){ UseMethod("plotLH_behavior") }

#' plotLH_behavior
#' 
#' @describeIn discreteABMSim
#' @param x 
#' @param groups 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plotLH_behavior.discreteABMSim<- function(x, groups=c("all", "habitat", "age", "habitat*age"), ...){
  if (is.integer(groups)) groups<- c("habitat", "age", "habitat*age", "all")[groups]
  else groups<- match.arg(groups)
  
  tmp<- x
  
  switch(groups,
         all={
           color<- grDevices::rainbow(length(grep("N1", colnames(x))), start=0, end=1/7)
           color<- c(color, rev(grDevices::rainbow(length(grep("N2", colnames(x))), start=2/6, end=4/7)))
           ylab<- "N by state"
         },
         habitat={
           tmp<- array(dim=c(nrow(x), 2, dim(x)[3]), dimnames=list(replicate=NULL, state=c("Hab1", "Hab2"), t=0:(dim(x)[3] - 1)))
           tmp[, 1,]<-  apply(x[, grep("^N1", colnames(x)), ], 3, rowSums)
           tmp[, 2,]<-  apply(x[, grep("^N2", colnames(x)), ], 3, rowSums)
           
           color<- grDevices::rainbow(length(grep("1", colnames(tmp))), start=0, end=1/7)
           color<- c(color, rev(grDevices::rainbow(length(grep("2", colnames(tmp))), start=2/6, end=4/7)))
           ylab<- "N by Habitat"
         },
         age={
           tmp<- array(dim=c(nrow(x), 2, dim(x)[3]), dimnames=list(replicate=NULL, state=c("Ad", "subAd"), t=0:(dim(x)[3] - 1)))
           tmp[, 1,]<-  apply(x[, grep("^N[12]{1}[sbF]{1,2}$", colnames(x)), ], 3, rowSums)
           tmp[, 2,]<-  apply(x[, grep("^N[12]{1}sa", colnames(x)), ], 3, rowSums)
           
           color<- grDevices::rainbow(2, start=0, end=4/6)
           ylab<- "N by Age"
         },
         `habitat*age`={
            tmp<- array(dim=c(nrow(x), 4, dim(x)[3]), dimnames=list(replicate=NULL, state=c("AdHab1", "subAdHab1", "AdHab2", "subAdHab2"), t=0:(dim(x)[3] - 1)))
            tmp[, 1,]<-  apply(x[, grep("^N1[sbF]{1,2}$", colnames(x)), ], 3, rowSums)
            tmp[, 2,]<-  apply(x[, grep("^N1sa", colnames(x)), ], 3, rowSums)
            tmp[, 3,]<-  apply(x[, grep("^N2[sbF]{1,2}$", colnames(x)), ], 3, rowSums)
            tmp[, 4,]<-  apply(x[, grep("^N2sa", colnames(x)), ], 3, rowSums)
           
            color<- grDevices::rainbow(length(grep("1", colnames(tmp))), start=0, end=1/8)
            color<- c(color, rev(grDevices::rainbow(length(grep("2", colnames(tmp))), start=2/6, end=4/7)))
            ylab<- "N by Habitat and Age"
         }
  )
  
  # graphics::matplot(1:dim(tmp)[3], t(tmp[1,,]), pch=19, cex=0.1, col=color, bty="n", xlab="Time", ylab=ylab, lty=1, type="b", ...)
  # graphics::legend("topright", legend = colnames(tmp)[-1], bty = "y", pch = 19, col=color)
  
  # tmpMean<- colMeans(tmp)
  # tmpMean<- data.frame(t=as.numeric(colnames(tmpMean)), t(tmpMean))
  # tmpMean<- reshape2::melt(tmpMean, id.vars="t", value.name="meanN")
  # names(tmpMean)[2]<- "state"
  # tmpMean$id<- paste0(tmpMean$t, "_", tmpMean$state)
  
  tmpQuantile<- apply(tmp, 2:3, quantile)
  tmpQuantile<- apply(tmpQuantile, 3, function(x){
    out<- reshape2::melt(x, value.name="N")
    colnames(out)[1]<- "quantile"
   
    out
  })
  tmpQuantile<- lapply(seq_along(tmpQuantile), function(x){
    # out<- as.data.frame(t(tmpQuantile[[x]]), stringsAsFactors=FALSE)
    # names(out)<- paste0(out["state",], "-", out["quantile",])
    # out<- out["N",]
    out<- tmpQuantile[[x]]
    out$t<- x
    out
  })
  tmpQuantile<- do.call(rbind, tmpQuantile)
  
  tmpQuantile$color<- tmpQuantile$state
  
  # Avoid -Inf when scale_y_log10
  if (any(tmpQuantile$N == 0)) tmpQuantile$N<- tmpQuantile$N + .5
  
  tmpMean<- tmpQuantile[tmpQuantile$quantile == "50%",]
  tmpQuantile<- tmpQuantile[tmpQuantile$quantile %in% c("25%", "75%"),]
  tmpQuantile<- reshape2::dcast(tmpQuantile, t + state + color ~ quantile, value.var="N")
  # tmp<- reshape2::melt(tmpQuantile, id.vars="t", value.name="N")
  
  ggplot2::ggplot(tmpMean, ggplot2::aes(x=t, y=N, color=color, group=state)) + ggplot2::scale_y_log10() +
    ggplot2::geom_ribbon(data=tmpQuantile, mapping=ggplot2::aes(x=t, ymin=`25%`, ymax=`75%`, fill=color, color=color), alpha=.4, linetype=0, inherit.aes=FALSE) + ggplot2::geom_line() + 
    ggplot2::labs(x="Time", y=expression(N[t] * " (mean and 50%)"))
}


#' Plot Model
#'
#' @describeIn Model.ABM
#' @param x 
#' @param resultType 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plotLH_behavior.Model<- function(x, resultType=c("Pest_N0", "G", "N0_Pest", "Ntf"), facet_grid=breedFail ~ idHabDiff + idBehavior, ...){
  if (missing(resultType)) noType<- TRUE else noType<- FALSE
  
  resultType<- match.arg(resultType)
  
  stats<- nrow(x@sim) > 0
  N0_Pest<- nrow(x@sim@N0_Pest) > 0
  Ntf<- nrow(x@sim@Ntf) > 0
  
  # if no type is specified select a existing one with precedence stats > N0_Pest > Ntf
  if (noType){
    if (stats) resultType<- "Pest_N0"
    else if (N0_Pest) resultType<- "N0_Pest"
    else if (Ntf) resultType<- "Ntf"
    else {
      invisible(graphics::plot(x))
    }
  }
  
  out<- NA
  
  if (stats & resultType == "Pest_N0"){
    res<- result(x, type="stats")

    out<- plotPest_N0.LH_behavior(res, facet_grid=facet_grid, ...)
  }
  
  if (stats & resultType == "G"){
    res<- result(x, type="stats")

    out<- plotG.LH_behavior(res, facet_grid=facet_grid, ...)
  }
  
  if (N0_Pest & resultType == "N0_Pest"){  
    res<- result(x, type="N0_Pest")

    out<- plotN0_Pest.LH_behavior(res, facet_grid=facet_grid, ...)
  }
  
  if (Ntf & resultType == "Ntf"){
    res<- result(x, type="Ntf")

    out<- plotNtf.LH_behavior(res, facet_grid=facet_grid, ...)
  }
  
  # out<- out + ggplot2::scale_color_gradient2(low="red", mid="yellow", high="green") +
  #   ggplot2::scale_fill_gradient2(low="red", mid="yellow", high="green")
  
  # out<- out + ggplot2::scale_color_gradientn(colors=grDevices::terrain.colors(nrow(lh))) +
  #   ggplot2::scale_fill_gradientn(colors=grDevices::terrain.colors(nrow(lh)))
  
  return(out)
}

#' PlotLH_behavior
#' 
#' @describeIn plotLH_behavior
#' @param x a \code{data.frame} from \code{\link{result}}.
#' @param resultType 
#' @param facet_grid 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plotLH_behavior.data.frame<- function(x, resultType=c("Pest_N0", "G", "N0_Pest", "Ntf"), facet_grid=breedFail ~ idHabDiff + idBehavior, ...){
  if (missing(resultType)) noType<- TRUE else noType<- FALSE
  
  resultType<- match.arg(resultType)
  
  stats<- all(c("idScenario", "N0", "GL", "colorLH", "lambda", "seasonAmplitude", "varJ", "varA", "breedFail") %in% names(x))
  N0_Pest<- all(c("idScenario", "N0interpoled", "colorLH", "lambda", "seasonAmplitude", "varJ", "varA", "breedFail") %in% names(x))
  Ntf<- all(c("idScenario", "N0", "25%", "50%", "75%", "colorLH", "lambda", "seasonAmplitude", "varJ", "varA", "breedFail") %in% names(x))
  
  # if no type is specified select a existing one with precedence stats > N0_Pest > Ntf
  if (noType){
    if (stats) resultType<- "Pest_N0"
    else if (N0_Pest) resultType<- "N0_Pest"
    else if (Ntf) resultType<- "Ntf"
    else {
      invisible(plot(x))
    }
  }
  
  out<- NA
  
  if (stats & resultType == "Pest_N0"){
    x$Pest<- 1 - x$extinct
    
    out<- plotPest_N0.LH_behavior(x, facet_grid=facet_grid, ...)
  }
  
  if (stats & resultType == "G"){
    out<- plotG.LH_behavior(x, facet_grid=facet_grid, ...)
  }
  
  if (N0_Pest & resultType == "N0_Pest"){  
    out<- plotN0_Pest.LH_behavior(x, facet_grid=facet_grid, ...)
  }
  
  if (Ntf & resultType == "Ntf"){
    out<- plotNtf.LH_behavior(x, facet_grid=facet_grid, ...)
  }
  
  # out<- out + ggplot2::scale_color_gradient2(low="red", mid="yellow", high="green") +
  #   ggplot2::scale_fill_gradient2(low="red", mid="yellow", high="green")
  
  # out<- out + ggplot2::scale_color_gradientn(colors=grDevices::terrain.colors(nrow(lh))) +
  #   ggplot2::scale_fill_gradientn(colors=grDevices::terrain.colors(nrow(lh)))
  
  return(out)
}

plotPest_N0.LH_behavior<- function(x, facet_grid=breedFail ~ idHabDiff + idBehavior, ...){
  x$Pest<- 1 - x$extinct
  
  ggplot2::ggplot(data=x, ggplot2::aes(x=N0, y=Pest, group=idScenario, color=colorLH)) + 
    ggplot2::geom_line(mapping=ggplot2::aes(size=lambda), alpha=0.5) + ggplot2::geom_point() +
    ggplot2::facet_grid(facet_grid, labeller=ggplot2::label_both) +
    ggplot2::scale_size(breaks=unique(x$lambda), range=c(0.5, 2))
}

plotG.LH_behavior<- function(x, facet_grid=breedFail ~ idHabDiff + idBehavior, ...){
  ggplot2::ggplot(data=x, ggplot2::aes(x=N0, y=GL, group=idScenario, color=colorLH)) + 
    ggplot2::geom_hline(yintercept=1) + ggplot2::geom_line(mapping=ggplot2::aes(size=lambda), alpha=0.5) + # ggplot2::geom_point() +
    ggplot2::geom_line(mapping=ggplot2::aes(y=meanL), linetype=2) +
    ggplot2::facet_grid(facet_grid, labeller=ggplot2::label_both) +
    ggplot2::labs(x=expression(N[0]), y=expression(lambda * " geometric mean and mean (dashed) for all transitions")) +
    ggplot2::scale_size(breaks=unique(x$lambda), range=c(0.5, 2))
}

plotN0_Pest.LH_behavior<- function(x, facet_grid=breedFail ~ idHabDiff + idBehavior, ...){
  ggplot2::ggplot(data=x, ggplot2::aes(x=lambda, y=N0interpoled, group=colorLH, color=colorLH)) + ggplot2::scale_y_log10() +
    ggplot2::geom_point() + ggplot2::geom_line(size=0.5) + ggplot2::facet_grid(facet_grid, labeller=ggplot2::label_both)  
}

plotNtf.LH_behavior<- function(x, facet_grid=breedFail ~ idHabDiff + idBehavior, ...){
  abline<- 1:max(x$N0)
  abline<- data.frame(x=abline, y=abline)
  
  # Avoid -Inf when scale_y_log10
  if (any(x[, grep("%", names(x), value=T)] == 0)){
    x[, grep("%", names(x), value=T)]<- x[, grep("%", names(x), value=T)] + 0.5
  }
  
  ggplot2::ggplot(data=x, ggplot2::aes(x=N0, color=colorLH, fill=colorLH, group=idScenario)) + ggplot2::scale_y_log10() + 
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=`25%`, ymax=`75%`), alpha=.4, linetype=0) + ggplot2::geom_line(mapping=ggplot2::aes(y=`50%`, size=lambda)) + 
    ggplot2::facet_grid(facet_grid, labeller=ggplot2::label_both) + ggplot2::labs(x=expression(N[0]), y=expression(N[tf] * " (mean and 50%)")) +
    ggplot2::scale_size(breaks=unique(x$lambda), range=c(0.5, 2)) + ggplot2::geom_line(mapping=ggplot2::aes(x, y), data=abline, inherit.aes=FALSE, show.legend=FALSE, lty=2)
}
