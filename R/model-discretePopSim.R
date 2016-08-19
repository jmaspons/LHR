#' Discrete Time Simulation
#' 
#' \code{discretePopSim} class represent the result of discrete time simulations with
#' replicates on rows and time on columns. The class inherits from \code{matrix} for 
#' subsetting by timesteps or replicates.
#' 
#' @name discretePopSim
#' @importFrom stats rbinom rbeta
#' @exportClass discretePopSim
NULL

# When the population gets extinct it fills results with NAs.
extinctNA<- function(pop){
  extT<- apply(pop, 1, function(x) match(0, x))
  extT[extT == ncol(pop)]<- NA 
  extPop<- which(!is.na(extT))
  for (i in extPop){
    pop[i,(extT[i]+1):ncol(pop)]<- NA
  }
  return(pop)
}

## S4 model dispatcher ----

#' Discrete time and population models
#'
#' TODO: add varA
#'
#' @rdname discretePopSim
#' @param broods 
#' @param b 
#' @param j 
#' @param a 
#' @param breedFail 
#' @param varJ 
#' @param varBreedFail 
#' @param sexRatio 
#' @param matingSystem 
#' @param N0 
#' @param replicates 
#' @param tf 
#' @param maxN 
#'
#' @return a \code{discretePopSim} object.
#' @examples discretePopSim(b=2, j=.5, a=.5)
#'
#' @export
setGeneric("discretePopSim", function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0, 
                                      sexRatio=.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"), 
                                      N0=10, replicates=100, tf=10, maxN=100000) standardGeneric("discretePopSim"))
#### Common parameters:
# j=j, a=a, N0=N0, replicates=replicates, tf=tf, maxN=maxN
# 
#### Parameters by path:
# if(any(breedFail != 0)) # brood mortality: 
#   b=b, broods=broods, breedFail=breedFail
# else
#   b=b * broods
# 
# if (!is.na(matingSystem)) # two sexes: 
#   sexRatio=sexRatio, matingSystem=matingSystem
# 
# if (varJ != 0 | varBreedFail !=0) # environmental variation: 
#   varJ=varJ, varBreedFail=varBreedFail
# 
# if (length(j) > 1 | length(a) > 1) # environmental seasonality:
  
setMethod("discretePopSim", 
          signature(broods="ANY", b="ANY", j="numeric", a="numeric", breedFail="ANY", varJ="ANY", varBreedFail="ANY",
                    sexRatio="ANY", matingSystem="ANY", N0="ANY", replicates="ANY", tf="ANY", maxN="ANY"),
          function(broods=1, b, j, a, breedFail=0, varJ=0, varBreedFail=0,
                   sexRatio=.5, matingSystem=c(NA, "monogamy", "polygyny", "polyandry"),
                   N0=10, replicates=100, tf=10, maxN=100000){
            
            if (!all(is.na(matingSystem))){
              matingSystem<- match.arg(matingSystem)
            }          
            
            broodMortality<- sex<- varEnv<- seasonalEnv<- FALSE

            if(any(breedFail != 0)) broodMortality<- TRUE         # differential brood mortality
            if (varJ != 0 | varBreedFail !=0) varEnv<- TRUE       # environmental variation
            if (length(j) > 1 | length(a) > 1) seasonalEnv<- TRUE # environmental seasonality
            if (!is.na(matingSystem)) sex<- TRUE                  # two sexes

            if (!broodMortality & !varEnv & !seasonalEnv & !sex){
              out<- mFit.t(fecundity=b, j=j, a=a, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & !seasonalEnv &  sex){
              out<- mFitSex.t(fecundity=b, j=j, a=a, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & !varEnv & seasonalEnv  & !sex){
              ## TODO: out<- mFit.tseason
            }
            if (!broodMortality & !varEnv & seasonalEnv  &  sex){
              ## TODO: out<- mFitSex.tseason
            }
            if (!broodMortality & varEnv  & !seasonalEnv & !sex){
              out<- mFit.tvar(fecundity=b, j=j, a=a, varJ=varJ, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & !seasonalEnv &  sex){
              out<- mFitSex.tvar(fecundity=b, j=j, a=a, varJ=varJ, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (!broodMortality & varEnv  & seasonalEnv  & !sex){
              ## TODO: out<- mFit.tvarseason
            }
            if (!broodMortality & varEnv  & seasonalEnv  &  sex){
              ## TODO: out<- mFitSex.tvarseason
            }
            
            if (broodMortality  & !varEnv & !seasonalEnv & !sex){
              out<- mSurvBV.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & !seasonalEnv &  sex){
              out<- mSurvBVSex.t(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  & !sex){
              out<- mSurvBV.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & !varEnv & seasonalEnv  &  sex){
              out<- mSurvBVSex.tseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv & !sex){
              out<- mSurvBV.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & !seasonalEnv &  sex){
              out<- mSurvBVSex.tvar(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  & !sex){
              out<- mSurvBV.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            if (broodMortality  & varEnv  & seasonalEnv  &  sex){
              out<- mSurvBVSex.tvarseason(broods=broods, b=b, j=j, a=a, breedFail=breedFail, varJ=varJ, varBreedFail=varBreedFail, sexRatio=sexRatio, matingSystem=matingSystem, N0=N0, replicates=replicates, tf=tf, maxN=maxN)
            }
            
            return (out)
          }
)
