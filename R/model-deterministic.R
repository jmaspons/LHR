## Parameters
# AFR: Age at First Reproduction
# a, j, s: survival for adults (A), juveniles (J) and subadults (S)
# b: fecundity
# bj=F: net fecundity
# e.g.
# b<-3; a<-.9; s=.7; j=.5; bj=F= b * j

## Matricial models ----
## Pre-breeding census
# AFR=1
# b*j b*j
#   a   a
# 
# AFR=2
# 0 b*j
# s a


#' Lefkovich matrix
#' 
#' @name lefkovitch 
#' @param a 
#' @param s
#' @param j
#' @param b 
#' @param bj 
#' @param AFR 
NULL

#' @rdname lefkovitch
#' @export
LefkovitchPre<- function(a, s, bj, AFR=1){
  rows<- max(AFR, 2)
  mat<- matrix(0, nrow=rows, ncol=rows)
  mat[1, rows]<- bj
  mat[rows, rows]<- a
  stages<- rep("A", 2)
  if (AFR > 1){
    stages<- c(rep("S", rows-1), "A")
    for (i in 2:(rows)){
      mat[i,i-1]<- s
    }
  }else{
    mat[,1]<- mat[,2]
  }
  dimnames(mat)<- list(stages, stages)
  class(mat)<- "leslieMatrix"
  return (mat)
}

## Post-breeding census
# AFR=1
# 0 b
# j a
# 
# AFR=2
# 0 0 b
# j 0 0
# 0 s a

#' @rdname lefkovitch
#' @export
LefkovitchPost<- function(a, s, j, b, AFR=1){
  rows<- AFR + 1
  stages<- c("J", rep("S", rows-2), "A")
  mat<- matrix(0, nrow=rows, ncol=rows, dimnames=list(stages, stages))
  mat[1, rows]<- b
  mat[rows, rows]<- a
  mat[2, 1]<- j
  if (AFR >= 2){
    for (i in 3:(rows)){
      mat[i,i-1]<- s
    }
  }
  class(mat)<- "leslieMatrix"
  return (mat)
}

#' @export
lambda<- function(x, ...){
  UseMethod("lambda")  
}

#' Lambda
#'
#' @param mat 
#'
#' @return
#' @export
#'
#' @examples
lambda.leslieMatrix<- function(x, ...){ # Ted J. Case. "An Illustrated Guide to Theoretical Ecology" pag. 64-73
  max(Re(eigen(x, only.values=TRUE)$values))
}

# Use eigen.analysis function to analyse matrices. It is imported from popbio package
# popbio::eigen.analysis(mat)

## Inverse eigenvalue problem: find a, j, b or F from lambda and the rest of parameters ----
# Not so easy, some parameters fall outside the domain (e.g. survival <0 or >1): 

## CORRESPONDS TO A PRE-BREEDING CENSUS MATRIX (see test/Euler-Lotka.r)
# Euler-Lotka equation fails for some cases with no clear pattern (see test/Euler-Lotka.r)
# http://en.wikipedia.org/wiki/Euler-Lotka_equation
# http://darwin.eeb.uconn.edu/eeb310/lecture-notes/spotted-owl/node3.html
# see exec/Euler-Lotka.txt & inst/Euler-lotka.cws (cantor)
# alpha= AFR
#(l_{alpha}*F)*lambda^-alpha sumatori(from x=alpha to infinite){lambda^(-x+alpha)*s^(x-alpha)} = 1
#(l_{alpha}*F)*(1/(1-(a/lambda))) = lambda^alpha
# lambda^alpha*(1-(a/lambda)) = (l_{alpha}*F)
# l_{alpha}:  survival until reproductive class -> j * a^(alpha-1)
# lambda^alpha*(1-(a/lambda)) = j * a^(alpha-1) * F

## Maxima
# lambda^alpha * (1 - (a / lambda)) = j * s^(alpha-1) * F;
#
# solve(lambda^alpha * (1 - (a / lambda)) = j * s^(alpha-1) * F, F);
# solve(lambda^alpha * (1 - (a / lambda)) = j * s^(alpha-1) * j * b, a);
# solve(lambda^alpha * (1 - (a / lambda)) = j * s^(alpha-1) * j * b, j);
# solve(lambda^alpha * (1 - (a / lambda)) = j * s^(alpha-1) * j * b, b);
# solve(lambda^alpha * (1 - (a / lambda)) = j * s^(alpha-1) * j * b, alpha);

#' Euler-Lotka equations
#'
#' @name EulerLotka
#' @param lambda 
#' @param a adult survival
#' @param s subadult survival
#' @param j juvenile survival
#' @param b litter size
#' @param AFR age at first reproduction
NULL

#' @rdname EulerLotka
#' @export
findF_EulerLotka<- function(lambda, a, s=a, AFR){
  F<- -(a * lambda^AFR - lambda^(AFR+1)) / (s^(AFR-1) * lambda)
  F[F < 0]<- NA
  return (F)
}

#' @rdname EulerLotka
#' @export
findJ_EulerLotka<- function(lambda, b, a, s=a, AFR){
  F<- findF_EulerLotka(lambda=lambda, a=a, s=s, AFR=AFR)
  j<- F / b
  j[j < 0 | j > 1]<- NA
  return (j)
}

#' @rdname EulerLotka
#' @export
findA_EulerLotka<- function(lambda, b, j, s=ifelse(AFR == 1, 0, NA), AFR){
  a<- -( (b * j^2 * s^(AFR-1) * lambda - lambda^(AFR + 1)) / (lambda^AFR) )
  a[a < 0 | a > 1]<- NA
  warning("findA_EulerLotka() produce wrong results.")
  return (a)
}

#' @rdname EulerLotka
#' @export
findS_EulerLotka<- function(lambda, b, j, a, AFR){
  s<- ((lambda^AFR) / b - a * lambda^(AFR-1) / b)^(1/(AFR-1)) / (j^(2 / (AFR - 1)))
  s[s < 0 | s > 1]<- NA
  warning("findS_EulerLotka() produce wrong results.")
  return (s)
}

#' @rdname EulerLotka
#' @export
findB_EulerLotka<- function(lambda, j, a, s=a, AFR){
  F<- findF_EulerLotka(lambda=lambda, a=a, s=s, AFR=AFR)
  b<- F / j
  return (b)
}


## Rmax ----
cole<- function(x, fecundity, ageFirstBreeding, lifespan){
  return(exp(-x) + fecundity * exp(-(x * ageFirstBreeding)) - fecundity * exp(-(x * lifespan)) - 1)
}

#' Cole's Rmax
#'
#' @param fecundity 
#' @param ageFirstBreeding 
#' @param lifespan 
#'
#' @return Rmax
#' @examples Rmax(2, 2, 5)
#' @references Cole 1954. The population consequences of life history phenomena. The Quarterly Review of Biology 29 (2) pp: 103-137 
#' @export
Rmax<- function(fecundity, ageFirstBreeding, lifespan){ # values in years
  return (stats::uniroot(cole, c(0.01,30), fecundity=fecundity, ageFirstBreeding=ageFirstBreeding, lifespan=lifespan)$root)
}



## popbio matrix analysis ----
# Population matrix anaysis with popbio
#
# @param x a Leslie or Lefkovitch matrix
# @return a data.frame with the characteristics of the matrix.
eigen.analisys2df<- function(x){
  tmpRes<- popbio::eigen.analysis(x)

  ## Elasticities of survival for reproductive stages
  # Reproduction reclute individuals in the first class only (first row)
  selReproClass<- which(x[1,] != 0)
  selSurvReproClass<- which(x[-1, selReproClass] != 0)
  elasSurvRepro<- sum(tmpRes$elasticities[-1, selReproClass][selSurvReproClass]) # Survival
  
  ## Elasticities of survival for non reproductive stages
  selNonReproClass<- which(x[1,] == 0)
  selSurvNonReproClass<- which(x[-1, selNonReproClass] != 0)
  elasSurvNonRepro<- sum(tmpRes$elasticities[-1, selNonReproClass][selSurvNonReproClass]) # Survival
  
  ## Elasticities of net fecundity (for pre-breeding census matrix, fecundity * juvSurvival) or fecundity (for post-breeding census matrix)
  elasFecundity<- sum(tmpRes$elasticities[1,])
  
  tmp<- popbio::fundamental.matrix(x)
  if (is.list(tmp)){
    # Life expectancy for the first repoductive class
    matureLifeExpectancy<- as.numeric(tmp$N[selReproClass[1], selReproClass[1]])
  }else{
    matureLifeExpectancy<- NA
  }
  
  generation.time<- popbio::generation.time(x)
  net.reproductive.rate<- popbio::net.reproductive.rate(x)
  
  survCurve<- cumprod(colSums(x[-1, , drop=FALSE]))
  ## how to classify the curves??

  data.frame(elasFecundity, elasSurvRepro, elasSurvNonRepro,
             generation.time, net.reproductive.rate,
             matureLifeExpectancy, damping.ratio=tmpRes$damping.ratio)
}


