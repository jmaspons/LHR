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
#'
#' @exportClass leslieMatrix
setOldClass("leslieMatrix")

#' @describeIn lefkovitch
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

#' @describeIn lefkovitch
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


#' Lambda
#'
#' @param mat 
#'
#' @return
#' @export
#'
#' @examples
lambda.leslieMatrix<- function(mat){ # Ted J. Case. "An Illustrated Guide to Theoretical Ecology" pàg. 64-73
  return (max(Re(eigen(mat, only.values=TRUE)$values)))
}

# Use eigen.analysis function to analyse matrices. It is imported from popbio package
# popbio::eigen.analysis(mat)

## Inverse eigenvalue problem: find a, j, b or F from lambda and the rest of parameters ----
# Not so easy, some parameters fall outside the domain (e.g. survival <0 or >1): 

## CORRESPONDS TO A PRE-BREEDING CENSUS MATRIX (see test/Euler-Lotka.r)
# Euler-Lotka equation fails for some cases with no clear pattern (see test/Euler-Lotka.r)
# http://en.wikipedia.org/wiki/Euler–Lotka_equation
# http://darwin.eeb.uconn.edu/eeb310/lecture-notes/spotted-owl/node3.html
# see Euler-Lotka.txt
# alpha= AFR
#(l_{alpha}*F)*lambda^-alpha sumatori(from x=alpha to infinite){lambda^(-x+alpha)*s^(x-alpha)} = 1
#(l_{alpha}*F)*(1/(1-(a/lambda))) = lambda^alpha
# lambda^alpha*(1-(a/lambda)) = (l_{alpha}*F)
# l_{alpha}:  survival until reproductive class -> j * a^(alpha-1)
# lambda^alpha*(1-(a/lambda)) = j * a^(alpha-1) * F

## Maxima
# solve(lambda^alpha * (1 - (s / lambda)) = j * s^(alpha-1) * F, F);
# solve(lambda^alpha * (1 - (a / lambda)) = j * a^(alpha-1) * j * b, a);
# solve(lambda^alpha * (1 - (a / lambda)) = j * a^(alpha-1) * j * b, j);
# solve(lambda^alpha * (1 - (a / lambda)) = j * a^(alpha-1) * j * b, b);
# solve(lambda^alpha * (1 - (a / lambda)) = j * a^(alpha-1) * j * b, alpha);

#' Euler-Lotka equations
#'
#' @name EulerLotka
#' @param lambda 
#' @param a 
#' @param AFR
NULL

#' @describeIn EulerLotka
#' @export
findF_EulerLotka<- function(lambda, a, AFR){
  F<- -(a * lambda^AFR - lambda^(AFR+1)) / (a^(AFR-1) * lambda)
  F[F < 0]<- NA
  return (F)
}

#' @describeIn EulerLotka
#' @export
findJ_EulerLotka<- function(lambda, b, a, AFR){
  F<- findF_EulerLotka(lambda, a, AFR)
  j<- F / b
  j[j < 0 | j > 1]<- NA
  return (j)
}

#' @describeIn EulerLotka
#' @export
findA_EulerLotka<- function(lambda, b, j, AFR){
  a<- -(b * j^2 * a^(AFR-1) * lambda - lambda^(AFR + 1)) / (lambda^AFR)
  a[a < 0 | a > 1]<- NA
  return (a)
}

#' @describeIn EulerLotka
#' @export
findB_EulerLotka<- function(lambda, j, a, AFR){
  F<- findF_EulerLotka(lambda, a, AFR)
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
  return (uniroot(cole, c(0.01,30), fecundity=fecundity, ageFirstBreeding=ageFirstBreeding, lifespan=lifespan)$root)
}



## TODO: check and reference ----

# lifeExpectancy<- - 1 / log(GdemoLH$sA) ## SRC: http://www.countrysideinfo.co.uk/bird_lifespan.htm

# TempsGeneracio<- function(sensitivities){##Article matrius multiestats
#   return (1 / sensitivities[1, ncol(sensitivities)])
# }

# Gen.time<- function(A, peryear=5){ ##T  #(Ted J. Case. "An Illustrated Guide to Theoretical Ecology" pàg 89-90 #mirar el peryear!!
#   ro <- calc.ro(A)
#   ea <- eigen.analysis(A)
#   T <- peryear*log(ro)/log(ea$lambda)
#   
#   T
# }

