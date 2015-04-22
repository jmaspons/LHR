## TODO: translate 
# see popbio package https://github.com/cstubben/popbio/wiki/Creating-projection-matrices
# LefkovitchPre(sA=0.7, sSubA=0.5, FA=1, alpha=2)
# LefkovitchPost(sA=.7, sSubA=0.5, sJ=.5, bA=2, alpha=2)
# LefkovitchPre(sA=0.7, sSubA=0.5, FA=1, alpha=1)
# LefkovitchPost(sA=.7, sSubA=0.5, sJ=.5, bA=2, alpha=1)
# lambda.leslie.matrix(matriu=LefkovitchPre(sA=0.7, sSubA=0, FA=1, alpha=2))
# findF(lambda=1.1, sA=0.6, alpha=1) #Retorna la fertilitat neta en funció de la lambda, supervivència dels adults i edat de primera reproducció
# LifeExpectancy()
##revisar:
# CalculaValorReproduccio(estrategia, numPostes=1)
# Eigen.analysis(A)
# TempsGeneraciofunction(sensitivities)
# Gen.time(A,peryear=5)
# CalculaRmax()

## PARÀMETRES ##
# alpha: Age at First Breeding

LefkovitchPre<- function(sA, sSubA=0, FA, alpha=1){
  nClasses<- alpha
  if (nClasses == 1){nClasses<- 2}
  matriu<- matrix(0, nrow=nClasses, ncol=nClasses, dimnames=list(c("J", rep("sA", nClasses-2), "A"), c("J", rep("sA", nClasses-2), "A")))
  matriu[1, nClasses]<- FA
  matriu[nClasses, nClasses]<- sA
  if (alpha > 1){
    for (i in 2:(nClasses)){
      matriu[i,i-1]<- sSubA
    }
  }else{
    matriu[1, 1]<- FA
    matriu[2, 1]<- sA
  }
  class(matriu)<- "leslie.matrix"
  return (matriu)
}

## TODO: check for errors...
LefkovitchPost<- function(sA, sSubA, sJ, bA, alpha=1){
  nClasses<- alpha + 1
  matriu<- matrix(0, nrow=nClasses, ncol=nClasses, dimnames=list(c("J", rep("sA", nClasses-2), "A"), c("J", rep("sA", nClasses-2), "A")))
  matriu[1, nClasses]<- bA
  matriu[nClasses, nClasses]<- sA
  matriu[2, 1]<- sJ
  if (nClasses > 2){
    for (i in 3:(nClasses)){
      matriu[i,i-1]<- sSubA
    }
  }else{
    matriu[1, 1]<- 0
    matriu[2, 1]<- sJ
  }
  class(matriu)<- "leslie.matrix"
  return (matriu)
}

lambda.leslie.matrix<- function(mat){ # Ted J. Case. "An Illustrated Guide to Theoretical Ecology" pàg. 64-73
  return (max(Re(eigen(mat, only.values=TRUE)$values)))
}


##Funció d'Euler-Lotka per trobar la fertilitat neta (F) a partir de lambda i sA
# http://en.wikipedia.org/wiki/Euler–Lotka_equation
#http://darwin.eeb.uconn.edu/eeb310/lecture-notes/spotted-owl/node3.html

#(l_{alpha}*F)*lambda^-alpha sumatori(from x=alpha to infinite){lambda^(-x+alpha)*s^(x-alpha)} = 1
#(l_{alpha}*F)*(1/(1-(s/lambda))) = lambda^alpha
#lambda^alpha*(1-(s/lambda)) = (l_{alpha}*F)
##Maxima
# solve(lambda^alpha*(1-(s/lambda))=l*f, f)
# \begin{eqnarray*}\left[ f={{-\left(s\,\lambda^{\alpha}-\lambda^{\alpha+1}\right)}\over{l\,\lambda}} \right] \end{eqnarray*}
# f = (-(sA * lamda^alpha - lambda^(alpha+1))) / (l_{alpha} * lambda)
# l_{alpha} = supervivència fins a l'edat reproductora alpha -> (sJ*sA^(alpha-1))
# f = annual reproductive rate, i.e., the average number of fledged offspring per individual ->b
##Test de la funció trobaF()
# lambda<- 1.123456
# sA<-.9
# alpha<- 2
# F<- trobaF(lambda, sA, alpha)
# mat<- construeixMatriuLefkovitchPre(sA, sSubA=sA, F, alpha)
# if (calculaLambda(mat) == lambda) {cat("OK!, les lambdes coincideixen i F=", F, "\n")}else{cat("KO! això falla\n")}

findF<- function(lambda, adultSurv, alpha){ #Fertilitat neta en funció dels paràmetres demogràfics
  F<- (-(adultSurv*lambda^alpha - lambda^(alpha+1))) / (adultSurv^(alpha-1)*lambda)
  return (F)
}

# lifeExpectancy<- - 1 / log(GdemoLH$sA) ## SRC: http://www.countrysideinfo.co.uk/bird_lifespan.htm

##Funció d'Euler-Lotka per trobar la fertilitat neta (F) a partir de lambda i sA
# http://en.wikipedia.org/wiki/Euler–Lotka_equation
#http://darwin.eeb.uconn.edu/eeb310/lecture-notes/spotted-owl/node3.html

#(l_{alpha}*F)*lambda^-alpha sumatori(from x=alpha to infinite){lambda^(-x+alpha)*s^(x-alpha)} = 1
#(l_{alpha}*F)*(1/(1-(s/lambda))) = lambda^alpha
#lambda^alpha*(1-(s/lambda)) = (l_{alpha}*F)
##Maxima
# solve(lambda^alpha*(1-(s/lambda))=l*f, f)
# \begin{eqnarray*}\left[ f={{-\left(s\,\lambda^{\alpha}-\lambda^{\alpha+1}\right)}\over{l\,\lambda}} \right] \end{eqnarray*}
# f = (-(sA * lamda^alpha - lambda^(alpha+1))) / (l_{alpha} * lambda)
# l_{alpha} = supervivència fins a l'edat reproductora alpha -> (sJ*sA^(alpha-1))
# f = annual reproductive rate, i.e., the average number of fledged offspring per individual ->b

CalculaValorReproduccio<- function(estrategia, numPostes=1){ ##Bókony et al 2009. Stress response and the value of reproduction
  anysReproduccio<- 1 / (1 - estrategia$sA)
  valorReproduccio<- log10((estrategia$b /numPostes) / (estrategia$b * anysReproduccio))
  return (valorReproduccio)
}

Eigen.analysis<- function(A){
  ev <- eigen(A)
  lmax <- which(Re(ev$values)==max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[,lmax]))
  V <- Conj(solve(W))
  v <- abs(Re(V[lmax,]))
  
  s <- v%o%w
  s[A == 0] <- 0
  class(s) <- "leslie.matrix"
  e <- s*A/lambda
  rho <- lambda/abs(Re(ev$values[2]))
  
  
  eigen.analysis <- list(lambda1=lambda, rho=rho, sensitivities=s,
                         elasticities=e, stable.age=w/sum(w),
                         repro.value=v/v[1])
  
  eigen.analysis
  
}

TempsGeneracio<- function(sensitivities){##Article matrius multiestats
  return (1 / sensitivities[1, ncol(sensitivities)])
}

Gen.time<- function(A, peryear=5){ ##Temps de generació  #(Ted J. Case. "An Illustrated Guide to Theoretical Ecology" pàg 89-90 #mirar el peryear!!
  ro <- calc.ro(A)
  ea <- eigen.analysis(A)
  T <- peryear*log(ro)/log(ea$lambda)
  
  T
}

cole<- function(x, fecundity, ageFirstBreeding, lifespan){
  return(exp(-x) + fecundity * exp(-(x * ageFirstBreeding)) - fecundity * exp(-(x * lifespan)) - 1)
}

Rmax<- function(fecundity, ageFirstBreeding, lifespan, cole){ # values in years
  return (uniroot(cole, c(0.01,30), fecundity=fecundity, ageFirstBreeding=ageFirstBreeding, lifespan=lifespan)$root)
}

# F <- Fecundity
# A <- d$AgeFirstBreeding/12
# L <- d$Lifespan
# cole<- function(x) exp(-x) + d2$F[i]*exp(-(x*(d2$A[i]))) - d2$F[i]*exp(-(x*(d2$L[i]))) - 1
# #craete an empty vector
# R_max <- c(1:length(d2$F))
# #fill the vector with the Rmax values.
# for (i in 1:length(d2$F)){R_max[i] <- uniroot(cole, c(0.01,30))$root}

