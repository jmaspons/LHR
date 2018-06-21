context("Deterministic models")

lambda<- seq(.9, 1.5, length.out=3)
b<- 3^(0:2)
a<- seq(0.4, 0.95, length.out=3)
j<- seq(0.4, 0.7, length.out=3)
AFR<- 2^(0:2)
pars<- expand.grid(b=b, j=j, a=a, AFR=AFR) # free lambda
pars$s<- pars$a * .8




test_that("Euler-Lotka inverse eigenvalue", {
  matPost<- apply(pars, 1, function(x){
      mat<- with(as.list(x), LefkovitchPost(a=a, s=s, j=j, b=b, AFR=AFR)) 
    })

  matPre<- apply(pars, 1, function(x){
      mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=b*j, AFR=AFR))
    })
  
  ## Euler-Lotka equations corresponds to a pre-breding census matrix
  pars$lambdaPre<- sapply(matPre, lambda)
  
  res<- apply(pars, 1, function(x){
    with(as.list(x), expect_equal(findF_EulerLotka(lambda=lambdaPre, a=a, s=s, AFR=AFR), b*j))
  })
  res<- apply(pars, 1, function(x){
    with(as.list(x), expect_equal(findJ_EulerLotka(lambda=lambdaPre, b=b, a=a, s=s, AFR=AFR), j))
  })
  res<- apply(pars, 1, function(x){
    with(as.list(x), expect_equal(findB_EulerLotka(lambda=lambdaPre, j=j, a=a, s=s, AFR=AFR), b))
  })
  
  ## TODO: ERRORS finding a & s
  # expect_equal(a - findA_EulerLotka(lambda=lambdaPre, b=b, j=j, s=s, AFR=AFR), 0)
  # expect_equal(a - findA_EulerLotka(lambda=lambdaPre, b=b, j=j, s=s, AFR=AFR), NA_real_) # Error a > 1
  # a - findA_EulerLotka(lambda=lambdaPost, b=b, j=j, AFR=AFR)
  
  # pars$aEL<- apply(pars, 1, function(x){
  #   # with(as.list(x), expect_equal(findA_EulerLotka(lambda=lambdaPre, b=b, j=j, s=s, AFR=AFR), a))
  #   with(as.list(x), findA_EulerLotka(lambda=lambdaPre, b=b, j=j, s=s, AFR=AFR))
  # })
  # selS<- pars$AFR > 1
  # pars$sEL[selS]<- apply(pars[selS,], 1, function(x){
  #   # with(as.list(x), expect_equal(findS_EulerLotka(lambda=lambdaPre, b=b, j=j, a=a, AFR=AFR), s))
  #   with(as.list(x), findS_EulerLotka(lambda=lambdaPre, b=b, j=j, a=a, AFR=AFR))
  # })
  # 
  # pars[, c("a", "aEL", "s", "sEL")]
  # pars[selS, c("a", "aEL", "s", "sEL")]
  # table(selErrA<- pars$a == pars$aEL, useNA="alw"); selErrA<- !selErrA | is.na(selErrA)
  # table(selErrS<- pars$s == pars$sEL, useNA="alw"); selErrS<- !selErrS | is.na(selErrS)
  # 
  # table(selErrA & selErrS)
  # 
  # View(pars[selErrA, ])
  # pars[selErrA, ]
  
})


test_that("Sample LH", {
  broods<- 2^(0:2)
  
  parsL<- sampleLH(free="lambda", method="regular")
  # warnings: Some parameter combinations produce NAs and are discarded (eg. probabilities > 1)
  suppressWarnings(parsJ<- sampleLH(free="j", AFR=AFR, method="regular"))
  suppressWarnings(parsA<- sampleLH(free="a", method="regular"))
  
  parsL$lambdaMat<- apply(parsL[, sapply(parsL, function(x) inherits(x, "numeric"))], MARGIN=1, function(x){
    mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=j * fecundity, AFR=AFR))
    lambda(mat)
  })
  
  parsJ$lambdaMat<- apply(parsJ[, sapply(parsJ, function(x) inherits(x, "numeric"))], MARGIN=1, function(x){
    mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=j * fecundity, AFR=AFR))
    lambda(mat)
  })
  
  parsA$lambdaMat<- apply(parsA[, sapply(parsA, function(x) inherits(x, "numeric"))], MARGIN=1, function(x){
    mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=j * fecundity, AFR=AFR))
    lambda(mat)
  })
  
  parsL$errLambda<- abs(parsL$lambda - parsL$lambdaMat)
  parsJ$errLambda<- abs(parsJ$lambda - parsJ$lambdaMat)
  parsA$errLambda<- abs(parsA$lambda - parsA$lambdaMat)
  
  expect_lt(max(parsL$errLambda), .Machine$double.eps)
  expect_lt(max(parsJ$errLambda), .Machine$double.eps * 10)
  # expect_lt(max(parsA$errLambda), .Machine$double.eps) # ERRORS in findA_EulerLotka()

  # plot(parsA, col=(abs(parsA$errLambda) >= .1) +1) # In red wrong estimates of lambda
  # plot(parsA[abs(parsA$errLambda) >= .01,], col="red")
  # plot(parsA[abs(parsA$errLambda) < .01,], col="blue")
  # hist(parsA$errLambda[abs(parsA$errLambda) < .01])
  
  # abline(0,1, col="red")
  
  err<- apply(parsL[, sapply(parsL, function(x) inherits(x, "numeric"))], 1, function(x){
    # Net fecundity: F = fecundity * juvenile survival
    errF<- with(as.list(x), fecundity * j - findF_EulerLotka(lambda=lambda, a=a, s=s, AFR=AFR)) #Error on findF
    # errA<- with(as.list(x), a - findA_EulerLotka(lambda=lambda, b=fecundity, j=j, AFR=AFR)) # Error on findA returns 4 values instead of 1
    errJ<- with(as.list(x), j - findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, s=s, AFR=AFR))
    errB<- with(as.list(x), fecundity - findB_EulerLotka(lambda=lambda, j=j, a=a, s=s, AFR=AFR))
    
    c(errF=errF, errJ=errJ, errB=errB) # errA=errA, 
  })
  
  # matplot(t(err))
  
  err<- abs(err)
  expect_equal(max(err[1,]), 0)
  expect_equal(max(err[2,]), 0)
  expect_equal(max(err[3,]), 0)
})
