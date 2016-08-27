context("Euler-Lotka")

test_that("find deterministic parameters f(lambda, ...) ", {
  a=.7; s=a; j=.3; b=5; bj=b * j; AFR=4;

  matPre<- LefkovitchPre(a=a, s=s, bj=bj, AFR=AFR)
  matPost<- LefkovitchPost(a=a, s=s, j=j, b=b, AFR=AFR)
  lambdaPre<- lambda(matPre)
  lambdaPost<- lambda(matPost)
  
  expect_equal(bj - findF_EulerLotka(lambda=lambdaPre, a=a, AFR=AFR), 0) # Corresponds to a pre-breding census matrix
  # bj - findF_EulerLotka(lambda=lambdaPost, a=a, AFR=AFR)
  
  expect_equal(j - findJ_EulerLotka(lambda=lambdaPre, b=b, a=a, AFR=AFR), 0) # Corresponds to a pre-breding census matrix
  # j - findJ_EulerLotka(lambda=lambdaPost, b=b, a=a, AFR=AFR)
  
  expect_equal(b - findB_EulerLotka(lambda=lambdaPre, j=j, a=a, AFR=AFR), 0) # Corresponds to a pre-breding census matrix
  # b - findB_EulerLotka(lambda=lambdaPost,j=j, a=a, AFR=AFR)
  
  ## ERRORS
  # expect_equal(a - findA_EulerLotka(lambda=lambdaPre, b=b, j=j, AFR=AFR), 0)
  # expect_equal(a - findA_EulerLotka(lambda=lambdaPre, b=b, j=j, AFR=AFR), NA_real_) # Error a > 1
  # a - findA_EulerLotka(lambda=lambdaPost, b=b, j=j, AFR=AFR)
  
  
  ## sampleLH
  lambda=seq(.8, 2, by=0.2); broods=2^(0:2); b=c(1, seq(2, 20, by=2)); j=seq(0.2, 0.8, by=0.2); a=seq(0.3, 0.9, by=0.2); AFR=1
  parsL<- sampleLH(free="lambda", census="pre-breeding")
  parsJ<- sampleLH(free="j", AFR=1:5)
  # parsA<- sampleLH(free="a")

  parsL$lambdaMat<- apply(parsL, MARGIN=1, function(x){
    mat<- with(as.list(x), LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
    lambda(mat)
  })

  parsJ$lambdaMat<- apply(parsJ, MARGIN=1, function(x){
    mat<- with(as.list(x), LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
    lambda(mat)
  })
  
  # parsA$lambdaMat<- apply(parsA, MARGIN=1, function(x){
  #   mat<- with(as.list(x), LefkovitchPre(a=a, s=s, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
  #   lambda(mat)
  # })
  
  parsL$errLambda<- abs(parsL$lambda - parsL$lambdaMat)
  parsJ$errLambda<- abs(parsJ$lambda - parsJ$lambdaMat)
  # parsA$errLambda<- abs(parsA$lambda - parsA$lambdaMat)
  
  expect_lt(max(parsL$errLambda), 1e-10)
  expect_lt(max(parsJ$errLambda), 1e-10)
  # expect_lt(max(parsA$errLambda), 1e-10)
  
  # plot(parsJ, col=(abs(parsJ$errLambda) >= .1) +1) # In red wrong estimates of lambda
  # plot(parsJ[abs(parsJ$errLambda) >= .01,], col="red")
  # plot(parsJ[abs(parsJ$errLambda) < .01,], col="blue")
  # hist(parsJ$errLambda[abs(parsJ$errLambda) < .01])
  
  # abline(0,1, col="red")

  err<- apply(parsL, 1, function(x){
    # Net fecundity: F = fecundity * juvenile survival
    errF<- with(as.list(x), fecundity * j - findF_EulerLotka(lambda=lambda, a=a, AFR=AFR)) #Error on findF
    # errA<- with(as.list(x), a - findA_EulerLotka(lambda=lambda, b=fecundity, j=j, AFR=AFR)) # Error on findA returns 4 values instead of 1
    errJ<- with(as.list(x), j - findJ_EulerLotka(lambda=lambda, b=fecundity, a=a, AFR=AFR))
    errB<- with(as.list(x), fecundity - findB_EulerLotka(lambda=lambda, j=j, a=a, AFR=AFR))
    
    c(errF=errF, errJ=errJ, errB=errB) # errA=errA, 
  })
  
  # matplot(t(err))

  expect_lt(max(err[1,]), 1e-10)
  expect_lt(max(err[2,]), 1e-10)
  expect_lt(max(err[3,]), 1e-10)
})
