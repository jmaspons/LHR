a=.9; s=a; j=.1; b=20; bj=b * j; AFR=5
matPre<- LefkovitchPre(a=a, s=s, bj=bj, AFR=AFR)
matPost<- LefkovitchPost(a=a, s=s, j=j, b=b, AFR=AFR)
lambdaPre<- lambda(matPre)
lambdaPost<- lambda(matPost)

bj - findF_EulerLotka(lambda=lambdaPre, a=a, AFR=AFR) # Corresponds to a pre-breding census matrix
bj - findF_EulerLotka(lambda=lambdaPost, a=a, AFR=AFR)

j - findJ_EulerLotka(lambda=lambdaPre, b=b, a=a, AFR=AFR) # Corresponds to a pre-breding census matrix
j - findJ_EulerLotka(lambda=lambdaPost, b=b, a=a, AFR=AFR)

b - findB_EulerLotka(lambda=lambdaPre, j=j, a=a, AFR=AFR) # Corresponds to a pre-breding census matrix
b - findB_EulerLotka(lambda=lambdaPost,j=j, a=a, AFR=AFR)

a - findA_EulerLotka(lambda=lambdaPre, b=b, j=j, AFR=AFR) # Error
a - findA_EulerLotka(lambda=lambdaPost, b=b, j=j, AFR=AFR)


## sampleLH
lambda=seq(.8, 2, by=0.2); broods=2^(0:2); b=c(1, seq(2, 20, by=2)); j=seq(0.2, 0.8, by=0.2); a=seq(0.3, 0.9, by=0.2); AFR=1
parsL<- sampleLH(free="lambda", census="pre-breeding")

parsJ<- sampleLH(free="j", AFR=1:5)
jOut<- parsJ[is.na(parsJ$j),] ## FILTERED inside the sampleLH function
plot(jOut[,])

parsJ<- parsJ[!is.na(parsJ$j),]
parsJ$lambdaMat<- NA
for (i in 1:nrow(parsJ)){
  mat<- with(parsJ[i,], LefkovitchPre(a=a, s=a, bj=j * fecundity, AFR=AFR)) # subadult survival equal to adult survival
  parsJ$lambdaMat[i]<- lambda(mat)
}
parsJ$errLambda<- abs(parsJ$lambda - parsJ$lambdaMat)
parsJ<- parsJ[order(abs(parsJ$errLambda), decreasing=TRUE),]

plot(parsJ, col=(abs(parsJ$errLambda) >= .1) +1) # In red wrong estimates of lambda
plot(parsJ[abs(parsJ$errLambda) >= .01,], col="red")
plot(parsJ[abs(parsJ$errLambda) < .01,], col="blue")
hist(parsJ$errLambda[abs(parsJ$errLambda) < .01])

abline(0,1, col="red")
# Errors!! Not vectorizable??
for (i in 1:nrow(parsL)){
  errF[i]<- parsL$fecundity[i]  - with(parsL[i,], findF(lambda=lambda, a=a, AFR=AFR))
  errj[i]<- parsL$j[i] - with(parsL[i,], findJ(lambda=lambda, b=b, a=a, AFR=AFR))
}
plot(data.frame(errF, parsL))

with(parsL, findJ(lambda=labmda, b=b, a=a, AFR=AFR))
with(parsL, findB(lambda=lambda, j=j, a=a, AFR=AFR))