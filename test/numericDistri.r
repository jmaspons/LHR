gctorture(on=TRUE)
gctorture(on=FALSE)
distri<- distriBinom(2, .6)
distriC<- distriBinom(distri, .3)

res<- resC<- resS<- resP<- numeric()
for (i in 1:1000){
  res[i]<- cumsum(distriBinom(2, .6))$cump[3]
  resP[i]<- cumsum(distri * 2)$cump[3]
  resC[i]<- cumsum(distriBinom(distri, .3))$cump[3] ## Fixed
  resS[i]<- cumsum(distri + distriC)$cump[5] ## Fixed
  # print(resS[i])
}

table(res)
table(resP)
table(resC)
table(resS)

## logP
distri<- distriBinom(2, .6, log=TRUE)
distriC<- distriBinom(distri, .3, log=TRUE)

res<- resC<- resS<- resP<- numeric()
for (i in 1:1000){
  res[i]<- cumsum(distriBinom(2, .6, log=TRUE))$cump[3]
  resP[i]<- cumsum(distri * 2)$cump[3]
  resC[i]<- cumsum(distriBinom(distri, .3, log=TRUE))$cump[3] ## Fixed
  resS[i]<- cumsum(distri + distriC)$cump[5] ## Fixed
  # print(resS[i])
}

table(res)
table(resP)
table(resC)
table(resS)
