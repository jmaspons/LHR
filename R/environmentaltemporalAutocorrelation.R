## TEMPORAL AUTOCORRELATED ENVIRONMENT
## FIXME for testing pourposes. Keep it without errors to avoid installation fails
# library(fArma)
# armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100, innov = NULL, n.start = 100, start.innov = NULL, rand.gen = rnorm, rseed = NULL, addControl = FALSE, ...)

# arima.sim(model=list(ar=c(), ma=c(), order=c(length(ar), int>0 degree of differencing (no stationary), length(ma)), n=, rand.gen=rbeta, shape1=, shape2=)
# arima.sim(model, n, rand.gen = rnorm, innov = rand.gen(n, ...), n.start = NA, start.innov = rand.gen(n.start, ...), ...)
# require(truncnorm)

  
# ar<- c(.5, .5, -.6, -.4)
# f<- function(x, ar){
#   ans<- 1
#   for (i in 1: length(ar)){
#     ans<- ans + x * ar[i]^i
# # print(ans)
#   }
#   return (ans)
# }
# f(seq(-5, 10, .05), ar)
# curve(f(x, ar), -10, 10)
# 
# minroots <- min(Mod(polyroot(c(1, -ar))))
# if (minroots <= 1) 
#   cat("'ar' part of model is not stationary")
# ma<- c(1, .5)
# mean<- 0.6
# cv<- 0.1
# ar<- c(.5, -0.4)
# plot(arima.sim(model=list(ar=c(ar), ma=c(ma), order=c(length(ar), 0, length(ma))), n=50, n.start=1000, rand.gen=rbeta, shape1=2, shape2=2))
# plot(arima.sim(model=list(ar=c(ar), ma=c(ma), order=c(length(ar), 0, length(ma))), n=50, rand.gen=rtruncnorm, n.start=6, start.innov=rep(0.6, 50), a=0, b=1,mean=mean, sd=cv * mean), ylab="Parameter")
# range(arima.sim(model=list(ar=c(ar), ma=c(ma), order=c(length(ar), 0, length(ma))), n=200, rand.gen=rbeta, shape1=2, shape2=2))
