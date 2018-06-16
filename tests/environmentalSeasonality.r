nAnys<- 2

x<- seq(0, 2*pi*nAnys - 2*pi/12, by=2*pi/12) # mesos de l'any (2*pi = 1 any), - 2*pi/12 pq comença a 0 i si no hi ha nAnys + 1 mes

amplitude<- c(1, .5, .1)
mean<- c(0.5, 0.75, .95)
s<- data.frame(seasonality=numeric(length(amplitude) * nAnys * 12), month=rep(1:(length(x)), times=length(amplitude)), amplitude=rep(amplitude, each= nAnys * 12), mean=rep(mean, each=nAnys * 12))
for (i in 1:nrow(s)){
  s$seasonality[i]<- sin(x[s$month[i]]) * s$amplitude[i] / 2 + s$mean[i]
}

library(lattice)
xyplot(seasonality ~ month, group= amplitude, data=s, ylab="Seasonality", xlab="Month", type="l", main="Seasonality",
       par.settings=list(par.main.text=list(cex=4), par.xlab.text=list(cex=3), par.ylab.text=list(cex=3)) , scales=list(cex=3), lwd=2)



s<- sin(x) * 1 / 2 + 0.5 # Valors entre 0 i 1

y1<- sin(x) + rnorm(length(x),0,0.1)
y2<- (sin(x) + rnorm(length(x),0,0.1))/2 + 0.3
y3<- (sin(x) + rnorm(length(x),0,0.1))/10 + 0.5
y4<- sin(x)

plot(y1, type="l", ylab="Produccio neta", xlab="Temps (mesos)")
points(y2, type="l", col="green")
points(y3, type="l", col="red")
points(y4, type="l", col="blue")

abline(h=0, lty=2)
abline(h=0.3, col="green", lty=2)
abline(h=0.5, col="red", lty=2)



# Values between 0 and 1 to modify aprameters
seasonPars<- seasonal.range(min=c(.1,.5,1),max=1)
env<- list()
for (i in 1:nrow(seasonPars)){
  env[[i]]<- seasonality(years=2, seasonPars$mean[i], seasonPars$amplitude[i])
}
env<- data.frame(env)
names(env)<- seasonPars$mean

# Stochasticity
y1<- env + rnorm(length(x),0,0.1)
y2<- (sin(x) + rnorm(length(x),0,0.1))/2 + 0.3
y3<- (sin(x) + rnorm(length(x),0,0.1))/10 + 0.5
y4<- sin(x)

plot(y1, type="l", ylab="Produccio neta", xlab="Temps (mesos)")
points(y2, type="l", col="green")
points(y3, type="l", col="red")
points(y4, type="l", col="blue")

abline(h=0.3, col="green", lty=2)
abline(h=0.5, col="red", lty=2)
abline(h=0, lty=2)



