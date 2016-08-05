env<- Env(mean=1, var=0, seasonAmplitude=0, breedFail=c(0, .5, .75, 1))
sim<- Sim.numericDistri(raw = FALSE) # warning
sim<- Sim.numericDistri(raw = TRUE) # TODO FIXME: error on raw
lh<- LH()
model<- Model(lh=lh, env=env, sim=sim)
res<- run(model)
sapply(res@sim@raw, function(x) sum(x[[1]]$p))
sapply(res@sim@raw, function(x) sum(x[[2]]$p))

#Load a strategy
lh<- LH()[1,]


# Using the basic functions
a=.6; j=.3; broods=2; b=1; fecundity=broods * b;
N0<- 2
breedFail<- .5
sexRatio<- .5
tf<- 1
matingSystem<- "monogamy"

x<- tDistri(j=j, a=a, b=b, N=N0, tf=tf)
y<- tDistri(j=j, a=a, b=b, broods=broods, breedFail=breedFail, N=N0, tf=tf)

x<- tDistri(j=j, a=a, b=b, N=N0, tf=tf)
y<- tDistri(j=j, a=a, b=b, broods=broods, breedFail=0, N=N0, tf=tf)


# gctorture(on=TRUE)
mFit.distri(fecundity, j, a, N0)
mSurvBV.distri(broods, b, j, a, breedFail, N=N0)

## Check: should be identical
x<- mSurvBV.distri(broods, b, j, a, 0, N0)
y<- mSurvBV.distri(broods/2, b*2, j, a, 0, N0)
z<- mFit.distri(broods*b, j, a, N0)
plot(x); points(y); points(z)

mFitSex.distri(fecundity, j, a, sexRatio, matingSystem, N0)
mSurvBVSex.distri(broods, b, j, a, breedFail, sexRatio, matingSystem, N0) # TODO

## Multiple time steps
tf<- 20
system.time(rFit<- mFit.tdistri(fecundity, j, a, N0, tf, log=FALSE))
system.time(rFitL<- mFit.tdistri(fecundity, j, a, N0, tf, log=TRUE))
rBV<- mSurvBV.tdistri(broods, b, j, a, breedFail, N0, tf)
