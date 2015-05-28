env<- Env()
env<- Env(env[env$var == 0 & env$seasonAmplitude == 0,]) #stable environment
model<- Model(lh=LH(), env=env, sim=Sim.numericDistri())
res<- run(model)
sapply(res@sim@raw, function(x) sum(x[[1]]$p))
sapply(res@sim@raw, function(x) sum(x[[2]]$p))

#Load a strategy
lh<- LH()[1,]


# Using the basic functions
a=.6; j=.3; broods=2; b=2; fecundity=4;
N0<- 5
breedFail<- .5
sexRatio<- .5
matingSystem<- "monogamy"

mFit.trans(fecundity, j, a, N0)
mSurvBV.trans(broods, b, j, a, breedFail, N0)
mFitSex.trans(fecundity, j, a, sexRatio, matingSystem, N0)
mSurvBVSex.trans(broods, b, j, a, breedFail, sexRatio, matingSystem, N0)

