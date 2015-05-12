#Load a strategy
LH<- LH()[1,]
a=.6; j=.3; broods=2; b=2; fecundity=4;

# Run discrete timemodels
replicates<- 10
tf<- 10
N0<- 5
breedFail<- .5
sexRatio<- .5
matingSystem<- "monogamy"

mFit.t(fecundity, j, a, N0, replicates, tf)
mSurvBV.t(broods, b, j, a, breedFail, N0, replicates, tf)
mFitSex.t(fecundity, j, a, sexRatio, matingSystem, N0, replicates, tf)
mSurvBVSex.t(broods, b, j, a, breedFail, sexRatio, matingSystem, N0, replicates, tf)

discretePopSim(j=j, a=a, b=b, N0=N0, replicates=replicates, tf=tf, maxN=maxN)

mFit.t(LH, N0, replicates, tf)
mSurvBV.t(LH, breedFail, N0, replicates, tf)
mFitSex.t(LH, sexRatio, matingSystem, N0, replicates, tf)
mSurvBVSex.t(LH, breedFail, sexRatio, matingSystem, N0, replicates, tf)