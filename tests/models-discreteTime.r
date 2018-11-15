#Load a strategy
LH<- LH()[1,]
a=.6; s=.5;  j=.3; broods=2; b=2; fecundity=4; AFR=1

# Run discrete timemodels
replicates<- 10
tf<- 10
N0<- 5
breedFail<- .5
sexRatio<- .5
matingSystem<- "monogamy"

LHR:::mFit.t(fecundity, j, s, a, AFR, N0, replicates, tf)
LHR:::mSurvBV.t(broods, b, j, s, a, AFR, breedFail, N0, replicates, tf)
LHR:::mFitSex.t(fecundity, j, s, a, AFR, sexRatio, matingSystem, AFR, N0, replicates, tf)
LHR:::mSurvBVSex.t(broods, b, s, j, a, AFR, breedFail, sexRatio, matingSystem, N0, replicates, tf)

discretePopSim(j=j, a=a, b=b, N0=N0, replicates=replicates, tf=tf)

# mFit.t(LH, N0, replicates, tf)
# mSurvBV.t(LH, breedFail, N0, replicates, tf)
# mFitSex.t(LH, sexRatio, matingSystem, N0, replicates, tf)
# mSurvBVSex.t(LH, breedFail, sexRatio, matingSystem, N0, replicates, tf)