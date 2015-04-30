# Constructor
env<- Env()
Env(.5, var=.1)
Env(.5, seasonAmplitude=.3)
Env(seasonRange=c(0,1))
Env(.5, var=.1, seasonAmplitude=.3)
Env(.5, var=1) # parameters out of the Beta distribution domain
Env(.5, seasonRange=c(0,1)) # mean and range are redundant parameters on a sinusoidal function

## Seasonal pattern
seasonalPattern(env)
seasonalPattern(env, resolution=365)
seasonalPattern(env, cicles=2)

seasonOptimCal(env)
seasonOptimCal(env, resolution=12, nSteps=3, interval=1, criterion="maxFirst")
seasonOptimCal(env, resolution=12, nSteps=3, interval=3, criterion="maxFirst")
seasonOptimCal(env, resolution=12, nSteps=3, interval=1, criterion="maxMean")
seasonOptimCal(env, resolution=12, nSteps=3, interval=3, criterion="maxMean")

# # nSteps<- c(1,3, 5)
# # interval<- c(1,3)
# # criterion<- c("maxFirst", "maxMean")
# # cals<- expand.grid(nSteps, interval, criterion, stringsAsFactors=FALSE)
# # names(cals)<- c("nSteps", "interval", "criterion")
# # cals$resolution<- 12
# #
# # envCals<- merge(S3Part(env), cals)
# # seasonOptimCal(envCals)
