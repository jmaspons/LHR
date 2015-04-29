# Constructor
Environment(1)
Environment(.5, var=.1)
Environment(.5, seasonAmplitude=.3)
Environment(seasonRange=c(0,1))
Environment(.5, var=.1, seasonAmplitude=.3)
Environment(.5, var=1) # parameters out of the Beta distribution domain
Environment(.5, seasonRange=c(0,1)) # mean and range are redundant parameters on a sinusoidal function

## Seasonal pattern
seasonalPattern(Environment(.5, seasonAmplitude=1))
seasonalPattern(Environment(.5, seasonAmplitude=1), 365)
seasonalPattern(Environment(.5, seasonAmplitude=1), nSteps=3, interval=1)

seasonalPattern(Environment(.5, seasonAmplitude=1), nSteps=3, interval=1)
seasonalPattern(Environment(.5, seasonAmplitude=1), nSteps=3, interval=1, criterion="maxMean")
seasonalPattern(Environment(.5, seasonAmplitude=1), nSteps=3, interval=2, criterion="maxMean")
seasonalPattern(Environment(.5, seasonAmplitude=1), nSteps=4, interval=2, criterion="maxMean")