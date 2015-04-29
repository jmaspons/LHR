## SEASONAL EVIRONMENT
# modify a parameter following a seasonal pattern
par.seasonality<- function(broods, breedInterval, mean, amplitude, years=1, 
                           criterion=c("maxFirst", "maxMean"), resolution=c("months", "days")){
  seasons<- seasonality(years=years, mean=mean, amplitude=amplitude, resolution=resolution[1])
  seasons<- par.modifier(pattern=seasons, events=broods, interval=breedInterval, criterion=criterion[1])
  seasons<- sort(seasons, decreasing=TRUE)

  return (seasons)
}

# Returns a vector following a sinusoidal pattern. range = mean +- amplitude / 2  
seasonality<- function(years, mean, amplitude, resolution=c("months", "days")){
  switch (resolution[1],
      months = t<- seq(0, 2*pi*years - 2*pi/12, by=2*pi/12),
      days = t<- seq(0, 2*pi*years - 2*pi/365, by=2*pi/365)
  )

  return (sin(t)*amplitude/2 + mean) # max amplitude=1 & mean=0.5
}

# Parameter modification factor (par = par * par.seasonality())
# pattern vector with the annual pattern
# events events per year
# interval number of units between events (units are the same as in pattern)
## TODO criterion: optimize the timing of reproduction according to interbreeding intervals or some reasonable rules.
par.modifier<- function(pattern, events, interval, criterion=c("maxFirst", "maxMean")){
  n<- length(pattern)
  if (n < events * interval) stop("Pattern too short for the events number and interval parameters")
  
  if (criterion[1] == "maxFirst"){
    firstDate<- which.max(pattern)
    dates<- seq(firstDate, firstDate + events * interval, length=events)
    if (firstDate + events * interval > n){
      nextPeriod<- which(firstDate + 1:events * interval > n)
      dates[nextPeriod]<- dates[nextPeriod] - n
    }
    
  }else if (criterion[1] == "maxMean"){
    seasonLength<- events * interval
    dates<- seq(1, seasonLength, by=interval)
    ##TODO Optimization
#     ?optimize(function(i, pattern, dates) sum(pattern[dates + round(i)]), interval=0:n, maximum=TRUE tol=1)
    pattern<- rep(pattern,2)
    objective<- -Inf
    for (i in 0:n){
      newObjective<- sum(pattern[dates + i])
      if (objective < newObjective){
	objective<- newObjective
	bestFirstDate<- i
      }
    }
    dates<- dates + bestFirstDate
  }
  
  return (pattern[dates])
}

# Calculates the mean and amplitude parameters of seasonality() from the range
seasonAmplitude<- function(seasonRange){
  mean<- (min(seasonRange) + max(seasonRange))/2
  seasonAmplitude<- max(seasonRange) - min(seasonRange)
  
  return (data.frame(mean, seasonAmplitude))
}

seasonRange<- function(mean, seasonAmplitude){
  return (data.frame(min=mean - amplitude, max=mean + amplitude))
}