# Geometric mean
# gr<- r(pop) | lambda(pop) # grow rate
## TODO: check paper about geometric mean:
# Habib, E. A. (2012). Geometric mean for negative and zero values. International Journal of Research and Reviews in Applied Sciences, 11, 419-432.
Gmean<- function(gr){
  if (all(gr != 0)){
    res<- prod(gr)^(1/length(gr)) # buffer overflow for large numbers
    if (is.finite(res))  return (res)
    res<- exp(mean(log(gr))) # NaN for negative numbers
    if (!is.na(res)) return (res)
  }else{
    return (G(mean(gr), var(as.numeric(gr))))
  }
}


# Geometric mean approximation from mean and variance
# Stearns, S. C. (2000). Daniel Bernoulli (1738): evolution and economics under risk. Journal of Biosciences, 25(3), 221–228
G<- function(...){
  UseMethod("G")
}

G.numeric<- function(mu, sigma2){
  return (mu - sigma2 / (2 * mu))
}

G.discretePopSim<- function(pop, growRate=c("r", "lambda")[1]){
  Gres<- switch(growRate,
         r={rPop<- as.numeric(r(pop)); # avoid var() returns a var/cov matrix
            G(mean(rPop), var(rPop))},
         lambda={lambdaPop<- as.numeric(lambda(pop));
            GLambda<- G(mean(lambdaPop), var(lambdaPop))})
  return (Gres)
} 