library(LHR)
# library(betaBinomial) ???
pochhammer<- function(x, n){
  return(gamma(x+n)/gamma(x))
}

size<- 5
shape1<- 1
shape2<- 9
x<- 2
## Beta-binomial distribution
choose(size,x) * beta(x+shape1, size-x+shape2) / beta(shape1, shape2) # WP
choose(shape1 + x - 1, x) * choose(shape2 + size - x - 1, size - x) / choose(shape1 + shape2 + size - 1, size) #JAGS
# LOG
lchoose(size,x) + lbeta(x+shape1, size-x+shape2) - lbeta(shape1, shape2) # WP
lchoose(shape1 + x - 1, x) + lchoose(shape2 + size - x - 1, size-x) - lchoose(shape1 + shape2 + size - 1, size) #JAGS

## Beta-negative binomial distribution
gamma(shape1 + size) / gamma(shape1) * gamma(size + x) / gamma(size) * gamma(shape2 + x) / gamma(shape2) / (gamma(x+1) * gamma(shape1 + shape2 + size) / gamma(shape1 + shape2) * gamma(size + shape1 + shape2 + x) / gamma(size + shape1 + shape2))

# LOG
log(gamma(shape1 + size) / gamma(shape1) * gamma(size + x) / gamma(size) * gamma(shape2 + x) / gamma(shape2) / (gamma(x+1) * gamma(shape1 + shape2 + size) / gamma(shape1 + shape2) * gamma(size + shape1 + shape2 + x) / gamma(size + shape1 + shape2)))
lgamma(shape1 + size) - lgamma(shape1) + lgamma(size + x) - lgamma(size) + lgamma(shape2 + x) - lgamma(shape2) - (lgamma(x+1) + lgamma(shape1 + shape2 + size) - lgamma(shape1 + shape2) + lgamma(size + shape1 + shape2 + x) - lgamma(size + shape1 + shape2))