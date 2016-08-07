#' LHR.
#'
#' Implements three different types of stochastic demographic models to explore the influence
#'  of different mechanisms on the population dynamics.
#'  
#' \enumerate{
#'   \item Discrete time random simulations whith diferent underlaying probabilities distributions.
#'   \item Compound probability distribution of N_t+1 for a given N calculated numerically following the desired stochastic structure.
#'   \item IBM implemented using the stochastic simulation algorithm or Gillespie algorithm from the \link[adaptivetau]{adaptivetau} package.
#' }
#' @name LHR
#' @details 
#' The discrete random simulations and the compound probability distributions are two different approaches to explore the same underlaying model (Melbourne and Hastings 2008). The implemented underlaying structures combines individual mortality for juveniles and adults, complete failure of the reproductive attempt (e.g. nest predation) and stochasticity in the sex ratio.
#' Adult mortality + offspring mortality
#' N_t+1 = B(n=N_t * fecundity, p=juvSurv) + B(n=N_t, p=adultSurv)
#' 
#' Adult mortality + Brood mortality + offspring mortality
#' N_t+1 = B(n=B(n=N_t * broods, p=1-breedFail) * clutch, p=juvSurv) + B(n=N_t, p=adultSurv)
#' 
#' Adult mortality + offspring mortality + sex ratio + different mating systems (polygamy and monogamy)
#' N_t+1 = B(n=B(n=N_t * fecundity, p=juvSurv), p=sexRatio) + B(n=N_t, p=adultSurv)
#' 
#'Adult mortality + Brood mortality + offspring mortality + sex ratio (assuming that males are not a limiting factor) + different mating systems (polygamy and monogamy)
#' Nfem_t+1 = B(n=B(n=B(n=N_t * broods, p=1-breedFail) * clutch, p=juvSurv), p=sexRatio) + B(n=N_t, p=adultSurv)
#' 
#' where B(n,p) is a binomial distribution with a size parameter n a probability p.
#' 
#' @docType package
NULL


#' @importFrom methods setClass setGeneric setMethod setRefClass setOldClass
NULL