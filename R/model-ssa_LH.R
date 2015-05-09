# Propensity functions for each transition ----
rateFunc.LH<- compiler::cmpfun( # byte-compile the function
  function(x, params, t){
    with(params, {
      return(c(breed=b * x["Na"],
               deadA=da * (1 + x["Na"] / K) * x["Na"],
               deadJ=dj * (1  + (x["Na"] + x["Nj"]) / K) * x["Nj"],
               grow=g * x["Nj"])
      )}
    )
  }
)

# State-change matrix for each transition ----
## Notes: individuals from NxbF move to N-xb. Individuals which failed on the last reproductive event have
# higher probability to move to another habitat (c < cF) but they relax after an habitat change to evaluate the reproduction outcome in the new habitat.
transitionMat.LH<- function(params){
#                                   repro     dead    grow
  transMat<- with(params, matrix(c(      0,  -1, 0,   1,  # a
                                    clutch,   0,-1,  -1), # j
                                  ncol=4, byrow=TRUE,
                                  dimnames=list(state=c("Na", "Nj"), event=c("breed",  "deadA", "deadJ",   "grow")))
  )

  return (transMat)
}

