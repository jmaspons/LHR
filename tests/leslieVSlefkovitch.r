a<- .7
s<- .6
j<- .5
bj<- 2
fecundity<- 4


## Compare leslie VS lefkovitch matrices ----
LesliePre<- function(a, s, bj, AFR=1){
  rows<- max(AFR, 2)
  mat<- matrix(0, nrow=rows, ncol=rows)
  mat[1, rows]<- bj
  mat[rows, rows]<- a
  stages<- rep("A", 2)
  if (AFR > 1){
    stages<- c(rep("S", rows-1), "A")
    for (i in 2:(rows)){
      mat[i,i-1]<- s
    }
  }else{
    mat[,1]<- mat[,2]
  }
  dimnames(mat)<- list(stages, stages)
  class(mat)<- "leslieMatrix"
  return (mat)
}

LefkovitchPre<- function(a, s, bj, AFR=1){
  if (AFR < 3) return(LesliePre(a, s, bj, AFR))
  
  mat<- matrix(0, nrow=2, ncol=2)
  mat[1, 2]<- bj
  mat[2, 2]<- a
  stages<- NA_character_
  if (AFR > 1){
    stages<- c("S", "A")
    mat[1, 1]<- s - s^(AFR) # survival - transition to adults
    mat[2, 1]<- s^(AFR)
    if (sum(mat[, 1]) != s) warning("The probability of subadult transitions doesn't match survival + growth to adult.")
  }else{
    stages<- c("A", "A")
    mat[, 1]<- mat[, 2]
  }
  dimnames(mat)<- list(stages, stages)
  class(mat)<- "leslieMatrix"
  return (mat)
}

m<- list()
m[[1]]<- LesliePre(a, s, bj, AFR=1)
m[[2]]<- LefkovitchPre(a, s, bj, AFR=1)

m[[3]]<- LesliePre(a, s, bj, AFR=2)
m[[4]]<- LefkovitchPre(a, s, bj, AFR=2)

m[[5]]<- LesliePre(a, s, bj, AFR=3)
m[[6]]<- LefkovitchPre(a, s, bj, AFR=3)

m[[7]]<- LesliePre(a, s, bj, AFR=4)
m[[8]]<- LefkovitchPre(a, s, bj, AFR=4)

sapply(m, eigen.analisys2df)
sapply(m, lambda)

m[5:8]

## CONCLUSION: there are differences in the deterministic models.
# see Fujiwara, M., & Diaz-Lopez, J. (2017). Constructing stage-structured matrix population models from life tables: comparison of methods. PeerJ, 5, e3971. https://doi.org/10.7717/peerj.3971

## Simulations -----
mFit.t.leslie<- function(fecundity, j, s, a, AFR, N0, replicates, tf, maxN=100000){
  if (AFR == 1){
    stages<- "A"
  } else stages<- c(paste0("S", 1:(AFR-1)), "A")
  
  pop<- array(0, dim=c(replicates, tf+1, AFR), dimnames=list(replicate=NULL, t=0:tf, age=stages))
  pop[, 1, AFR]<- N0
  for (t in 1:tf){
    # Juvenile survivors
    pop[, t+1, 1]<- rbinom(replicates, pop[, t, AFR] * fecundity, j)
    # Add adult survivors
    pop[, t+1, AFR]<- pop[, t+1, AFR] + rbinom(replicates, pop[, t, AFR], a)
    # Subadults (growth transitions)
    if (AFR > 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR] + apply(pop[, t, -AFR], 2, function(x) rbinom(replicates, x, s))
    }else if (AFR == 2){
      pop[, t+1, 2:AFR]<- pop[, t+1, 2:AFR, drop=TRUE] + apply(pop[, t, -AFR, drop=FALSE], 2, function(x) rbinom(replicates, x, s))
    }
    
    pop[which(pop[,t+1, AFR] > maxN), t+1, AFR]<- maxN
  }
  
  pop<- pop[order(pop[, ncol(pop), AFR]),, AFR] # Keep adults only and drop subadults
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

pop<- mFit.t.leslie(fecundity, j, s, a, AFR=3, N0=5, replicates=10000, tf=10, maxN=100000)
pop<- mFit.t.leslie(fecundity, j, s, a, AFR=2, N0=5, replicates=10000, tf=10, maxN=100000)
pop<- mFit.t.leslie(fecundity, j, s, a, AFR=1, N0=5, replicates=10000, tf=10, maxN=100000)
hist(pop)
  
mFit.t.lefkovitch<- function(fecundity, j, s, a, AFR, N0, replicates, tf, maxN=100000){
  if (AFR == 1){
    stages<- "A"
    pop<- array(0, dim=c(replicates, tf+1, 1), dimnames=list(replicate=NULL, t=0:tf, age=stages))
    pop[, 1, 1]<- N0
  } else {
    stages<- c("S", "A")
    pop<- array(0, dim=c(replicates, tf+1, 2), dimnames=list(replicate=NULL, t=0:tf, age=stages))
    pop[, 1, 2]<- N0
  }
  
  for (t in 1:tf){
    if (AFR > 1){
      # Recruitment
      pop[, t+1, 1]<- rbinom(replicates, pop[, t, 2] * fecundity, j)
      # Add adult survivors
      pop[, t+1, 2]<- pop[, t+1, 2] + rbinom(replicates, pop[, t, 2], a)
      # Subadults
      pop[, t+1, 2]<- pop[, t+1, 2] + rbinom(replicates, pop[, t, 1], s^AFR) # S -> A
      pop[, t+1, 1]<- pop[, t+1, 1] + rbinom(replicates, pop[, t, 1], s - s^AFR) # S -> S
      
      pop[which(pop[,t+1, 2] > maxN), t+1, 2]<- maxN
    } else {
      # Recruitment
      pop[, t+1, 1]<- rbinom(replicates, pop[, t, 1] * fecundity, j)
      # Add adult survivors
      pop[, t+1, 1]<- pop[, t+1, 1] + rbinom(replicates, pop[, t, 1], a)
      
      pop[which(pop[,t+1, 1] > maxN), t+1, 1]<- maxN
    }
  }
  if (AFR > 1){
    pop<- pop[order(pop[, ncol(pop), 2]),, 2] # Keep adults only and drop subadults
  } else {
    pop<- pop[order(pop[, ncol(pop), 1]),, 1] # Keep adults only and drop subadults
  }
  pop<- extinctNA(pop)
  class(pop)<- c("discretePopSim", "matrix")
  
  return(pop)
}

pop<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=3, N0=5, replicates=10000, tf=10, maxN=100000)
pop<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=2, N0=5, replicates=10000, tf=10, maxN=100000)
pop<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=1, N0=5, replicates=10000, tf=10, maxN=100000)
hist(pop)


## Compare leslie VS lefkovitch simulations ----
replicates<- 99999
tf<- 5

les<- list()
les[[1]]<- mFit.t.leslie(fecundity, j, s, a, AFR=1, N0=5, replicates=replicates, tf=tf, maxN=100000)
les[[2]]<- mFit.t.leslie(fecundity, j, s, a, AFR=2, N0=5, replicates=replicates, tf=tf, maxN=100000)
les[[3]]<- mFit.t.leslie(fecundity, j, s, a, AFR=3, N0=5, replicates=replicates, tf=tf, maxN=100000)
les[[4]]<- mFit.t.leslie(fecundity, j, s, a, AFR=10, N0=5, replicates=replicates, tf=tf, maxN=100000)

lef<- list()
lef[[1]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=1, N0=5, replicates=replicates, tf=tf, maxN=100000)
lef[[2]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=2, N0=5, replicates=replicates, tf=tf, maxN=100000)
lef[[3]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=3, N0=5, replicates=replicates, tf=tf, maxN=100000)
lef[[4]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=10, N0=5, replicates=replicates, tf=tf, maxN=100000)

sapply(les, summary)
sapply(lef, summary)

sapply(les, summary) - sapply(lef, summary)

par(mfrow=c(2, 2))
kk<- lapply(les, hist)
kk<- lapply(lef, hist)


## Profiling ----
profvis::profvis({
  les<- list()
  les[[1]]<- mFit.t.leslie(fecundity, j, s, a, AFR=1, N0=5, replicates=replicates, tf=tf, maxN=100000)
  les[[2]]<- mFit.t.leslie(fecundity, j, s, a, AFR=2, N0=5, replicates=replicates, tf=tf, maxN=100000)
  les[[3]]<- mFit.t.leslie(fecundity, j, s, a, AFR=3, N0=5, replicates=replicates, tf=tf, maxN=100000)
  les[[4]]<- mFit.t.leslie(fecundity, j, s, a, AFR=10, N0=5, replicates=replicates, tf=tf, maxN=100000)
})

profvis::profvis({
  lef<- list()
  lef[[1]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=1, N0=5, replicates=replicates, tf=tf, maxN=100000)
  lef[[2]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=2, N0=5, replicates=replicates, tf=tf, maxN=100000)
  lef[[3]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=3, N0=5, replicates=replicates, tf=tf, maxN=100000)
  lef[[4]]<- mFit.t.lefkovitch(fecundity, j, s, a, AFR=10, N0=5, replicates=replicates, tf=tf, maxN=100000)
})

profvis::profvis({
  ori<- list()
  ori[[1]]<- mFit.t(fecundity, j, a, N0=5, replicates=replicates, tf=tf, maxN=100000)
  ori[[2]]<- mFit.t(fecundity, j, a, N0=5, replicates=replicates, tf=tf, maxN=100000)
  ori[[3]]<- mFit.t(fecundity, j, a, N0=5, replicates=replicates, tf=tf, maxN=100000)
  ori[[4]]<- mFit.t(fecundity, j, a, N0=5, replicates=replicates, tf=tf, maxN=100000)
})

