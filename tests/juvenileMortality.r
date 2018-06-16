## How to decompose the mortality in two independent but sequential events 
## accounting a fixed proportion of the total mortality?
# breedFail is a proportion of juvenile mortality correlated at the brood level
# http://en.wikipedia.org/wiki/Probability#Mathematical_treatment
#           /\
#          /  \	
#         /    \
#        /     jbr
#       /      /\
#      /      /  \
#  1-jbr  1-jind jind
#       \  /      |
#        \/       |
#       1-j       j

# survival: P(j) = P(jbr ∩ jind) = P(jbr) * P(jind)
# death:    P(-j) = P(-jbr ∪ P(-jind | jbr))
# j= jbr * jind

## d = 1 - j
## dbr = d * breedFail
## d = dbr + dind * (1-dbr) = dbr + dind - dind * dbr
# dind = (dbr - d) / (dbr - 1)
# dind = (d * breedFail - d) / (d * breedFail - 1)

######################################################
# dbr = d * breedFail                                #
# dind = (d * breedFail - d) / (d * breedFail - 1)   #
######################################################

## 1 - jbr = (1 - j) * breedFail
#  jbr = j * breedFail - breedFail + 1 = breedFail * (j-1) + 1
## j = jbr - 1 - (1 - jind) * (1 - (1 - jbr)) + 1 = jbr - (1 - jind) * jbr = jbr - (jbr - jind*jbr) 
# j = jbr * jind
# jind = -1 - ((1 - jbr) - (1 - j)) / ((1 - jbr) - 1) = -1 - (j - jbr) / jbr = -1 - (j / jbr -1)
# jind = j / jbr

######################################
# jbr = breedFail * (j-1) + 1        #
# jind = j / (breedFail * (j-1) + 1) #
######################################


j<- seq(.1, .9, by=.1) # juvenile survival
d<- 1 - j
breedFail<- .7 # Proportion of j due to breeding fail

## SCHEMA
dbr = d * breedFail                                #
dind = (d * breedFail - d) / (d * breedFail - 1)   #
data.frame(j, jCalc=(1-dbr) * (1-dind), dCalc=dbr + dind * (1-dbr), jbr, jind)

jbr = breedFail * (j-1) + 1
jind = j / (breedFail * (j-1) + 1) 
tmp<- data.frame(j, jCalc=jbr * jind, jbr, jind, dbr=1-jbr, dind=1-jind)

## Compare deterministic and probability distribution models ----

a=.6; j=.3; broods=2; b=4; fecundity=broods * b;
N0<- 5

i<- 1
j<- tmp$j[i]
jind<- tmp$jind[i]
jbr<- tmp$jbr[i]
sdistri(l<- lambda(mSurvBV.distri(broods=broods, b=b, j=jind, a=a, breedFail=1-jbr, N=N0), N0=N0))
lambda(LefkovitchPre(a=a, bj=fecundity * j))


