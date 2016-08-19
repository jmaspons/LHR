context("discretePopSim")

test_that("discretePopSim generics", {
  # pop<- with(scenario, discretePopSim_dispatch(broods=broods, b=b, j=jindSeason, a=a, breedFail=1 - jbrSeason,
  #                                              varJ=ifelse(pars$envVar$j, var, 0), varBreedFail=ifelse(pars$envVar$breedFail, var, 0),
  #                                              sexRatio=pars$sexRatio, matingSystem=pars$matingSystem, N0=N0, replicates=pars$replicates, tf=pars$tf))
  pop<- discretePopSim(broods=1, b=1, j=.5, a=.5, breedFail=0,
                                               varJ=0, varBreedFail=0,
                                               sexRatio=NA, matingSystem=NA, N0=2, replicates=15, tf=10)

  
  # env<- Env()
  # env<- env[env$var == 0 & env$breedFail == 0,]
  # lh<- LH()
  # lh<- lh[lh$lambda == 1,]
  # model<- Model(lh=lh, env=env)
  # model<- model[1,]
  # res<- run(model)
  # pop<- res@sim@raw[[1]][[1]]
  
  expect_is(summary(pop, dt = 1), "data.frame")
  expect_is(r(pop, dt = 1), "matrix")
  expect_is(lambda(pop, dt = 1), "matrix")
  expect_is(trendsProp(pop, dt = 1), "data.frame")
  
  expect_is(G(pop), "numeric")
  expect_is(Gmean(pop), "numeric")
})
