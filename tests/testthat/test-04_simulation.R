test_that("simulation for DM, DTM, DMLM", {
   data_dm <- simulate_DM(n_obs = 100,
                          n_vars = 30,
                          n_taxa = 50,
                          n_relevant_vars = 5,
                          n_relevant_taxa = 10,
                          beta_min = 1,
                          beta_max = 1.5,
                          signoise = 0.5,
                          n_reads_min = 5000,
                          n_reads_max = 10000,
                          theta0 = 0.01,
                          rho = 0.5)

   str(data_dm)
   # View(data_dm$X)
   # View(data_dm$Y)
   # View(data_dm$betas)
})




test_that("simulation for Xsim and DTM", {
  Xmat <- Xsim(subject_sim = 100,
               tree = NULL,
               num_leaf = 10,
               covariates_sim = 50,
               rho = 0.5,
               Sigma = NULL,
               num_branch = 3,
               num_cov = 5,
               seed = 555)
  str(Xmat)
  View(Xmat$Sigma)

  data_dtm <- simulate_DTM(subject_sim = 100,
                           tree = Xmat$tree,
                           num_leaf = 10,
                           covariates_sim = 50,
                           rho = 0.5,
                           X = Xmat$X,
                           zeta_sim = Xmat$zeta,
                           rep = 5,
                           num_branch = 3,
                           num_cov = 5,
                           phi_min = 0.9,
                           phi_max = 1.2,
                           seed = 555)

  str(data_dtm$Y)
  str(data_dtm)
  # View(data_dtm$X)
  View(data_dtm$Y)
  View(data_dtm$Y[[1]])
  View(data_dtm$X)
  # View(data_dtm$phi_sim)
  # View(data_dtm$zeta_sim)
  # View(data_dtm$phi_sim %*% t(data_dtm$zeta_sim))
  phytools::plotTree(Xmat$tree, node.numbers = T)
  phytools::plotTree(data_dtm$tree, node.numbers = T)
})
