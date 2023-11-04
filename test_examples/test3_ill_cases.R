library(pvLRT)

data <- readRDS("C:/Users/Xinwei/Documents/test_data.RDS")

res1 <- pvbayes(
  contin_table = data$contin_table,
  model = "zip_indep",
  starting = "random",
  stan_chains = 1
)

res2 <- pvbayes(
  contin_table = data$contin_table,
  model = "zip_indep",
  starting = "LRT",
  stan_chains = 1
)

res3 <- pvbayes(
  contin_table = data$contin_table,
  model = "zip_indep",
  starting = "MAP",
  stan_chains = 1
)



res1$MCMC_convergence


