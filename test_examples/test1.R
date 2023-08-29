library(pvLRT)
signal_mat <- matrix(1, nrow(statin46), ncol(statin46))

signal_mat[1, 1] <- 10
signal_mat[2, 1] <- 5
signal_mat[3, 1] <- 2
data <- r_contin_table_zip(
  n = 1,
  row_marginals = rowSums(statin46),
  col_marginals = colSums(statin46),
  signal_mat = signal_mat,
  no_zi_idx = list(
    c(1, 1),
    c(nrow(statin46), 0),
    c(0, ncol(statin46))
  )
)[[1]]

lambda_draws <- pvBayes::pvbayes(
  data,
  "poisson",
  stan_chains = 1
)

temp <- pvBayes:::pFDR(
  lambda_draws$draws$lambda,
  test_stat = function(x){quantile(x, .05)}
)
temp$k
temp$optim
