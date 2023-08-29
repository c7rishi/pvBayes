library(tidyverse)
library(pvLRT)

signal_mat <- matrix(1, nrow(statin46), ncol(statin46))

signal_mat[1, 1] <- 10
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

res <- pvBayes::pvbayes(
  data,
  "poisson",
  stan_chains = 1
)


data
res$E %>% round(2)
res$draws$n_pred  %>% mean() %>% round(2)
