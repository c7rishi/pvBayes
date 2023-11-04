library(tidyverse)
library(pvLRT)
library(posterior)


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
  "poisson_indep_2",
  stan_chains = 1
)

# saveRDS(res, file = "~/test_draws.RDS", compress = "xz")
res <- readRDS("~/test_draws.RDS")

res$draws$lambda_indep
res$draws$lambda_resid
res %>%
  extract_correlation_matrix(par = "lambda",
                             method = "spearman",
                             by_row = TRUE)


