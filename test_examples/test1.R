library(pvLRT)
library(tidyverse)
library(pvBayes)
signal_mat <- matrix(1, nrow(statin46), ncol(statin46))

signal_mat[30, 1] <- 10
signal_mat[30, 2] <- 5
# signal_mat[10, 1] <-

set.seed(1234)
data <- r_contin_table_zip(
  n = 1,
  row_marginals = rowSums(statin46),
  col_marginals = colSums(statin46),
  signal_mat = signal_mat,
  no_zi_idx = list(
    c(45, 1),
    c(45, 5),
    c(nrow(statin46), 0),
    c(0, ncol(statin46))
  )
)[[1]]

saveRDS(data, file = "~/test_contin_data.RDS", compress = "xz")
data <- readRDS("~/test_contin_data.RDS")


res_lrt1 <- pvlrt(data)
res_lrt1 %>%
  extract_p_value_matrix() %>%
  {ifelse(.<.05, 1, 0)} %>%
  unname()


lambda_draws <- pvBayes::pvbayes(
  data,
  "zip",
  stan_chains = 1
)

saveRDS(lambda_draws, file = "~/test_lambda_draws.RDS", compress = "xz")
# lambda_draws <- readRDS("~/test_lambda_draws.RDS")


temp <- pvBayes:::pFDR(
  lambda_draws$draws$lambda,
  test_stat = function(x){mean(x>1)}
)
temp$test_stat %>% View()
temp$sig_pfdr %>% unname()
temp$k
temp$BayesTIE %>% View()
temp$pFDR
lambda_draws$draws$lambda[47,1] %>% {mean(.>1)}
lambda_draws$draws$lambda[47,1] %>%
  posterior::as_draws_list() %>%
  unlist() %>%
  unname() %>%
  summary()


model_name()
pvBayes_setup()
