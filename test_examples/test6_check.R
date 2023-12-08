library(pvLRT)
library(tidyverse)

# signal_mat <- matrix(1, nrow(statin46), ncol(statin46))
#
# signal_mat[30, 1:2] <- 3
# signal_mat[45, 1:2] <- 2.5
#
# for (ii in seq_len(length(signal_mat)) ){
#   signal_mat[ii] <- signal_mat[ii] + runif(n = 1, min = 0, max = 0.05)
#
# }
#
# set.seed(1234)
# data <- r_contin_table_zip(
#   n = 1,
#   row_marginals = rowSums(statin46),
#   col_marginals = colSums(statin46),
#   signal_mat = signal_mat,
#   no_zi_idx = list(
#     c(30, 1),
#     c(30, 2),
#     c(nrow(statin46), 0),
#     c(0, ncol(statin46))
#   )
# )[[1]]
set.seed(1234)
sig_mat <- signal_generation(
  I = 47,
  J = 7,
  sig_position = list(
    c(1,0)
  ),
  sig_strength = list(
    2
  ),
  non_sig_corr_mat = matrix(0.4, nrow = 6, ncol = 6) %>% `diag<-`(1),
  #sig_mean_vec = rep(lambda_true, 6),
  sig_corr_mat = matrix(0.6, nrow = 6, ncol = 6) %>% `diag<-`(1)
)


data <- data_generation(
  row_marginals = rowSums(pvLRT::statin46),
  col_marginals = colSums(pvLRT::statin46),
  lambda_mat = sig_mat$lambda,
  signal_position = sig_mat$signal_true,
  omega_vec = rep(0, 7)
)



res_list <- c("poisson",
              "zip",
              "poisson_indep",
              "poisson_indep2",
              "zip_indep",
              "poisson_correlated",
              "poisson_LKJ",
              "poisson_ridge") %>%
  sapply(
    FUN = function(mod){
      cat("-----------------",mod,"-----------------\n")
      temp$data %>%
        pvbayes(
          model = mod,
          stan_chains = 1,
          stan_iter_sampling = 1000
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )


lambda_list <- res_list %>%
  lapply(
    FUN = function(x) x$draws$lambda %>% unname()
  )

est_list <- res_list %>%
  lapply(
    FUN = function(x) x$draws$lambda
  ) %>%
  lapply(
    FUN = function(x) {
      x %>% pvbayes_est(
        test_stat = function(x)quantile(x, 0.05),
        alpha = .05
      )
    }

  )


est_list %>% lapply(
  function(x) {x$sig_pfdr %>% unname()}
)


res_lrt <- data$contin_table %>% pvlrt()
res_lrt %>% extract_p_value_matrix() %>%
  unname() %>%
  {ifelse(.<0.05, 1, 0)}

