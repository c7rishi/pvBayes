library(pvLRT)
library(tidyverse)

sig_mat <- signal_generation(
    I = 47,
    J = 7,
    sig_position = list(
      c(33,0),
      c(45,0)
    ),
    sig_strength = list(
      3,
      2
    ),
    non_sig_corr_mat = matrix(0.4, nrow = 6, ncol = 6) %>% `diag<-`(1),
    sig_corr_mat = matrix(0.6, nrow = 6, ncol = 6) %>% `diag<-`(1)
)


cor(t(sig_mat$lambda), method = "pearson")


sig_mat <- matrix(NA, nrow = 47, ncol = 7)

for(j in 1:7) {

  sig_mat[,j] <- runif(47, min = 0, max = 0.99)

}

sig_mat[33,1] <- 3
sig_mat[45,1] <- 2.5

e <- runif(47, min = 0, max = 0.05)
sig_mat[,2] <- sig_mat[,1] + e

cor(sig_mat, method = "pearson")

sig_true <- matrix(0, nrow = 47, ncol = 7 )
sig_true[33,1] <- 1
sig_true[45,1] <- 1
sig_true[33,2] <- 1
sig_true[45,2] <- 1


data <- data_generation(
    lambda_mat = sig_mat,
    signal_position = sig_true,
    omega_vec = rep(0, 7)
)


res_list <- c("poisson_indep",
              "poisson_correlated",
              "poisson_LKJ") %>%
  sapply(
    FUN = function(mod){
      data$contin_table %>%
        pvbayes(
          model = mod,
          stan_chains = 1,
          stan_iter_sampling = 1000
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )


sig_mat %>% cor(method = "pearson") %>% round(2)

cor_list <- res_list %>% lapply(
  FUN = function(res){
    res %>%
      extract_correlation_matrix(
        par = "lambda",
        method = "kendall"
      ) %>% round(2) %>% unname()
  }
)

res_list %>% map(
  function(x) {x$draws$lambda %>%
      unname() }
)

res_poisson <- data$contin_table %>%
  pvbayes(
    model = "poisson",
    stan_chains = 1,
    stan_iter_sampling = 500
  )

res_poisson %>%
  extract_correlation_matrix(
    par = "lambda",
    method = "kendall"
  ) %>% round(2) %>% unname()

res_poisson_indep <- data$contin_table %>%
  pvbayes(
    model = "poisson_indep",
    stan_chains = 1,
    stan_iter_sampling = 500
  )

res_poisson_indep %>%
  extract_correlation_matrix(
    par = "lambda",
    method = "pearson"
  ) %>% round(2) %>% unname()

res_poisson_correlated <- data$contin_table %>%
  pvbayes(
    model = "poisson_correlated",
    stan_chains = 1,
    stan_iter_sampling = 500
  )

cor_poisson_correlated <- res_poisson_correlated %>%
  extract_correlation_matrix(
    par = "lambda",
    method = "pearson"
  ) %>% round(2) %>% unname()

res_poisson_correlated_row <- data$contin_table %>%
  pvbayes(
    model = "poisson_correlated_row",
    stan_chains = 1,
    stan_iter_sampling = 500
  )

cor_poisson_correlated_row <- res_poisson_correlated_row %>%
  extract_correlation_matrix(
    par = "lambda",
    method = "pearson",
    by_row = TRUE
  ) %>% round(2) %>% unname()

cor_poisson_correlated_row[33,45]
