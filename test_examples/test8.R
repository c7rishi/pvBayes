simu_df %>%
  filter(lambda>2)

temp <- try(readRDS(file_names[[501]]))

res <- temp$data %>%
  pvBayes::pvbayes(
    model = "zip_horseshoe_correlated",
    stan_chains = 1,
    stan_seed = 1234,
    stan_iter_sampling = 1000,
    starting = "LRT"
  )

res %>%
  pvBayes::extract_correlation_matrix("lambda") %>%
  unname() %>%
  round(3)


signal_true <- signal_generation(
  I = 47,
  J = 7,
  sig_position = sig_pos,
  sig_strength = {lambda_true %>%
      str_split("_") %>%
      unlist() %>%
      as.numeric() %>%
      as.list()},
  non_sig_corr_mat = matrix(0.4, nrow = 6, ncol = 6) %>% `diag<-`(1),
  sig_corr_mat = matrix(0.6, nrow = 6, ncol = 6) %>% `diag<-`(1)
)



sigma_signal <- 0.01
sigma_noise <- 0.005
log_lambda_col1 <- rnorm(47, 0, sigma_signal)
log_lambda_col1[c(1,15,30,45)] <- rnorm(4, log(3), sigma_signal)
log_lambda_col2 <- log_lambda_col1 + rnorm(47, 0, sigma_noise)
log_lambda_col3_6 <- rnorm(47*4, 0, sigma_noise) %>% matrix(nrow = 47)

log_lambda_signal <- cbind(
  log_lambda_col1,
  log_lambda_col2,
  log_lambda_col3_6,
  rep(0, 47)
)

lambda_signal <- log_lambda_signal %>% exp() %>% round(3)

data <- r_contin_table_zip(
  n = 1,
  row_marginals = rowSums(statin46),
  col_marginals = colSums(statin46),
  signal_mat = lambda_signal,
  # omega_vec = omega_tru,
  no_zi_idx = list(
    c(1, 1),
    c(nrow(statin46), 0), # the entire last row
    c(0, ncol(statin46)) # the entire last column
  )
)[[1]]

res <- data %>%
  pvBayes::pvbayes(
    model = "zip_horseshoe_LKJ_new3",
    stan_chains = 1,
    stan_seed = 1234,
    stan_iter_sampling = 1000,
    starting = "LRT"
  )

res$draws$log_lambda <- res$draws$lambda %>% log()

res_est <- res$draws$lambda %>%
  pvbayes_est(
    test_stat = function(x)quantile(x, 0.05),
    thresh = 1.05
  )

log_lambda_signal %>% unname() %>% cor(method = "pearson") %>% round(5)


lambda_cor_mat <- array(NA, dim = c(1000,7,7) )

for( i in 1:1000) {

  lambda_cor_mat[i,,] <- res$draws$lambda %>%
    posterior::draws_of() %>%
    .[i,,] %>%
    # na_if(0) %>%
    log() %>%
    {ifelse(is.infinite(.), NA, .)} %>%
    cor(use = "complete.obs", method = "pearson")

}

lambda_cor_mat %>% posterior::rvar()


res$draws$lambda[c(1,15,30,45),] %>%
  unname() %>%
  round(3)

res_est$k

res_est$sig_pfdr %>% unname() %>% .[c(1,15,30,45),]

res_lrt <- data %>% pvlrt()
res_lrt %>% extract_p_value_matrix() %>% unname() %>%
  .[c(1,15,30,45),]




