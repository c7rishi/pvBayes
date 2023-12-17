library(tidyverse)
library(posterior)
library(pvLRT)

folder <- "simulation_case1"
file_path0 <- "D:/Documents/UB/Research/Code/CCR/Download"
file_path1 <- file.path(file_path0, folder)

file_names <- list.files(path = file_path1,
                         pattern = "\\.RDS",
                         full.names = TRUE)

simu_df <- readRDS(paste(file_path1,"_full.RDS", sep = ""))

simu_df %>%
  select(index, model, lambda, Sensitivity) %>%
  pivot_wider(
    names_from = model,
    values_from = Sensitivity
  ) %>%
  select(
    index, lambda, LRT, zip_horseshoe_LKJ
  ) %>%
  filter(
    lambda == 3
  ) %>%
  filter(
    LRT < zip_horseshoe_LKJ
  )







temp <- try(readRDS(file_names[[817]]))
temp$res_pvbayes_est$zip_horseshoe_LKJ$sig_pfdr %>% unname()
temp$res_pvlrt$sig %>% unname()




log_lambda_signal <-
  matrix(rnorm(n = 47*7, mean = 0, sd = 0.00001), nrow = 47, ncol = 7)

rho0 <- 0.4
sigma_signal0 <- 0.005
x <- rnorm(47, mean = 0, sd = sigma_signal0)
x[c(1:5)] <- log(2) + x[c(1:5)]

sigma_signal <- sd(x)

sigma_noise <- sigma_signal * sqrt( (1 - rho0^2)/rho0^2 ); sigma_noise

# x <- rnorm(47, mean = 0, sd = sigma_signal)
# x[c(1:5)] <- log(2) + x[c(1:5)]

y <- x + rnorm(47, 0, sigma_noise)

cor(x,y)

log_lambda_signal <- rnorm(47*7, mean = 0, sd = sigma_signal0) %>%
  matrix(nrow = 47)

common_sig_pos <- c(1, 10, 15, 30, 45)
common_sig_pos <- c(1, 10)
# sig_pos0 <- sig_pos

sig_pos <- lapply(1:6, function(x){
  sample(setdiff(c(1:46), common_sig_pos),
         size = 6-length(common_sig_pos)) %>%
    c(common_sig_pos) %>%
    sort()})

log_lambda_signal <-
  matrix(rnorm(n = 47*7, mean = 0, sd = 0.00001), nrow = 47, ncol = 7)

for (j in 1:6){

  log_lambda_signal[sig_pos[[j]],j] <- log(1.5) + log_lambda_signal[sig_pos[[j]],j]

}

log_lambda_signal %>% round(3)

cor(log_lambda_signal) %>% round(3)

# sqrt(sigma_signal^2/ (sigma_signal^2 + sigma_noise^2))

#
# log_lambda_signal[,1] <- rnorm(47, 0, sigma_signal0)
#
#
# log_lambda_signal[c(1,15,30,45),1] <-  log(2) + log_lambda_signal[c(1,15,30,45),1]
#
# sigma_signal <- sd(log_lambda_signal[,1])
#
# sigma_noise <- sigma_signal * sqrt( (1 - rho0^2)/rho0^2 ); sigma_noise
#
# for(j in 2:6) {
#
#   log_lambda_signal[, j] <- log_lambda_signal[, 1]+ rnorm(47, 0, sigma_noise)
#
# }
#

# sigma_signal <- 0.01
# sigma_noise <- 0.005
# log_lambda_col1 <- rnorm(47, 0, sigma_signal)
# log_lambda_col1[c(1,15,30,45)] <- rnorm(4, log(3), sigma_signal)
# log_lambda_col2 <- log_lambda_col1 + rnorm(47, 0, sigma_noise)
# log_lambda_col3_6 <- rnorm(47*4, 0, sigma_noise) %>% matrix(nrow = 47)
#
#
#
# log_lambda_signal <- cbind(
#   log_lambda_col1,
#   log_lambda_col2,
#   log_lambda_col3_6,
#   rep(0, 47)
# )

lambda_signal <- log_lambda_signal %>% exp()# %>% round(3)

log_lambda_signal %>% cor()

lambda_signal %>% round(3)

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


res_lrt <- data %>% pvlrt()
res_lrt %>% extract_p_value_matrix() %>% unname() %>%
  {(.<0.05)*1 } %>%
  apply( 2, function(x){which(x==1)})

sig_pos

res <- c(
  "zip_horseshoe_LKJ"#,
  #"zip_horseshoe_LKJ_new2",
  #"zip_horseshoe_LKJ_new3"
) %>% sapply(

  function(mod){
    cat(glue::glue("----------{mod}-----------\n"))
    pvBayes::pvbayes(
      contin_table = data,
      model = mod,
      stan_chains = 1,
      stan_seed = 1234,
      stan_iter_sampling = 1000,
      starting = "LRT"
    )
  },
  USE.NAMES = TRUE,
  simplify = FALSE


)
lambda_signal[28,4] %>% quantile2(0.05)
res$zip_horseshoe_LKJ$draws$lambda[28,] %>% unname() %>% quantile2(0.05)

res_est <- res %>%
  lapply(
    function(r){
      tryCatch(
        pvbayes_est(
          lambda_draws = r$draws$lambda,
          test_stat = function(x)quantile(x, 0.05),
          thresh = 1.1
        )
        ,
        error = function(e){ e }
      )

    }
  )

# res_est$zip_horseshoe_LKJ %>%
#   lapply(
#     function(est) {est$sig_pfdr %>% unname()}
#   )

dis_lrt <- res_lrt %>% extract_p_value_matrix() %>% unname() %>%
  {(.<0.05)*1 } %>%
  apply( 2, function(x){which(x==1)})

dis_bayes <- res_est$zip_horseshoe_LKJ$sig_pfdr %>% unname() %>%
  apply( 2, function(x){which(x==1)})

sig_pos

map2(
  sig_pos, dis_lrt[-7], setdiff
)

map2(
  sig_pos, dis_bayes[-7], setdiff
)

map2(
  dis_lrt[-7], sig_pos,  setdiff
)

map2(
  dis_bayes[-7], sig_pos,  setdiff
)

temp$lambda_true_mat %>%
  cor()



res$zip_horseshoe_LKJ$draws$lambda
res$zip_horseshoe_LKJ %>%
  extract_correlation_matrix(
    "lambda",
    use = "complete.obs"
  ) %>% unname()
