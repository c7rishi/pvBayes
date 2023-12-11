library(tidyverse)
library(pvLRT)
library(posterior)


signal_mat <- matrix(1, nrow(statin46), ncol(statin46))

signal_mat[1, 1] <- 1.5
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

data <- temp$data

model_name()

res <-
  c(
    "zip_horseshoe",
    "zip_horseshoe_correlated",
    "zip_horseshoe_LKJ",
    "zip_horseshoe_LKJ_other"
  ) %>%
  sapply(
    function(mod){
      tryCatch(
        pvBayes::pvbayes(
          data,
          mod,
          stan_chains = 1,
          stan_iter_sampling = 1000
        ),
        error = function(e){ e }
      )

    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )


res_est <- res %>%
  lapply(
    function(r){
      r$draws$lambda %>% pvbayes_est(
        test_stat = function(x) quantile(x, 0.05)
      )
    }
  )


res_est %>%
  lapply(
    \(x){x$sig_pfdr %>% unname()}
  )





temp$res_pvbayes_est$zip_indep$sig_pfdr %>% unname()

# saveRDS(res, file = "~/test_draws.RDS", compress = "xz")
res <- readRDS("~/test_draws.RDS")

res$draws$rho_AE

cor_mat <- res %>%
  extract_correlation_matrix(
    "lambda"
  )




cor_mat %>% unname()


res$draws$lambda_indep
res$draws$lambda_resid
res %>%
  extract_correlation_matrix(par = "lambda",
                             method = "spearman",
                             by_row = TRUE)


