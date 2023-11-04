library(tidyverse)
library(pvBayes)
library(pvLRT)

temp <- readRDS(file_names[[800]])

temp$lambda_true

temp$signal_true

temp$res_pvlrt$sig %>% unname()

temp$res_pvbayes_est$zip_indep$test_stat %>% unname
temp$res_pvbayes_est$zip_indep$sig_pfdr %>% unname()

res_lrt <- temp$data %>% pvlrt()

res_lrt %>% extract_p_value_matrix() %>% unname() %>% round(2)

########################




Simulation <- function(job_id = NULL,
                       task_id = NULL,
                       seed_task = NULL,
                       lambda_true,
                       sig_row,
                       sig_col,
                       data_input = NULL,
                       omega = NULL,
                       mod = mod,
                       test_stat,
                       stan_chains,
                       stan_seed,
                       stan_iter_sampling,
                       path
){

  browser()
  seed_task <- is.null(seed_task) %>%
    ifelse(job_id + task_id, seed_task)


  set.seed(seed_task)

  signal_generation(
    I = 47,
    J = 7,
    sig_row = sig_row,
    sig_col = sig_col,
    #non_sig_prob = 0.95,
    sig_mean_vec = rep(lambda_true, 6),
    sig_corr_mat = matrix(0.2, nrow = 6, ncol = 6) %>% `diag<-`(1)
  ) %>% {.[[1]]} %>% round(2)

  data <- signal_generation(
    I = 47,
    J = 7,
    sig_row = sig_row,
    sig_col = sig_col,
    #non_sig_prob = 0.95,
    sig_mean_vec = rep(lambda_true, 6),
    sig_corr_mat = matrix(0.2, nrow = 6, ncol = 6) %>% `diag<-`(1)
  ) %>%
    {
      data_generation(
        lambda_mat = .$lambda,
        signal_position = .$signal_true,
        omega_vec = c(rep(omega, 6),0)
      )
    }

  if( !is.null(data_input) ){
    data <- data_input
  }

  temp_pvbayes <- tryCatch(
    {mod %>%
        lapply( function(m)pvbayes(data$contin_table,
                                   model = m,
                                   stan_chains = stan_chains,
                                   stan_seed = stan_seed,
                                   stan_iter_sampling = stan_iter_sampling,
                                   starting = "LRT")
        ) %>%
        `names<-`(mod) %>%
        lapply(function(x){
          x$draws[c("lambda", "n_pred", "omega")]
        })
    },
    error = function(e){ e }
  )

  temp_pvbayes_est <- tryCatch(
    {temp_pvbayes %>%
        lapply(
          function(res){
            pvbayes_est(
              lambda_draws = res$lambda,
              test_stat = test_stat
            )
          }
        )
    },
    error = function(e){ e }
  )

  temp_pvlrt <- list()
  temp_pvlrt$p_value <- tryCatch(
    {pvlrt(data$contin_table) %>%
        extract_p_value_matrix()},
    error = function(e){ e }
  )

  err <- ifelse(inherits(temp_pvbayes, "error")|
                  inherits(temp_pvbayes_est, "error")|
                  inherits(temp_pvlrt$p_value, "error"),
                1, 0)

  temp_pvlrt$sig <- temp_pvlrt$p_value %>%
    {ifelse(.<0.05, 1, 0)}

  temp_pvlrt$sig[,7] <- 0

  out_name <- glue::glue("res_\\
  lambda{lambda_true}_\\
  r{sig_row}_\\
  c{sig_col}_\\
  job{job_id}_\\
  task{task_id}")

  assign(
    out_name,
    list(
      job_id = job_id,
      task_id = task_id,
      seed_task = seed_task,
      lambda_true = lambda_true,
      sig_row = sig_row,
      sig_col = sig_col,
      omega = omega,
      data = data$contin_table,
      signal_true = data$signal_true,
      res_pvbayes = temp_pvbayes,
      res_pvbayes_est = temp_pvbayes_est,
      res_pvlrt = temp_pvlrt,
      err = err
    )
  )

  if (path == "Local_test") {
    return(get(out_name))
  } else {
    saveRDS(get(out_name),
            file = glue::glue(path,"/",out_name,".RDS"))
  }
}


res <- Simulation(
  job_id = 123,
  task_id = 456,
  # data_input = list(
  #   contin_table = temp$data,
  #   sigal_true = temp$signal_true
  # ),
  # seed_task = temp$seed_task,
  lambda_true = temp$lambda_true,
  sig_row = 45,
  sig_col = "row_signal",
  omega = 0.1,
  mod = c("zip_indep"),
  test_stat = function(x){quantile(x, 0.05)},
  stan_chains = 1,
  stan_seed = 123,
  stan_iter_sampling = 1000,
  path = ifelse(FALSE, path_output1, "Local_test")
)

res$res_pvlrt$sig %>% unname()
res$res_pvbayes_est$zip_indep$sig_pfdr %>% unname() %>% round(3)
res$res_pvbayes_est$zip_indep$test_stat %>% unname() %>% round(3)
res$res_pvbayes_est$zip_indep$k
