library(tidyverse)
library(pvBayes)
library(pvLRT)

# 01: Process arguments --------------

## Pattern
## Simulation_[casename].R #1{lambda_true} #2{omega} #3{n_chains} #4{folder}
## {lambda_true} a string consisted by real numbers separated by '_', e.g. '2_3_4'
## {omega} a real number of omega in zip
## {mod} (currently ignore this argument, make a change in following part)
## {mod} (1) "All" means all models will be ran;
##       (2) model names separated by '__';
##       (3) model names start with 'except_' and separated by '__'.
## {n_chains} number of stan chains
## {folder} the name of the folder to store outputs, will be created in the path
## All arguments are separated by a space ' '



###settings
mod <- c(
  "zip_horseshoe",
  "zip_horseshoe_LKJ"
)

###command lines
cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("2", #lambda
        "0.1", #omega
        "1", #chain number
        "simulation_case0")
    } else {
      .
    }
  }


lambda_true <- cmd_args[[1]]

omega <- cmd_args[[2]] %>%
  as.numeric %>%
  {
    if (is.na(.)) {
      NULL
    } else {
      .
    }
  }


n_chains <- cmd_args[[3]] %>%
  as.numeric()

folder <- cmd_args[[4]] %>%
  tryCatch(
    .,
    error = function(e) {}
  )


###directory
is.server <- {Sys.info()["user"] == "xinweihu" & Sys.info()["sysname"] == "Linux"} %>%
  as.vector()

path_output0 <- file.path("/projects", "academic", "chakrab2", "xinweihu", "output")

if( !is.null(folder) ){

  path_output1 <- file.path(path_output0, folder)

  if( is.server ){
    if( !dir.exists(path_output1) ){
      dir.create(path_output1)
      cat("Folder created \n")
    } else{
      cat("Folder already exists \n")
    }
  } else {
    cat("Local test")
  }

} else{
  path_output1 <- path_output0
}

path <- ifelse(is.server, path_output1, "Local_test")


### system
job_id <- Sys.getenv("SLURM_JOB_ID") %>%
  as.numeric() %>%
  {
    if (is.na(.)){
      ## fake input for local test
      13225877
    } else{
      .
    }
  }

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID") %>%
  as.numeric() %>%
  {
    if (is.na(.)){
      ## fake input for local test
      20
    } else{
      .
    }
  }

seed_task <- job_id + task_id

# 02 Simulation ----------------------------------

common_sig_pos <- c(1, 10, 15, 30, 45)
#common_sig_pos <- c(30, 45)

set.seed(seed_task)

sig_pos <-
  lapply(
    1:6,
    function(x){
      sample(setdiff(c(1:46), common_sig_pos),
             size = 6-length(common_sig_pos)) %>%
        c(common_sig_pos) %>%
        sort()}
  )

zi_protect <- list()
for (j in seq_along(sig_pos)) {
  for (i in sig_pos[[j]]) {
    zi_protect <- c(zi_protect, list(c(i, j)))
  }
}

sig_pos_mat <- matrix(0, nrow = 47, ncol = 7)
for (i in seq_along(zi_protect)) {
  sig_pos_mat[ zi_protect[[i]][1], zi_protect[[i]][2] ] <- 1
}

log_lambda_signal <-
  matrix(rnorm(n = 47*7, mean = 0, sd = 0.00001), nrow = 47, ncol = 7)

for (j in 1:6){

  log_lambda_signal[sig_pos[[j]],j] <- log(as.numeric(lambda_true)) + log_lambda_signal[sig_pos[[j]],j]

}

lambda_signal <- log_lambda_signal %>% exp()

data <- r_contin_table_zip(
  n = 1,
  row_marginals = rowSums(statin46),
  col_marginals = colSums(statin46),
  signal_mat = lambda_signal,
  omega_vec = rep(omega, 7),
  no_zi_idx =
    c(
      zi_protect,
      list(
        c(nrow(statin46), 0), # the entire last row
        c(0, ncol(statin46)) # the entire last column
    )
  )
)[[1]]

### if testing on specific data, replace data here
# data <- some_source_of_data

# if( !is.null(data_input) ){
#   data <- data_input
# }

temp_pvbayes <- mod %>%
  sapply( function(m){
    cat("-----------------",m,"-----------------\n")
    tryCatch(
      pvbayes(data,
              model = m,
              stan_chains = n_chains,
              stan_seed = 1234,
              stan_iter_sampling = 1000,
              starting = "LRT"),
      error = function(e){ e }
    )
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  ) %>%
  lapply(
    function(x){ x$draws[c("lambda", "n_pred", "omega")] }
  )


temp_pvbayes_est <- temp_pvbayes %>%
  lapply(
    function(res){
      tryCatch(
        pvbayes_est(
          lambda_draws = res$lambda
        )
        ,
        error = function(e){ e }
      )

    }
  )

temp_pvlrt <- list()
temp_pvlrt$p_value <- tryCatch(
  {pvlrt(data) %>%
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

res <- list(
  job_id = job_id,
  task_id = task_id,
  seed_task = seed_task,
  omega = omega,
  data = data,
  lambda_true = {lambda_true %>%
      str_split("_") %>%
      unlist() %>%
      as.numeric()},
  lambda_true_mat = lambda_signal,
  signal_pos_list = sig_pos,
  signal_pos_mat = sig_pos_mat,
  res_pvbayes = temp_pvbayes,
  res_pvbayes_est = temp_pvbayes_est,
  res_pvlrt = temp_pvlrt,
  sample_corr_mat = cor(lambda_signal),
  fixed_num = length(common_sig_pos),
  err = err
)

out_name <- glue::glue("res_\\
  lambda{lambda_true}_\\
  job{job_id}_\\
  task{task_id}")

if (path != "Local_test") {
  saveRDS(res,
          file = glue::glue(path,"/",out_name,".RDS"))
}

