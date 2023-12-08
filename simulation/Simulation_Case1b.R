library(tidyverse)
library(pvBayes)
library(pvLRT)

# 01: Process arguments --------------

## Pattern
## Simulation_[casename].R #1{lambda_true} #2{omega} #3{n_chains} #4{folder}
## {lambda_true} a string consisted by real numbers seperated by '_', e.g. '2_3_4'
## {omega} a real number of omega in zip
## {mod} (currently ignore this argument, make a change in following part)
## {mod} (1) "All" means all models will be ran;
##       (2) model names separated by '__';
##       (3) model names start with 'except_' and separated by '__'.
## {n_chains} number of stan chains
## {folder} the name of the folder to store outputs, will be created in the path
## All arguments are separated by a space ' '



###settings
mod <- c("poisson_test",
         "zip_test",
         "poisson_indep_test",
         "zip_indep_test",
         "poisson_correlated_test",
         "poisson_LKJ_test")

sig_pos <- list(
  c(45,0)
)

test_stat <- function(x){quantile(x, 0.05)}

sig_thresh <- 1

###command lines
cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("2",
        "0.1",
        # "poisson__poisson_indep__poisson_indep2__zip__zip_indep",
        "1",
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

set.seed(seed_task)

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

data <-  signal_true %>%
  {
    data_generation(
      lambda_mat = .$lambda,
      signal_position = .$signal_true,
      omega_vec = c(rep(omega, 6),0)
    )
  }

### if testing on specific data, replace data here
# data <- some_source_of_data

# if( !is.null(data_input) ){
#   data <- data_input
# }

temp_pvbayes <- tryCatch(
  {mod %>%
      sapply( function(m){
        cat("-----------------",m,"-----------------\n")
        pvbayes(data$contin_table,
                model = m,
                stan_chains = n_chains,
                stan_seed = 1234,
                stan_iter_sampling = 1000,
                starting = "LRT"
        )
      },
      simplify = FALSE,
      USE.NAMES = TRUE
      ) %>%
      lapply(function(x){
        x$draws[c("lambda",
                  "n_pred",
                  "omega")]
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
            test_stat = test_stat,
            thresh = sig_thresh
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

res <- list(
  job_id = job_id,
  task_id = task_id,
  seed_task = seed_task,
  omega = omega,
  data = data$contin_table,
  lambda_true = {lambda_true %>%
      str_split("_") %>%
      unlist() %>%
      as.numeric()},
  lambda_true_mat = signal_true$lambda,
  signal_pos_list = sig_pos,
  signal_pos_mat = data$signal_true,
  signal_threshold = sig_thresh,
  res_pvbayes = temp_pvbayes,
  res_pvbayes_est = temp_pvbayes_est,
  res_pvlrt = temp_pvlrt,
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


