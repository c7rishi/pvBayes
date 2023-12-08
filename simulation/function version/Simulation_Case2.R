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

# mod <- cmd_args[[3]] %>%
#   {
#     if (. == "All"){
#
#       model_name()
#
#     } else if (str_starts(., "except_")){
#
#       str_remove(., "except_") %>%
#         {str_split(., "\\__") %>%
#             unlist(.)} %>%
#         setdiff(model_name(), .)
#
#     } else {
#
#       {str_split(., "\\__") %>%
#           unlist(.)}
#
#     }
#   }

# mod <- c("poisson",
#          "zip",
#          "poisson_indep",
#          "zip_indep")

mod <- c("poisson",
         "zip",
         "poisson_indep",
         "poisson_indep2",
         "zip_indep",
         "poisson_correlated",
         "poisson_LKJ",
         "poisson_ridge")

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

# 02 Simulation ----------------------------------

Simulation <- function(job_id, #job id extracted from CCR
                       task_id, #task id extracted from CCR
                       seed_task = NULL, #specified seed (not use)
                       lambda_true, #string for true lambda value
                       sig_pos, #list of where signals exist
                       data_input = NULL, #specified data (not use)
                       omega = 0,
                       mod, #vector of model names
                       test_stat, #function of test statistics
                       stan_chains, #number of stan chains
                       stan_seed, #stan initial seeds
                       stan_iter_sampling, #number of stan sample
                       path #path of rds storage
){

  # browser()
  seed_task <- is.null(seed_task) %>%
    ifelse(job_id + task_id, seed_task)

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

  if( !is.null(data_input) ){
    data <- data_input
  }

  temp_pvbayes <- tryCatch(
    {mod %>%
        sapply( function(m){
          cat("-----------------",m,"-----------------\n")
          pvbayes(data$contin_table,
                  model = m,
                  stan_chains = stan_chains,
                  stan_seed = stan_seed,
                  stan_iter_sampling = stan_iter_sampling,
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
  job{job_id}_\\
  task{task_id}")

  assign(
    out_name,
    list(
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

# 03: run job -----------------------
res <- Simulation(
  job_id = job_id,
  task_id = task_id,
  lambda_true = lambda_true,
  sig_pos = list(
    c(5,0),
    c(45,0)
  ),
  omega = omega,
  mod = mod,
  test_stat = function(x){quantile(x, 0.05)},
  stan_chains = n_chains,
  stan_seed = 123,
  stan_iter_sampling = 1000,
  path = ifelse(is.server, path_output1, "Local_test")
)

