library(tidyverse)
library(pvBayes)
library(pvLRT)

# 01: Process arguments --------------

## Pattern
## Simulation_[teststatistics].R #1{lambda_true} #2{sig_col} #3{sig_row} #4{omega} #5{mod} #6{n_chains} #7{folder}
## {lambda_true} single real value
## {sig_col} vector of integer column index separated by '_', or 'row_signal'
## {sig_row} vector of integer row index separated by '_'
## {omega} a real number of omega in zip
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
      c("2", "row_signal", "45", "0.1", "poisson_indep__zip_indep", "1","single_cell_signal")
    } else {
      .
    }
  }

lambda_true <- cmd_args[[1]] %>%
  as.numeric()

signal_col <- cmd_args[[2]]

signal_row <- cmd_args[[3]]

omega <- cmd_args[[4]] %>%
  as.numeric %>%
  {
    if (is.na(.)) {
      NULL
    } else {
      .
    }
  }

mod <- cmd_args[[5]] %>%
  {
    if (. == "All"){

      model_name()

    } else if (str_starts(., "except_")){

      str_remove(., "except_") %>%
        {str_split(., "\\__") %>%
            unlist(.)} %>%
        setdiff(model_name(), .)

    } else {

      {str_split(., "\\__") %>%
          unlist(.)}

    }
  }

n_chains <- cmd_args[[6]] %>%
  as.numeric()

folder <- cmd_args[[7]] %>%
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


# 02 Data generation -----------------------------
simu_data <- function(lambda_true,
                      sig_row,
                      sig_col,
                      omega = NULL){

  n_row <- nrow(statin46)
  n_col <- ncol(statin46)

  row_marginals = rowSums(statin46)
  col_marginals = colSums(statin46)

  signal_mat <- matrix(1, n_row, n_col)

  sig_col_index <- sig_col %>%
    {
      if (. == "row_signal"){
        c(1 : (n_col-1))
      } else {
        str_split(., "_") %>%
          unlist() %>%
          as.numeric()
      }
    }

  sig_row_index <- sig_row %>%
    str_split("_") %>%
    unlist() %>%
    as.numeric()

  signal_mat[sig_row_index, sig_col_index] <- lambda_true

  if (is.null(omega)) {

    sim_data <- r_contin_table_zip(
      n = 1,
      row_marginals = row_marginals,
      col_marginals = col_marginals,
      signal_mat = signal_mat,
    )[[1]]
  } else {
    if (sig_col == "row_signal"){

      no_zi_idx_list <- lapply(1:length(sig_row_index), function(i) c(sig_row_index[i],0))

    } else{

      no_zi_idx_list <- mapply(c, sig_row_index, sig_col_index, SIMPLIFY = FALSE)

    }
    # browser()
    no_zi_idx_list <- no_zi_idx_list %>%
      c(.,
        list(
          c(n_row, 0),
          c(0, n_col)
        )
      )
    sim_data <- r_contin_table_zip(
      n = 1,
      row_marginals = row_marginals,
      col_marginals = col_marginals,
      signal_mat = signal_mat,
      omega_vec = c(rep(omega, n_col -1), 0),
      no_zi_idx = no_zi_idx_list
    )[[1]]
  }

  return(
    list(contin_table = sim_data, signal_true = ifelse(signal_mat <= 1, 0, 1))
  )

}


# simu_data(lambda_true = 1.5,
#           sig_row = signal_row,
#           sig_col = signal_col,
#           omega = 0.5)


# 03 Simulation ----------------------------------

Simulation <- function(job_id,
                       task_id,
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

  # browser()
  seed_task <- is.null(seed_task) %>%
    ifelse(job_id + task_id, seed_task)


  set.seed(seed_task)

  data <- simu_data(lambda_true = lambda_true,
                    sig_row = sig_row,
                    sig_col = sig_col,
                    omega = omega)

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

# 04: run job -----------------------
res <- Simulation(
  job_id = job_id,
  task_id = task_id,
  lambda_true = lambda_true,
  sig_row = signal_row,
  sig_col = signal_col,
  omega = omega,
  mod = mod,
  test_stat = function(x){quantile(x, 0.05)},
  stan_chains = n_chains,
  stan_seed = 123,
  stan_iter_sampling = 1000,
  path = ifelse(is.server, path_output1, "Local_test")
)

