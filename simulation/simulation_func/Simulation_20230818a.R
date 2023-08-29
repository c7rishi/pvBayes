library(tidyverse)
library(pvBayes)
library(pvLRT)

# 01: Process arguments --------------

## Pattern
## Simulation.R {lambda_true} {sig_col} {sig_row} {mod}
## {lambda_true} single real value
## {sig_col} vector of integer column index separated by '_', or 'row_signal'
## {sig_row} vector of integer row index separated by '_'
## {mod} (1) "All" means all model will be ran;
##       (2) model names separated by '_';
##       (3) model names start with 'except_' and separated by '_'.
## All arguments separated by a space ' '

## All models run under true signal 1.2 on the first and the
## third AEs and the first drug.
## Simulation.R 1.2 1_3 1 All

## Horseshoe model and zip mode run under true signal 2
## on the first AEs and all drugs.
## Simulation.R 2 1 row_signal horseshoe_zip

## All models except beta prime model run under true signal 2
## on the first AEs and all drugs.
## Simulation.R 2 1 row_signal except_beta_prime

###command lines
cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("1", "row_signal", "5_45", "poisson_indep__zip_indep")
    } else {
      .
    }
  }

lambda_true <- cmd_args[[1]] %>%
  as.numeric()

sig_col <- cmd_args[[2]]

sig_row <- cmd_args[[3]]

mod <- cmd_args[[4]] %>%
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

# 02 Data generation -----------------------------

simu_data <- function(lambda_true,
                      sig_row = sig_row,
                      sig_col = sig_col){

  n_row <- nrow(statin46)
  n_col <- ncol(statin46)

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

  sim_data <- r_contin_table_zip(
    n = 1,
    row_marginals = rowSums(statin46),
    col_marginals = colSums(statin46),
    signal_mat = signal_mat,
    no_zi_idx = list(
      c(1, 1),
      c(n_row, 0),
      c(0, n_col)
    )
  )[[1]]

  signal_true <- ifelse(signal_mat <= 1, 0, 1)

  return(
    list(contin_table = sim_data, signal_true = signal_true)
  )

}

# 03 Simulation ----------------------------------

Simulation <- function(job_id = job_id,
                       task_id = task_id,
                       lambda_true = lambda_true,
                       sig_row = sig_row,
                       sig_col = sig_col,
                       mod = mod,
                       stan_chains = 4,
                       stan_parallel_chains = 4,
                       stan_seed = 123,
                       path = "/projects/academic/chakrab2/xinweihu/output/"
){

  set.seed(seed_task)

  data <- simu_data(lambda_true = lambda_true,
                    sig_row = sig_row,
                    sig_col = sig_col)

  temp_pvbayes <- tryCatch(
    {mod %>%
        lapply( function(m)pvbayes(data$contin_table,
                                   model = m,
                                   stan_chains = stan_chains,
                                   stan_seed = stan_seed,
                                   stan_parallel_chains = stan_parallel_chains)
        ) %>%
        `names<-`(mod) %>%
        lapply(function(x){
          x$draws[c("lambda", "n_pred")]
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
              test_stat = function(x){
                quantile(x, 0.95)
              }
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


  out_name <- glue::glue("res_\\
  lambda{lambda_true}_\\
  r{sig_row}_\\
  c{sig_col}_\\
  job{job_id}_\\
  task{task_id}")
  #?list.files

  assign(
    out_name,
    list(
      job_id = job_id,
      task_id = task_id,
      lambda_true = lambda_true,
      sig_row = sig_row,
      sig_col = sig_col,
      res_pvbayes = temp_pvbayes,
      res_pvbayes_est = temp_pvbayes_est,
      res_pvlrt = temp_pvlrt,
      err = err
    )
  )

  if (is.null(path)) {
    return(get(out_name))
  } else {
    saveRDS(get(out_name),
            file = glue::glue(path,out_name,".RDS"))
  }
}


# 04: run job -----------------------
res <- Simulation(
  job_id = job_id,
  task_id = task_id,
  lambda_true = lambda_true,
  sig_row = sig_row,
  sig_col = sig_col,
  mod = mod,
  stan_chains = 1,
  stan_parallel_chains = 4,
  stan_seed = 123,
  path = "/projects/academic/chakrab2/xinweihu/output/simulation_sigle_cell_signal/"
)


