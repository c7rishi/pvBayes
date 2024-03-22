## Simulation Case 3a

library(tidyverse)
library(pvBayes)
library(pvLRT)
##source("D:/Documents/UB/Research/Code/pvBayes/simulation/Simulation_function.R")
source("/projects/academic/chakrab2/xinweihu/R_scripts/Simulation_function.R")


job_id <- Sys.getenv("SLURM_JOB_ID") %>%
  as.numeric() %>%
  replace_na(1356784)

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID") %>%
  as.numeric() %>%
  replace_na(20)

task_max <- Sys.getenv("SLURM_ARRAY_TASK_MAX") %>%
  as.numeric() %>%
  replace_na(100)


cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("1",
        "Case1",#
        "2",
        "3"
      )
    } else {
      .
    }
  }

## Simulation_Case1.R #1{lambda_true} #2{folder} #3{n_chains} #4{n_loop}

lambda_true <- cmd_args[[1]] %>% as.numeric()

folder <- cmd_args[[2]]

n_chains <- cmd_args[[3]] %>% as.numeric()

n_loop <- cmd_args[[4]] %>% as.numeric()

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

output <- tibble(
  lambda = numeric(),
  seed = integer(),
  signal_pos = list(),
  signal_mat = list(),
  data = list(),
  discovery = list(),
  metrics = list(),
  n_retry = integer(),
  stan_div = integer(),
  error = integer()
)

signal_row_common <- c(1, 15, 30, 45)

for (i in 1:n_loop) {

  n_retry <- 0



  while( n_retry <= 3 ) {

    seed_task <- (task_max * task_id * n_loop + i) * (n_retry + 1)

    set.seed(seed_task)

    signal_row_pos <-
      lapply(
        1:6, function(x){
          sample(
            setdiff(c(1:(nrow(statin46)-1)), signal_row_common),
            size = 6-length(signal_row_common)
          ) %>%
            c(signal_row_common) %>%
            sort()
        }
      )

    signal_row_pos <-
      lapply(
        1:6, function(x){
          sample(
            setdiff(c(1:(nrow(statin46)-1)), signal_row_common),
            size = 6-length(signal_row_common)
          ) %>%
            c(signal_row_common) %>%
            sort()
        }
      )

    zi_protect <- list()
    for (r in seq_along(signal_row_pos)) {
      for (j in signal_row_pos[[r]]) {
        zi_protect <- c(zi_protect, list(c(j, r)))
      }
    }

    signal_pos <- matrix(0, nrow = 47, ncol = 7)
    for (r in seq_along(zi_protect)) {
      signal_pos[ zi_protect[[r]][1], zi_protect[[r]][2] ] <- 1
    }

    log_signal_mat <-
      matrix(rnorm(n = 47*7, mean = 0, sd = 0.00001), nrow = 47, ncol = 7)

    for (j in 1:2){

      log_signal_mat[signal_row_pos[[j]],j] <-
        log(as.numeric(lambda_true)) + log_signal_mat[signal_row_pos[[j]],j]

    }
    for (j in 3:4){

      log_signal_mat[signal_row_pos[[j]],j] <-
        log(as.numeric(lambda_true*1.2)) + log_signal_mat[signal_row_pos[[j]],j]

    }

    signal_mat <- log_signal_mat %>% exp()

    data <- r_contin_table_zip(
      n = 1,
      row_marginals = rowSums(statin46),
      col_marginals = colSums(statin46),
      signal_mat = signal_mat,
      omega_vec = c(rep(0.1, (ncol(statin46)-1)), 0),
      no_zi_idx =
        c(
          zi_protect,
          list(
            c(nrow(statin46), 0), # the entire last row
            c(0, ncol(statin46)) # the entire last column
          )
        )
    )[[1]]


    temp <- Simulation(
      data = data,
      mod = c(
        "zip_horseshoe",
        "zip_horseshoe_LKJ"
      ),
      stan_chains = n_chains,
      stan_core = getOption("mc.cores", n_chains),
      stan_cov_rate = 0.1,
      stan_retry = 1
    )

    if( temp$stan_div == 1) {
      n_retry <- n_retry + 1
    } else{
      break;
    }

  }



  output <-
    add_row(
      output,
      lambda = lambda_true,
      seed = seed_task,
      signal_pos = list(signal_pos),
      signal_mat = list(signal_mat),
      data = list(data),
      discovery = temp$discovery,
      metrics = temp$metrics,
      n_retry = n_retry,
      stan_div = temp$stan_div,
      error = temp$error
    )


}


out_name <- glue::glue("{folder}_\\
  job{job_id}_\\
  task{task_id}")

if (path != "Local_test") {
  saveRDS(output,
          file = glue::glue("{path}/{out_name}.RDS"))
}

