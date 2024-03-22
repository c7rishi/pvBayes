library(tidyverse)

cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("Case1", ##job/folder
        "2", ##number of chains and parallel
        "1000", ##number of total replications
        "10", ##loops in each task,
        "25" ##number of simultaneously job run
        )
    } else {
      .
    }
  }

folder <- cmd_args[[1]]

function_name <- glue::glue("Simulation_{folder}.R")

n_chains <- cmd_args[[2]] %>% as.numeric()

N_total <- cmd_args[[3]] %>% as.numeric()

n_per_task <- cmd_args[[4]] %>% as.numeric()

N_duplicate <- N_total / n_per_task

n_task_simul <- cmd_args[[5]] %>% as.numeric()

lambda_true <- c(1.1, 1.3, 1.4, 1.5, 1.7, 1.9, 2, 2.5, 3, 4)

## Pattern
## Simulation_Case1.R #1{lambda_true} #2{folder} #3{n_chains} #4{n_per_task}
settings <- expand_grid(
  lambda_true = lambda_true
) %>%
  mutate(
    sbatchR_txt = glue::glue(
      "sbatchR -m 8 -n {n_chains} -t 10 \\
      -a 1-{N_duplicate}%{n_task_simul} \\
      --constrain=AVX512 \\
      {function_name} \\
      {lambda_true} \\
      {folder} \\
      {n_chains} \\
      {n_per_task} \\
      "
    )
  )

settings$sbatchR_txt
settings$sbatchR_txt %>% walk(system)



