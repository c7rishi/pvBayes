library(tidyverse)

cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("Case1")
    } else {
      .
    }
  }

folder <- cmd_args[[1]]

function_name <- glue::glue("Simulation_{folder}.R")

N_total_rep <- 100
n_task_simul <- tryCatch(
  cmd_args[[2]],
  error = function(e){
    30
  }
  )

lambda_true <- c(1, 1.3, 1.4, 1.5, 1.7, 1.9, 2, 2.5, 3, 4)

## Pattern
## Simulation_Case1.R #1{lambda_true} #2{folder}

settings <- expand_grid(
  lambda_true = lambda_true,
  N_total_rep = N_total_rep,
  n_task_simul = n_task_simul,
  folder = folder
) %>%
  mutate(
    sbatchR_txt = glue::glue(
      "sbatchR -m 8 -t 10 \\
      -a 1-{N_total_rep}%{n_task_simul} \\
      --constrain=AVX512 \\
      {function_name} \\
      {lambda_true} \\
      {folder}"
    )
  )

settings$sbatchR_txt
settings$sbatchR_txt %>% walk(system)



