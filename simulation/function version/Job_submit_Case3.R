library(tidyverse)

### pvBayes simulation setting

function_name <- "Simulation_Case3.R"

folder <- "simulation_case3"

N_total_rep <- 100
n_task_simul <- 20

model <- "except_beta_prime"
lambda_true_1 <- c(1, 1.2, 1.2, 1.5, 1.5, 1.6, 2, 2)
lambda_true_2 <- c(1, 1.5, 2, 1.5, 2, 1.6, 1.5, 2)
omega <- 0.1
n_chains <- 1

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


###
# settings <- expand_grid(
#   model = model,
#   lambda_true_1 = lambda_true_1,
#   lambda_true_2 = lambda_true_2,
#   N_total_rep = N_total_rep,
#   n_task_simul = n_task_simul,
#   n_chains = n_chains,
#   folder = folder
# ) %>%
#   mutate(
#     sbatchR_txt = glue::glue(
#       "sbatchR -m 8 -t 10 \\
#       -a 1-{N_total_rep}%{n_task_simul} \\
#       --constrain=AVX512 \\
#       {function_name} \\
#       {lambda_true_1}_{lambda_true_2} \\
#       {omega} \\
#       {model} \\
#       {n_chains} \\
#       {folder}"
#     )
#   )

settings <- tibble(
  model = model,
  lambda_true_1 = lambda_true_1,
  lambda_true_2 = lambda_true_2,
  N_total_rep = N_total_rep,
  n_task_simul = n_task_simul,
  n_chains = n_chains,
  folder = folder
) %>%
  mutate(
    sbatchR_txt = glue::glue(
      "sbatchR -m 8 -t 10 \\
      -a 1-{N_total_rep}%{n_task_simul} \\
      --constrain=AVX512 \\
      {function_name} \\
      {lambda_true_1}_{lambda_true_2} \\
      {omega} \\
      {model} \\
      {n_chains} \\
      {folder}"
    )
  )


settings$sbatchR_txt
settings$sbatchR_txt %>% walk(system)



