library(tidyverse)

### pvBayes simulation setting

function_name <- "Simulation_Case1.R"

folder <- "simulation_case1"

N_total_rep <- 100
n_task_simul <- 20

# model <- "except_beta_prime__poisson_correlated_row"
lambda_true <- c(1, 1.3, 1.5, 1.7, 1.9, 2, 2.5, 3, 3.5, 4)
omega <- "0.1"
n_chains <- 1

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


settings <- expand_grid(
  #model = model,
  lambda_true = lambda_true,
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
      {lambda_true} \\
      {omega} \\
      {n_chains} \\
      {folder}"
    )
  )

settings$sbatchR_txt
settings$sbatchR_txt %>% walk(system)



