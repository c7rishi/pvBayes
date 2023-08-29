library(tidyverse)

### pvBayes simulation setting

function_name <- "Simulation.R"

folder <- "single_cell_signal"

N_total_rep <- 200
n_task_simul <- 20

model <- "except_beta_prime"
lambda_true <- c(1.3, 1.5, 1.7, 2.5, 3.5)
sig_col <- "1"
sig_row <- "45"
n_chains <- 2


## Pattern
## Simulation.R {lambda_true} {sig_col} {sig_row} {mod} {n_chains} {folder}
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

###

settings <- expand_grid(
  model = model,
  lambda_true = lambda_true,
  sig_col = sig_col,
  sig_row = sig_row,
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
      {sig_col} \\
      {sig_row} \\
      {model} \\
      {n_chains} \\
      {folder}"
    )
  )

settings$sbatchR_txt
settings$sbatchR_txt %>% walk(system)



