library(tidyverse)
library(pvLRT)
# library(pvBayes)
library(ggplot2)


# tool functions

signal.table <- function(signal_detect,
                         signal_true){

  tn <- sum((signal_true == 0) & (signal_detect == 0))
  tp <- sum((signal_true == 1) & (signal_detect == 1))
  fp <- sum((signal_true == 0) & (signal_detect == 1))
  fn <- sum((signal_true == 1) & (signal_detect == 0))

  level <- ifelse( fp != 0, 1, 0)
  fdr <- ifelse( fp+tp == 0, 0, fp/(fp+tp) )
  power <- ifelse( (tp+fp) != 0, 1, 0)
  sen <- ifelse( tp+fn == 0, 0, tp/(tp+fn) )
  f <- ifelse( tp+1/2*(fp+fn) == 0, 0, tp/(tp+1/2*(fp+fn)) )

  return(
    c(Level = level,
      FDR = fdr,
      Power = power,
      Sensitivity = sen,
      F_score = f)
  )
}

simu_data <- function(lambda_0,
                      sig_row = NULL,
                      sig_col = NULL){


  n_row <- nrow(statin46)
  n_col <- ncol(statin46)

  signal_mat <- matrix(1, n_row, n_col)

  if (!is.null(sig_row)) {
    if (!is.null(sig_col)){
      signal_mat[sig_row, sig_col] <- lambda_0
    } else {
      signal_mat[sig_row, ] <- lambda_0
    }
  } else {
    if (!is.null(sig_col)){
      signal_mat[, sig_col] <- lambda_0
    }
  }

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

create_row <- function(model, method, lambda, repl, sig){

  metrics <- sig %>%
    sapply(
      function(x)signal.table(x, sig_true)
    )

  out <- tibble(
    model = model,
    method = method,
    lambda = lambda,
    repl = repl,
    Level = metrics[1,],
    FDR = metrics[2,],
    Power = metrics[3,],
    Sensitivity = metrics[4,],
    F_score = metrics[5,]
  )

  return(out)

}

# true signal matrix
sig_true <- simu_data(lambda_0 = 2,
                      sig_row = c(5, 45),
                      sig_col = c(1:6))$signal_true
# rds file name

file_names <- list.files(path = "D:/Documents/UB/Research/Code/CCR/Download/single_cell_signal",
                         pattern = "\\.RDS",
                         full.names = TRUE)

# create empty tibble
simu_df <- tibble(
  model = character(),
  method = character(),
  lambda = numeric(),
  repl = numeric(),
  Level = numeric(),
  FDR = numeric(),
  Power = numeric(),
  Sensitivity = numeric(),
  F_score = numeric()
)

# arrange data

for (i in 1:length(file_names)) {

  temp <- try(readRDS(file_names[[i]]))

  if( inherits(temp, "try-error")  ){
    next
  }

  if(temp$err == 1) {
    next
  }

  temp$res_pvlrt$sig[,7] <- 0

  simu_df <- simu_df %>%
    bind_rows(
      names(temp$res_pvbayes_est) %>%
        lapply(
          function(x){
            create_row(
              model = x,
              method = c("pFDR", "Naive"),
              lambda = temp$lambda_true,
              repl = temp$task_id,
              sig = list(temp$res_pvbayes_est[[x]]$sig_pfdr,
                         temp$res_pvbayes_est[[x]]$sig_naive)
            )
          }
        )
    ) %>% bind_rows(
      create_row(
        model = "LRT",
        method = "None",
        lambda = temp$lambda_true,
        repl = temp$task_id,
        sig = list(temp$res_pvlrt$sig)
      )
    )

}


simu_df_final <- simu_df %>%
  drop_na() %>%
  pivot_longer(
    cols = c(Level, FDR, Power, Sensitivity, F_score),
    names_to = "measures",
    values_to = "value"
  ) %>%
  group_by(model, method, lambda, measures) %>%
  summarise(
    value = mean(value)
  ) %>%
  mutate(
    method = factor(method, levels = c("None",
                                       "pFDR",
                                       "Naive")),
    model = factor(model, levels = c("LRT",
                                     "poisson",
                                     "poisson_indep",
                                     "zip",
                                     "zip_indep"))
  )



plot_metrics <- function(metric = NULL,
                         lwd = 1.2){

  simu_df_final %>%
    {if (is.null(metric)) {.}
      else {filter( ., measures %in% metric )}} %>%
    ggplot(
      aes(x = lambda,
          y = value,
          colour = model,
          linetype = method,
          group = interaction(model, method)
      ), size = 2
    ) +
    facet_wrap( ~measures ) +
    ylim(0,1) +
    geom_line(size = lwd) +
    scale_linetype_manual(
      values = c(
        "None" = "solid",
        "pFDR" = "solid",
        "Naive" = "dotted"
      )
    ) +
    scale_color_manual(
      values = c(
        "LRT" = "black",
        "poisson" = "blue1",
        "poisson_indep" = "lightskyblue2",
        "zip" = "yellow1",
        "zip_indep" = "orange1"
      )
    ) +
    geom_hline(yintercept = 0.05,
               lty = "dashed",
               color = "red",
               linewidth = lwd) +
    annotate(geom="text",
             x=4.1,
             y=0.07,
             label="0.05",
             color="red") +
    theme(
      legend.position = "top"
    )

}


plot_metrics() %>%
  plotly::ggplotly()

















