library(tidyverse)
library(pvLRT)
library(ggplot2)


### Arguments

sig_row <- "45"
sig_col <- "1"

folder <- "single_cell_signal_0823"
file_path0 <- "D:/Documents/UB/Research/Code/CCR/Download"
file_path1 <- file.path(file_path0, folder)


### tool functions

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

# library(caret)
# confusion_table(temp$res_pvbayes_est$poisson$sig_pfdr, sig_true)
# confusionMatrix(as.factor(temp$res_pvbayes_est$poisson$sig_pfdr),as.factor(sig_true))
# Verified to be good, the calculation of tn, tp, fp and fn are correct.
# Power is different from specificity bcz we defined power to be 1 if we detect at least
# one signal in one repetition.
confusion_table <- function(discovery,
                            signal_true){

  tn <- sum((signal_true == 0) & (discovery == 0))
  tp <- sum((signal_true == 1) & (discovery == 1))
  fp <- sum((signal_true == 0) & (discovery == 1))
  fn <- sum((signal_true == 1) & (discovery == 0))

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


create_row <- function(model, method, lambda, repl, discovery_matrices){

  metrics <- discovery_matrices %>%
    sapply(
      function(x)confusion_table(x, sig_true)
    ) %>% t() %>% as_tibble()

  # metrics <- sig %>%
  #   lapply(
  #     function(x)confusion_table(x, sig_true)
  #   )

  # lapply is not so convenience because combining from several lists is not easy to do

  out <- tibble(
    model = model,
    method = method,
    lambda = lambda,
    repl = repl,
    Level = metrics$Level,
    FDR = metrics$FDR,
    Power = metrics$Power,
    Sensitivity = metrics$Sensitivity,
    F_score = metrics$F_score
  )

  return(out)

}

# true signal matrix
sig_true <- simu_data(lambda_true = 2,
                      sig_row = sig_row,
                      sig_col = sig_col)$signal_true
# rds file name

file_names <- list.files(path = file_path1,
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
              discovery_matrices = list(temp$res_pvbayes_est[[x]]$sig_pfdr,
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
        discovery_matrices = list(temp$res_pvlrt$sig)
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

















