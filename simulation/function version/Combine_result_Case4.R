library(tidyverse)
library(pvLRT)
library(ggplot2)


### Arguments

folder <- "simulation_case4"
file_path0 <- "D:/Documents/UB/Research/Code/CCR/Download"
file_path1 <- file.path(file_path0, folder)


### tool functions

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

  level <- ifelse( tn+fp == 0, 0, fp/(tn+fp))
  fdr <- ifelse( fp+tp == 0, 0, fp/(fp+tp) )
  power <- ifelse( tp != 0, 1, 0)
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


create_row <- function(model, lambda, repl, discovery_matrices, sig_true, index){

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
    index = index,
    model = model,
    #method = method,
    lambda_1 = lambda[1],
    lambda_2 = lambda[2],
    lambda_3 = lambda[3],
    repl = repl,
    Level = metrics$Level,
    FDR = metrics$FDR,
    Power = metrics$Power,
    Sensitivity = metrics$Sensitivity,
    F_score = metrics$F_score
  )

  return(out)

}


# rds file name

file_names <- list.files(path = file_path1,
                         pattern = "\\.RDS",
                         full.names = TRUE)

# create empty tibble
simu_df <- tibble(
  index = numeric(),
  model = character(),
  #method = character(),
  lambda_1 = numeric(),
  lambda_2 = numeric(),
  lambda_3 = numeric(),
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
  simu_df <- simu_df %>%
    bind_rows(
      names(temp$res_pvbayes_est) %>%
        map(
          function(x){
            create_row(
              index = i,
              model = x,
              #method = c("pFDR", "Naive"),
              lambda = temp$lambda_true,
              repl = temp$task_id,
              discovery_matrices = list(temp$res_pvbayes_est[[x]]$sig_pfdr),
                         # temp$res_pvbayes_est[[x]]$sig_naive),
              sig_true = temp$signal_true
            )
          },
          .progress = TRUE
        )
    ) %>% bind_rows(
      create_row(
        index = i,
        model = "LRT",
        #method = "None",
        lambda = temp$lambda_true,
        repl = temp$task_id,
        discovery_matrices = list(temp$res_pvlrt$sig),
        sig_true = temp$signal_true
      )
    )

}


simu_df %>%
  saveRDS(
    file = paste(file_path1,"_full.RDS", sep = "")
  )

simu_df <- readRDS(paste(file_path1,"_full.RDS", sep = ""))



simu_df_final <- simu_df %>%
  drop_na() %>%
  pivot_longer(
    cols = c(Level, FDR, Power, Sensitivity, F_score),
    names_to = "measures",
    values_to = "value"
  ) %>%
  group_by(model, lambda_1, lambda_2, lambda_3, measures) %>%
  summarise(
    value = mean(value)
  ) %>%
  mutate(
    model = factor(model, levels = c("LRT",
                                     "poisson",
                                     "poisson_indep",
                                     "zip",
                                     "zip_indep"))
  ) %>% pivot_wider(
    names_from = measures,
    values_from = value
  )

simu_df_final %>%
  saveRDS(
    file = paste(file_path1,".RDS", sep = "")
  )

simu_df_final <- readRDS(paste(file_path1,".RDS", sep = ""))

plot_metrics <- function(metric = NULL,
                         lwd = 1.2){

  simu_df_final %>%
    {if (is.null(metric)) {.}
      else {filter( ., measures %in% metric )}} %>%
    ggplot(
      aes(x = lambda,
          y = value,
          colour = model
      ), size = 2
    ) +
    facet_wrap( ~measures ) +
    ylim(0,1) +
    geom_line(size = lwd) +
    scale_color_manual(
      values = c(
        "LRT" = "black",
        #"poisson" = "blue1",
        "poisson_indep" = "lightskyblue2",
        #"zip" = "yellow1",
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

###########################

simu_df %>%
  select(index, model, FDR) %>%
  pivot_wider(
    names_from = model,
    values_from = FDR
  ) %>%
  filter(poisson_indep != 0 & poisson != 0) %>%
  print(n = 40)


temp <- try(readRDS(file_names[[273]]))

temp$data %>% unname()

temp$signal_true %>% unname()

temp$lambda_true

###poisson
res_poisson <- temp$data %>% pvbayes(
  model = "poisson",
  stan_chains = 1
)

res_poisson$draws$lambda %>% unname()

est_poisson <- res_poisson$draws$lambda %>% pvbayes_est(
  test_stat = function(x) quantile(x, 0.05),
  alpha = 0.05
)

est_poisson$sig_pfdr %>% unname()

###poisson_indep
res_poisson_indep <- temp$data %>% pvbayes(
  model = "poisson_indep",
  stan_chains = 1
)

res_poisson_indep$draws$lambda %>% unname()

est_poisson_indep <- res_poisson_indep$draws$lambda %>% pvbayes_est(
  test_stat = function(x) quantile(x, 0.05),
  alpha = 0.05
)

est_poisson_indep$sig_pfdr %>% unname()

###poisson_indep2
res_poisson_indep2 <- temp$data %>% pvbayes(
  model = "poisson_indep2",
  stan_chains = 1
)

res_poisson_indep2$draws$lambda %>% unname()



est_poisson_indep2 <- res_poisson_indep2$draws$lambda %>% pvbayes_est(
  test_stat = function(x) quantile(x, 0.05),
  alpha = 0.05
)

est_poisson_indep2$sig_pfdr %>% unname()

est_poisson_indep2$k
est_poisson_indep2$test_stat %>% unname() %>% round(2)
###zip
res_zip <- temp$data %>% pvbayes(
  model = "zip",
  stan_chains = 1
)

res_zip$draws$lambda %>% unname()

est_zip <- res_zip$draws$lambda %>% pvbayes_est(
  test_stat = function(x) quantile(x, 0.05),
  alpha = 0.05
)

est_zip$sig_pfdr %>% unname()

###zip_indep
res_zip_indep <- temp$data %>% pvbayes(
  model = "zip_indep",
  stan_chains = 1
)

res_zip_indep$draws$lambda %>% unname()

est_zip_indep <- res_zip_indep$draws$lambda %>% pvbayes_est(
  test_stat = function(x) quantile(x, 0.05),
  alpha = 0.05
)

est_zip_indep$sig_pfdr %>% unname()






