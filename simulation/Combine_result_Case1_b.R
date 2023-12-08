library(tidyverse)
library(pvLRT)
library(ggplot2)
library(posterior)
### Arguments

folder <- "simulation_case1b"
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

  #browser()
  tn <- sum((signal_true == 0) & (discovery == 0))
  tp <- sum((signal_true == 1) & (discovery == 1))
  fp <- sum((signal_true == 0) & (discovery == 1))
  fn <- sum((signal_true == 1) & (discovery == 0))

  level <- ifelse( fp != 0, 1, 0)
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


# rds file name

file_names <- list.files(path = file_path1,
                         pattern = "\\.RDS",
                         full.names = TRUE)

# create empty tibble
simu_df <- tibble(
  index = numeric(),
  model = character(),
  #method = character(),
  lambda = numeric(),
  repl = numeric(),
  Level = numeric(),
  FDR = numeric(),
  Power = numeric(),
  Sensitivity = numeric(),
  F_score = numeric()
)

# arrange data

pb <- txtProgressBar( max = length(file_names),style = 3)

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
        lapply(
          function(x){
            create_row(
              index = i,
              model = x,
              #method = c("pFDR", "Naive"),
              lambda = temp$lambda_true,
              repl = temp$task_id,
              discovery_matrices = list(temp$res_pvbayes_est[[x]]$sig_pfdr),
              # temp$res_pvbayes_est[[x]]$sig_naive),
              sig_true = temp$signal_pos_mat
            )
          }
        )
    ) %>% bind_rows(
      create_row(
        index = i,
        model = "LRT",
        #method = "None",
        lambda = temp$lambda_true,
        repl = temp$task_id,
        discovery_matrices = list(temp$res_pvlrt$sig),
        sig_true = temp$signal_pos_mat
      )
    )

  setTxtProgressBar(pb, i)
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
  group_by(model, lambda, measures) %>%
  summarise(
    value = mean(value)
  ) %>%
  mutate(
    model = factor(model, levels = c("LRT",
                                     "poisson_test",
                                     "zip_test",
                                     "poisson_indep_test",
                                     "zip_indep_test",
                                     "poisson_correlated_test",
                                     "poisson_LKJ_test"))
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
        "poisson_test" = "blue1",
        "poisson_indep_test" = "lightskyblue2",
        "zip_test" = "yellow1",
        "zip_indep_test" = "orange1",
        "poisson_LKJ_test" = "lightskyblue4",
        "poisson_correlated_test" = "red"
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



simu_df %>%
  filter(
    lambda == 1
  ) %>%
  select(index, model, FDR) %>%
  pivot_wider(
    names_from = model,
    values_from = FDR
  ) %>%
  filter(
    LRT == 0 &
      poisson_indep > 0 &
      zip_indep > 0 &
      poisson_indep2 > 0 &
      poisson_indep3 > 0
  ) %>%
  select(
    LRT,
    poisson_indep,
    poisson_indep2,
    poisson_indep3,
    poisson_LKJ
  )

temp <- try(readRDS(file_names[[402]]))
temp$lambda_true

temp$signal_pos_mat
temp$lambda_true_mat %>% round(2)

temp$res_pvbayes_est$poisson_indep3$sig_pfdr %>% unname
temp$res_pvbayes_est$poisson_indep3$test_stat %>% unname
temp$res_pvbayes_est$poisson_indep3$k

temp$res_pvbayes$poisson_indep3$lambda %>% unname()

res_est <- temp$res_pvbayes$poisson_indep2$lambda %>%
  pvBayes::pvbayes_est(
    test_stat = function(x) quantile(x, 0.05)
  )

res_est$sig_pfdr %>% unname

res_lrt <- temp$data %>% pvLRT::pvlrt()
res_lrt %>% pvLRT::extract_lrstat_matrix() %>% unname()

E <- temp$data %>%
  {tcrossprod(rowSums(.),colSums(.))/sum(.)}

temp$res_pvbayes$poisson_indep2$lambda2 <-
  with(
    temp$res_pvbayes$poisson_indep2,
    -(lambda-1)*E + temp$data * log(lambda)
  )

temp$res_pvbayes$poisson_indep2$lambda2 %>% unname()





temp$res_pvbayes$poisson_indep2$lambda_max <- with(
  temp$res_pvbayes$poisson_indep2,
  lambda %>%
    draws_of() %>%
    .[ , -47, ] %>%
    apply(c(1, 3), max) %>%
    as.matrix()
) %>%
  {
    arr <- array(NA, dim = c(nrow(.), 1, ncol(.)))
    arr[, 1, ] <- .
  } %>%
  posterior::rvar()


out <- temp$res_pvbayes$poisson_indep2$lambda_max %>%
  pvbayes_est(
    test_stat = function(x){
      quantile(x, 0.05)
    },
    thresh = 1.05
  )

res_est <- temp$res_pvbayes$poisson_indep2$lambda %>%
  pvbayes_est(
    test_stat = function(x){
      quantile(x, 0.05)
    },
    thresh = 1.05
  )
res_est$sig_pfdr %>% unname

temp$data %>% unname()



