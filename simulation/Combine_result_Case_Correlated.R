library(tidyverse)
library(pvLRT)
library(ggplot2)
library(posterior)
### Arguments

folder <- "simulation_correlated"
file_path0 <- "D:/Documents/UB/Research/Code/CCR/Download"
file_path1 <- file.path(file_path0, folder)

# rds file name

file_names <- list.files(path = file_path1,
                         pattern = "\\.RDS",
                         full.names = TRUE)

### tool functions

# library(caret)
# confusion_table(temp$res_pvbayes_est$poisson$sig_pfdr, sig_true)
# confusionMatrix(as.factor(temp$res_pvbayes_est$poisson$sig_pfdr),as.factor(sig_true))
# Verified to be good, the calculation of tn, tp, fp and fn are correct.
# Power is different from specificity bcz we defined power to be 1 if we detect at least
# one signal in one repetition.
confusion_table <- function(discovery,
                            discovery_global,
                            signal_true){

  #discovery_global <- max(discovery)

  signal_true_global <- max(signal_true)

  #browser()
  tn <- sum((signal_true == 0) & (discovery == 0))
  tp <- sum((signal_true == 1) & (discovery == 1))
  fp <- sum((signal_true == 0) & (discovery == 1))
  fn <- sum((signal_true == 1) & (discovery == 0))

  tn_global <- sum((signal_true_global == 0) & (discovery_global == 0))
  tp_global <- sum((signal_true_global == 1) & (discovery_global == 1))
  fp_global <- sum((signal_true_global == 0) & (discovery_global == 1))
  fn_global <- sum((signal_true_global == 1) & (discovery_global == 0))


  level <- ifelse( fp_global != 0, 1, 0)
  fdr <- ifelse( fp+tp == 0, 0, fp/(fp+tp) )
  power <- ifelse( tp_global != 0, 1, 0)
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


# create_row <- function(model, lambda, repl, discovery_matrices, sig_true, index, fix){
#
#   metrics <- discovery_matrices %>%
#     sapply(
#       function(x)confusion_table(x, sig_true)
#     ) %>% t() %>% as_tibble()
#
#   out <- tibble(
#     index = index,
#     model = model,
#     fix = fix,
#     lambda = lambda,
#     repl = repl,
#     Level = metrics$Level,
#     FDR = metrics$FDR,
#     Power = metrics$Power,
#     Sensitivity = metrics$Sensitivity,
#     F_score = metrics$F_score
#   )
#
#   return(out)
#
# }



# create empty tibble
# simu_df <- tibble(
#   index = numeric(),
#   model = character(),
#   fix = integer(),
#   lambda = numeric(),
#   repl = numeric(),
#   Level = numeric(),
#   FDR = numeric(),
#   Power = numeric(),
#   Sensitivity = numeric(),
#   F_score = numeric()
# )

# arrange data
#
# pb <- txtProgressBar( max = length(file_names),style = 3)
#
# for (i in 1:length(file_names)) {
#
#   temp <- try(readRDS(file_names[[i]]))
#
#   if( inherits(temp, "try-error")  ){
#     next
#   }
#
#   if(temp$err == 1) {
#     next
#   }
#
#   simu_df <- simu_df %>%
#     bind_rows(
#
#       lres <- names(temp$res_pvbayes_est) %>%
#         lapply(
#           function(x){
#             tibble(
#               index = i,
#               model = x,
#               lambda = temp$lambda_true,
#               repl = temp$task_id,
#               discovery_matrices = list(temp$res_pvbayes_est[[x]]$sig_pfdr),
#               sig_true = list(temp$signal_pos_mat),
#               fix = temp$fixed_num
#             )
#           }
#         ) %>%
#         bind_rows() %>%
#         bind_rows(
#           tibble(
#             index = i,
#             model = "LRT",
#             lambda = temp$lambda_true,
#             repl = temp$task_id,
#             discovery_matrices = list(temp$res_pvlrt$sig),
#             sig_true = list(temp$signal_pos_mat),
#             fix = temp$fixed_num
#           )
#         )
#
#
#
#
#
#     ) %>% bind_rows(
#       create_row(
#         index = i,
#         model = "LRT",
#         #method = "None",
#         lambda = temp$lambda_true,
#         repl = temp$task_id,
#         discovery_matrices = list(temp$res_pvlrt$sig),
#         sig_true = temp$signal_pos_mat,
#         fix = temp$fixed_num
#       )
#     )
#
#   setTxtProgressBar(pb, i)
# }


simu_df <- map(
  file_names,
  function(file_name){

    temp <- try(readRDS(file_name))

    if( inherits(temp, "try-error")  ){
      out <- NULL
    } else if(temp$err == 1) {
      out <- NULL
    } else {
      out <- names(temp$res_pvbayes_est) %>%
        lapply(
          function(x){

            global_discovery_ind <- tryCatch(
              (temp$res_pvbayes[[x]]$lambda %>%
                 draws_of() %>%
                 apply(1, function(x)max(x[,-7]) ) %>%
                 quantile(0.05)) > temp$res_pvbayes_est[[x]]$k,
              error = function(e) { e }
            )

            if (is(global_discovery_ind, "error") ){
              return(NULL)
            } else {
              tibble(
                index = i,
                model = x,
                lambda = temp$lambda_true,
                repl = temp$task_id,
                discovery_matrices = list(temp$res_pvbayes_est[[x]]$sig_pfdr),
                sig_true = list(temp$signal_pos_mat),
                fix = temp$fixed_num,
                discovery_global = global_discovery_ind * 1
              )
            }
          }
        ) %>%
        bind_rows() %>%
        bind_rows(
          tibble(
            index = i,
            model = "LRT",
            lambda = temp$lambda_true,
            repl = temp$task_id,
            discovery_matrices = list(temp$res_pvlrt$sig),
            sig_true = list(temp$signal_pos_mat),
            fix = temp$fixed_num,
            discovery_global = max(temp$res_pvlrt$sig)
          )
        )
    }

    return(out)
    # simu_df <- simu_df %>%
    #   bind_rows(
    #     names(temp$res_pvbayes_est) %>%
    #       lapply(
    #         function(x){
    #           create_row(
    #             index = i,
    #             model = x,
    #             lambda = temp$lambda_true,
    #             repl = temp$task_id,
    #             discovery_matrices = list(temp$res_pvbayes_est[[x]]$sig_pfdr),
    #             sig_true = temp$signal_pos_mat,
    #             fix = temp$fixed_num
    #           )
    #         }
    #       )
    #   ) %>% bind_rows(
    #     create_row(
    #       index = i,
    #       model = "LRT",
    #       #method = "None",
    #       lambda = temp$lambda_true,
    #       repl = temp$task_id,
    #       discovery_matrices = list(temp$res_pvlrt$sig),
    #       sig_true = temp$signal_pos_mat,
    #       fix = temp$fixed_num
    #     )
    #   )

  },
  .progress = TRUE
) %>% bind_rows()


simu_df <- simu_df %>%
  mutate(
    test_res = pmap(
      list(discovery_global, discovery_matrices, sig_true),
      function(ds_g, ds, s_t){
        confusion_table(
          discovery = ds,
          discovery_global = ds_g,
          signal_true = s_t
        ) %>%
          as.list() %>%
          as_tibble()
      }
    )
  ) %>% unnest(test_res)


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
  group_by(model, lambda, measures, fix) %>%
  summarise(
    value = mean(value)
  ) %>%
  mutate(
    model = factor(model, levels = c("LRT",
                                     "zip_horseshoe",
                                     "zip_horseshoe_LKJ"
    )),
    cor = factor(fix, levels = c(5,2), labels = c("Strong", "Weak"))
  ) %>%
  select(!fix) %>%
  group_by(model, lambda, measures, cor)

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
          colour = model,
          linetype = cor
      ), size = 2
    ) +
    facet_wrap( ~measures ) +
    ylim(0,1) +
    geom_line(size = lwd) +
    scale_color_manual(
      values = c(
        "LRT" = "black",
        "zip_horseshoe"  = "blue",
        "zip_horseshoe_LKJ" = "yellow"
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



mat_strong <- mat_weak <- matrix(0, nrow = 7, ncol = 7)

pb <- txtProgressBar( max = length(file_names),style = 3)

for (i in 1:length(file_names)) {

  temp <- try(readRDS(file_names[[i]]))

  if (temp$fixed_num == 5) {
    mat_strong <- mat_strong + temp$sample_corr_mat
  } else {
    mat_weak <- mat_weak + temp$sample_corr_mat
  }

  setTxtProgressBar(pb, i)
}

mat_strong <- mat_strong/2000
mat_weak <- mat_weak/2000

mat_strong %>% round(3)
mat_weak %>% round(3)

temp$res_pvbayes_est$zip_horseshoe_LKJ$k

temp$res_pvbayes$zip_horseshoe$lambda %>%
  draws_of() %>%
  apply(
    c(1),
    function(x){
      max(x[,-7])
    }
  ) %>%
  quantile(0.05)
posterior::rvar_max()

res_lrt <- temp$data %>% pvlrt()
res_lrt
