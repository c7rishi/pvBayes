library(tidyverse)
library(pvLRT)
case <- "Case1"

tb_all <- readRDS(
  glue::glue(
    "D:/Documents/UB/Research/Code/CCR/Download/{case}_combine.RDS"
  )

)
#
#
# tb_long %>%
#   filter(metric == "Level" & lambda == 1) %>%
#   select(seed, value, model) %>%
#   pivot_wider(
#     names_from = model,
#     values_from = value
#   ) %>%
#   select(seed, zip_horseshoe_LKJ, pvlrt) %>%
#   filter(
#     pvlrt == 0, zip_horseshoe_LKJ == 1
#   )
#
# data <- tb_all %>%
#   filter(
#     lambda == 1, seed == 50001
#   ) %>%
#   {.$data[[1]]}
#
# tb_all %>%
#   filter(
#     lambda == 1, seed == 50001
#   ) %>%
#   {.$discovery[[1]][["zip_horseshoe_LKJ"]]} %>%
#   unname()
#
# tb_all %>%
#   filter(
#     lambda == 1, seed == 50001
#   ) %>%
#   {.$discovery_global[[1]]} %>%
#   unname()
#
# tb_all %>%
#   filter(
#     lambda == 1, seed == 50001
#   ) %>%
#   {.$metrics[[1]][["zip_horseshoe_LKJ"]]}
#
# res_lrt <- data %>%
#   pvlrt()
#
# res_lrt %>% extract_p_value_matrix() %>% unname()
#
#
# res_lkj <- data %>% pvbayes(
#   "zip_horseshoe_LKJ",
#   stan_chains = 1,
#   stan_iter_sampling = 1000
# )
#
# res_lkj$draws$lambda %>% unname()
#
# est <- res_lkj %>% pvbayes_est()
# est$discovery %>% unname()
# est$discovery_global

tb_long <- tb_all %>%
  select(!c(signal_pos, signal_mat, data, discovery, error)) %>%
  unnest_longer(
    col = metrics,
    values_to = "metircs",
    indices_to = "model"
  ) %>%
  unnest_longer(
    col = metircs,
    values_to = "value",
    indices_to = "metric"
  )



tb_summary <- tb_long %>%
  group_by(model, metric, lambda) %>%
  filter(
    model != "zip_horseshoe_LKJ2"
  ) %>%
  summarise(
    value = mean(value)
  )

plot_metrics <- function(  lwd = 1.2){

  tb_summary %>%
    ggplot(
      aes(x = lambda,
          y = value,
          colour = model#,
          #linetype = cor
      ), size = 2
    ) +
    facet_wrap( ~metric ) +
    ylim(0,1) +
    geom_line(size = lwd) +
    scale_color_manual(
      values = c(
        "pvlrt" = "black",
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

plot.case1 <- plot_metrics()

tiff("D:/Documents/UB/Research/Code/pvBayes/simulation/Case1.tiff",
     units="in", width=7, height=5, res=300)
plot_metrics()
dev.off()

# tb_sum <- tb_all  %>%
#   select(seed, lambda, metrics) %>%
#   mutate(
#     metrics_df = map(metrics, \(x){
#       x %>%
#         map( \(y) y %>% as.list() %>% as_tibble() ) %>%
#         bind_rows(.id = "model")
#     }, .progress = TRUE)
#   ) %>% select(-metrics) %>%
#   unnest(metrics_df)
#
# temp2 <- tb_sum %>%
#   select(seed, lambda, model, Sensitivity) %>%
#   filter( lambda == 1.5) %>%
#   pivot_wider(
#     names_from = model,
#     values_from = Sensitivity
#   ) %>%
#   mutate(
#     across(
#       where(is.list),
#       \(x){sapply(x, \(y)y[1])}
#     )
#   ) %>% select(
#     !zip_horseshoe_LKJ2
#   )
#
#  temp2 %>% filter(
#     pvlrt > zip_horseshoe
#   ) %>%
#    print(n = 50)
#
 data <- tb_all %>%
   filter(seed == 1002&lambda==1.5) %>%
   {.$data[[1]]}

res1 <- res
res <- data %>%
  pvbayes(
    stan_chains = 1,model = "zip_horseshoe_LKJ_noridge",
    starting = "LRT"
  )

E <- res$E

res_lrt <- data %>% pvlrt()
res_lrt %>% extract_p_value_matrix() %>% {(.<0.05)*1} %>% unname()

est_pFDR <- res1 %>%  pvbayes_est(
  m = "pFDR",
  # test_stat = function(x){
  #   lambdahat <- posterior::quantile2(x, 0.5)
  #   -(lambdahat - 1) * E + data * log(lambdahat)
  # },
  # thresh = 1e-5,
  thresh = 1.01,
  test_stat = function(x) {
    posterior::quantile2(x, probs = 0.05)
  }
)
est_pFDR$discovery %>% unname()

est_FDR <- res %>%  pvbayes_est(m = "FDR", thresh = 1.01)

est_FDR$discovery %>% unname()
