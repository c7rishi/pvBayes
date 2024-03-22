library(tidyverse)

case <- "Case2a"

tb_all <- readRDS(
  glue::glue(
    "D:/Documents/UB/Research/Code/CCR/Download/{case}_combine.RDS"
  )

  )

tb_long <- tb_all %>%
  select(!c(signal_pos, signal_mat, data, discovery, error)) %>%
  unnest_longer(
    col = metrics,
    values_to = "measure",
    indices_to = "model"
  ) %>%
  unnest_longer(
    col = measure,
    values_to = "value",
    indices_to = "measure"
  )

tb_summary <- tb_long %>%
  group_by(model, measure, lambda) %>%
  summarise(
    value = mean(value)
  )

plot_metrics <- function(metric = NULL,
                         lwd = 1.2){

  tb_summary %>%
    # {if (is.null(metric)) {.}
    #   else {filter( ., measures %in% metric )}} %>%
    ggplot(
      aes(x = lambda,
          y = value,
          colour = model#,
          #linetype = cor
      ), size = 2
    ) +
    facet_wrap( ~measure ) +
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
plot_metrics()

tb_long %>% filter(lambda == 2) %>%
  filter(measure == "Sensitivity") %>%
  pivot_wider(
    values_from = value,
    names_from = model
  ) %>% select(seed, zip_horseshoe, zip_horseshoe_LKJ, pvlrt) %>%
  filter(pvlrt < zip_horseshoe_LKJ)

data <- tb_all %>% filter(seed == 10005, lambda==2) %>%
  pull(data) %>% {.[[1]]}
res <- data %>%
  pvbayes(
    stan_chains = 1,model = "zip_horseshoe_LKJ_noridge",
    starting = "LRT"
  )

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
