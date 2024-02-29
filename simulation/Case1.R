library(tidyverse)
library(pvLRT)
case <- "Case1"

tb_all <- readRDS(
  glue::glue(
    "D:/Documents/UB/Research/Code/CCR/Download/{case}_combine.RDS"
  )

)

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
        "zip_horseshoe_LKJ" = "yellow",
        "zip_horseshoe_LKJ2" = "orange"
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


tb_sum <- tb_all  %>%
  select(seed, lambda, metrics) %>%
  mutate(
    metrics_df = map(metrics, \(x){
      x %>%
        map( \(y) y %>% as.list() %>% as_tibble() ) %>%
        bind_rows(.id = "model")
    }, .progress = TRUE)
  ) %>% select(-metrics) %>%
  unnest(metrics_df)

temp2 <- tb_sum %>%
  select(seed, lambda, model, Sensitivity) %>%
  filter( lambda == 1.5) %>%
  pivot_wider(
    names_from = model,
    values_from = Sensitivity
  ) %>%
  mutate(
    across(
      where(is.list),
      \(x){sapply(x, \(y)y[1])}
    )
  ) %>% select(
    !zip_horseshoe_LKJ2
  )

 temp2 %>% filter(
    pvlrt > zip_horseshoe
  ) %>%
   print(n = 50)

 data <- tb_all %>%
   filter(seed == 1002&lambda==1.5) %>%
   {.$data[[1]]}


