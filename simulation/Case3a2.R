library(tidyverse)

case <- "Case3a2"

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

tb_summary_3a <- tb_summary
tb_summary_3a$case <- "Strong"
