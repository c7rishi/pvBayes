library(tidyverse)

case <- "case1"

tb_all <- readRDS(
  glue::glue(
    "D:/Documents/UB/Research/Code/CCR/Download/simulation_{case}_combine.RDS"
  )

  )

temp <- tb_all %>%
  select(!c(signal_pos, signal_mat, data, discovery, error)) %>%
  unnest_longer(
    col = metrics,
    values_to = "metrics",
    indices_to = "model"
  ) %>%
  unnest_longer(
    col = metrics,
    values_to = "value",
    indices_to = "metric"
  )


temp %>% filter(
  is.na(value)
)
temp %>%
  group_by(model, metric) %>%
  summarise(
    value = mean(value)
  )




