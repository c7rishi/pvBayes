library(tidyverse)

folder <- "simulation_case1"
path_output0 <- file.path("/projects", "academic", "chakrab2", "xinweihu", "output")
path_output1 <- file.path(path_output0, folder)

path_output1 <- "D:/Documents/UB/Research/Code/CCR/Download/test"



file_names <- list.files(path = path_output1,
                         pattern = "\\.RDS",
                         full.names = TRUE)


tb_list <- list()
for (i in seq_along(file_names)){
  tb_list[[i]] <- readRDS(file_names[[i]])

}

tb_all <- do.call(rbind, tb_list)

tb_all_long <- tb_all %>% pivot_longer(
  cols = c(pvlrt, horseshoe, LKJ),
  values_to = "discovery",
  names_to = "model"
)

output %>%
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





pb <- txtProgressBar( max = length(file_names),style = 3)

for (i in 1:length(file_names)){}


