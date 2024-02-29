library(tidyverse)

cmd_args <- commandArgs(trailingOnly = TRUE)

Case <- cmd_args[[1]]

path_output0 <- file.path("/projects", "academic", "chakrab2", "xinweihu", "output")
path_output1 <- file.path(path_output0, Case)

# path_output1 <- "D:/Documents/UB/Research/Code/CCR/Download"

file_names <- list.files(path = path_output1,
                         pattern = "\\.RDS",
                         full.names = TRUE)

tb_list <- list()
for (i in seq_along(file_names)){
  tb_list[[i]] <- readRDS(file_names[[i]])

}

tb_all <- do.call(rbind, tb_list)

saveRDS(
  tb_all,
  file = glue::glue("{path_output1}_combine.RDS")
)

