## Simulation Case 2b
## 2 row signals with random position
## lambda value is supplied from simulation command in CCR
## zero inflation parameter omega is set to be 0.1


library(tidyverse)
library(pvBayes)
library(pvLRT)

job_id <- Sys.getenv("SLURM_JOB_ID") %>%
  as.numeric() %>%
  {
    if (is.na(.)){
      ## fake input for local test
      13225877
    } else{
      .
    }
  }

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID") %>%
  as.numeric() %>%
  {
    if (is.na(.)){
      ## fake input for local test
      20
    } else{
      .
    }
  }


cmd_args <- commandArgs(trailingOnly = TRUE) %>%
  {
    if (length(.) == 0){
      ## fake input for local test
      c("2",
        "simulation_case1")
    } else {
      .
    }
  }

lambda_true <- cmd_args[[1]] %>% as.numeric()

folder <- cmd_args[[2]] %>%
  tryCatch(
    .,
    error = function(e) {}
  )

is.server <- {Sys.info()["user"] == "xinweihu" & Sys.info()["sysname"] == "Linux"} %>%
  as.vector()

path_output0 <- file.path("/projects", "academic", "chakrab2", "xinweihu", "output")

if( !is.null(folder) ){

  path_output1 <- file.path(path_output0, folder)

  if( is.server ){
    if( !dir.exists(path_output1) ){
      dir.create(path_output1)
      cat("Folder created \n")
    } else{
      cat("Folder already exists \n")
    }
  } else {
    cat("Local test")
  }

} else{
  path_output1 <- path_output0
}

path <- ifelse(is.server, path_output1, "Local_test")

# signal_row_common <- c(1, 15, 30, 45)
signal_row_common <- c(30, 45)


output <- tibble(
  lambda = numeric(),
  seed = integer(),
  signal_pos = list(),
  signal_mat = list(),
  data = list(),
  discovery = list(),
  metrics = list(),
  error = integer()
)

confusion_matrix <- function(signal, discovery, discovery_global){

  signal_global <- max(signal)

  tn <- sum((signal == 0) & (discovery == 0))
  tp <- sum((signal == 1) & (discovery == 1))
  fp <- sum((signal == 0) & (discovery == 1))
  fn <- sum((signal == 1) & (discovery == 0))

  tn_global <- ((signal_global == 0) & (discovery_global == 0))
  tp_global <- ((signal_global == 1) & (discovery_global == 1))
  fp_global <- ((signal_global == 0) & (discovery_global == 1))
  fn_global <- ((signal_global == 1) & (discovery_global == 0))


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

for (i in 1:10) {

  seed_task <- job_id + task_id + i

  set.seed(seed_task)

  signal_row_pos <-
    lapply(
      1:length(signal_row_common), function(x){
        sample(
          setdiff(c(1:(nrow(statin46)-1)), signal_row_common),
          size = 6-length(signal_row_common)
        ) %>%
          c(signal_row_common) %>%
          sort()
      }
    )

  zi_protect <- list()
  for (r in seq_along(signal_row_pos)) {
    for (j in signal_row_pos[[r]]) {
      zi_protect <- c(zi_protect, list(c(j, r)))
    }
  }

  signal_pos <- matrix(0, nrow = 47, ncol = 7)
  for (r in seq_along(zi_protect)) {
    signal_pos[ zi_protect[[r]][1], zi_protect[[r]][2] ] <- 1
  }

  log_signal_mat <-
    matrix(rnorm(n = 47*7, mean = 0, sd = 0.00001), nrow = 47, ncol = 7)

  for (j in 1:10){

    log_signal_mat[signal_row_pos[[j]],j] <-
      log(as.numeric(lambda_true)) + log_signal_mat[signal_row_pos[[j]],j]

  }

  signal_mat <- log_signal_mat %>% exp()

  data <- r_contin_table_zip(
    n = 1,
    row_marginals = rowSums(statin46),
    col_marginals = colSums(statin46),
    signal_mat = signal_mat,
    omega_vec = c(rep(0.1, (ncol(statin46)-1)), 0),
    no_zi_idx =
      c(
        zi_protect,
        list(
          c(nrow(statin46), 0), # the entire last row
          c(0, ncol(statin46)) # the entire last column
        )
      )
  )[[1]]

  temp_pvbayes <- c(
    "zip_horseshoe",
    "zip_horseshoe_LKJ"
  ) %>%
    sapply( function(m){
      cat("-----------------",m,"-----------------\n")
      tryCatch(
        pvbayes(data,
                model = m,
                stan_chains = 4,
                stan_seed = 1234,
                stan_iter_sampling = 1000,
                starting = "LRT"),
        error = function(e){ e }
      )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
    )

  temp_pvbayes_est <- temp_pvbayes %>%
    lapply(
      function(res){
        tryCatch(
          pvbayes_est(
            pvbayes_obj = res
          )
          ,
          error = function(e){ e }
        )

      }
    )

  discovery <- temp_pvbayes_est %>%
    lapply(
      function(est){
        tryCatch(
          extract_discovery_mat(est)
          ,
          error = function(e){ e }
        )
      }
    )

  discovery_global <- temp_pvbayes_est %>%
    lapply(
      function(est){
        tryCatch(
          est$discovery_global
          ,
          error = function(e){ e }
        )
      }
    )

  temp_pvlrt <- tryCatch(
    {pvlrt(data) %>%
        extract_p_value_matrix()},
    error = function(e){ e }
  )

  temp_pvlrt_est <- temp_pvlrt %>%
    {ifelse(.<0.05, 1, 0)}
  temp_pvlrt_est[,7] <- 0

  discovery$pvlrt <- temp_pvlrt_est
  discovery_global$pvlrt <- max(temp_pvlrt_est)

  signal_pos_list <- rep(list(signal_pos),3) %>%
    `names<-`(
      c("zip_horseshoe",
        "zip_horseshoe_LKJ",
        "pvlrt")
    )

  confusion <- mapply(confusion_matrix,
                      signal_pos_list,
                      discovery,
                      discovery_global,
                      SIMPLIFY = FALSE)

  err <- ifelse(inherits(temp_pvbayes, "error")|
                  inherits(temp_pvbayes_est, "error")|
                  inherits(temp_pvlrt, "error"),
                1, 0)

  output <-
    add_row(
      output,
      lambda = lambda_true,
      seed = seed_task,
      signal_pos = list(signal_pos),
      signal_mat = list(signal_mat),
      data = list(data),
      discovery = list(discovery),
      metrics = list(confusion),
      error = err
    )


}


out_name <- glue::glue("{folder}_\\
  job{job_id}_\\
  task{task_id}")

if (path != "Local_test") {
  saveRDS(output,
          file = glue::glue(path,"/",out_name,".RDS"))
}


