### Simulation function
### source("/projects/academic/chakrab2/"xinweihu/R_scripts/Simulation_function.r")

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



Simulation <- function(
    data,
    mod = c(
      "zip_horseshoe",
      "zip_horseshoe_LKJ"
    ),
    stan_chains = 4,
    stan_core = 4,
    stan_cov_rate = 0.1,
    stan_retry = 3
){

    temp_pvbayes <- mod %>%
    sapply( function(m){
      cat("-----------------",m,"-----------------\n")
      tryCatch(
        pvBayes::pvbayes(data,
                model = m,
                stan_chains = stan_chains,
                stan_seed = 1234,
                stan_iter_sampling = 1000,
                starting = "LRT",
                stan_parallel_chains = getOption("mc.cores", stan_core),
                stan_cov_rate = stan_cov_rate,
                retry = stan_retry
        ),
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

  signal_pos_list <- rep(list(signal_pos), (length(mod)+1)) %>%
    `names<-`(
      c(mod, "pvlrt")
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

  pvbayes_cov <- temp_pvbayes %>%
    lapply(
      function(res){
        tryCatch(
         res$MCMC_convergence[2]
          ,
          error = function(e){ e }
        )

      }
    ) %>% unlist()

  return(
    list(
      discovery = list(discovery),
      metrics = list(confusion),
      stan_div = ifelse(inherits(pvbayes_cov, "error"),1,0),
      error = err
    )

  )

}

