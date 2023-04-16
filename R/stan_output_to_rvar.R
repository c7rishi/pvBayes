#' Basis function for transform MCMC samples to rvar object
#' @importFrom data.table :=
#' @importFrom data.table %like%

stan_output_to_rvar <- function(obj,...){

  #posterior
  stan_samps <- obj$lambda_draws
  stan_samps[, draw := 1:.N]
  I <- obj$contin_table_long

  browser()

  tmp <- data.table::melt(
    stan_samps,
    id.vars = "draw",
    variable.name = "lambda",
    value.name = "sample"
  ) %>%
    data.table::merge.data.table(
      obj$contin_table_long,
      by = "lambda"
    )

  # tmp2 <- tmp[, .(draw, AE, Drug)] %>%
  #   as.array()
  #
  # .SD = tmp[draw == 1]

  samp_list <- tmp[
    ,
    {
      out <- data.table::dcast(
        .SD,
        AE ~ Drug,
        value.var = "sample"
      ) %>%
        {
          input <- .
          tmp_mat <- as.data.frame(.)
          tmp_mat$AE <- NULL
          tmp_mat %>%
            as.matrix() %>%
            `rownames<-`(input$AE)
        }
      .(sample_mat = .(out))
    },
    by = draw
  ]

  nsim <- samp_list$draw %>% unique() %>% length()
  I <- obj$contin_table_long$AE %>% unique() %>% length()
  J <- obj$contin_table_long$Drug %>% unique() %>% length()

  samp_array <- array(NA, c(nsim, I, J))

  for (draw in 1:nsim) {
    samp_array[draw, , ] <- samp_list$sample_mat[[draw]]
  }

  # browser()
  samp_rvar <- posterior::rvar(samp_array)
  dimnames(samp_rvar) <- dimnames(samp_list$sample_mat[[1]])

  samp_rvar

  # out = list(lambda = samp_rvar)
}
