#' Basis function for transform MCMC samples to rvar object
#' @importFrom data.table :=
#' @importFrom data.table %like%

stan_output_to_rvar <- function(obj,...){

  #posterior
  stan_samps <- obj$lambda_draws
  stan_samps[, draw := 1:.N]
  I <- obj$contin_table_long

 #  tmp_rvar <-
 #    obj %>%
 #    .[grepl("draws", names(.))] %>%
 #    lapply(
 #      function(x)
 #        posterior::as_draws_df(x) %>%
 #        posterior::as_rvar() %>%
 #        .[[1]]
 #    )
 # %>%

  out <- obj %>%
    .[grepl("draws", names(.))] %>%
    .[!sapply(., is.null)] %>%
    lapply(
      function(x) {
        posterior::as_draws_df(x) %>%
          posterior::as_draws_rvars() %>%
          .[[1]]
      }
    )

  #
  #
  # browser()
  #
  # tmp <- data.table::melt(
  #   stan_samps,
  #   id.vars = "draw",
  #   variable.name = "lambda",
  #   value.name = "sample"
  # ) %>%
  #   data.table::merge.data.table(
  #     obj$contin_table_long,
  #     by = "lambda"
  #   )
  #
  # # tmp2 <- tmp[, .(draw, AE, Drug)] %>%
  # #   as.array()
  # #
  # # .SD = tmp[draw == 1]
  #
  #
  # nsim <- stan_samps[[1]] %>%  length()
  # I <- obj$contin_table_long$AE %>% unique() %>% length()
  # J <- obj$contin_table_long$Drug %>% unique() %>% length()
  #
  #
  # tmp_env <- new.env()
  # tmp_env$out_lambda_array <- array(
  #   NA, dim = c(nsim, I, J),
  #   dimnames = list(
  #     paste0("draw_", seq_len(nsim)),
  #     obj$contin_table_long$AE %>% unique(),
  #     obj$contin_table_long$Drug %>% unique()
  #   )
  # )
  #
  # samp_list <- tmp[
  #   ,
  #   {
  #     this_lambda_nm <- .SD$lambda %>% unique()
  #     tmp_env$out_lambda_array[, AE, Drug] <- stan_samps[[this_lambda_nm]]
  #     # out <- data.table::dcast(
  #     #   .SD,
  #     #   AE ~ Drug,
  #     #   value.var = "sample"
  #     # ) %>%
  #     #   {
  #     #     input <- .
  #     #     tmp_mat <- as.data.frame(.)
  #     #     tmp_mat$AE <- NULL
  #     #     tmp_mat %>%
  #     #       as.matrix() %>%
  #     #       `rownames<-`(input$AE)
  #     #   }
  #     # .(sample_mat = .(out))
  #   },
  #   keyby = .(AE, Drug),
  #   .SDcols = c("AE", "Drug", "lambda")
  # ]
  #
  #
  # samp_array <- array(NA, c(nsim, I, J))
  #
  # for (draw in 1:nsim) {
  #   samp_array[draw, , ] <- samp_list$sample_mat[[draw]]
  # }
  #
  # # browser()
  # samp_rvar <- posterior::rvar(samp_array)
  # dimnames(samp_rvar) <- dimnames(samp_list$sample_mat[[1]])
  #
  # samp_rvar

  # out = list(lambda = samp_rvar)

  out
}
