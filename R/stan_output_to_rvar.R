#' Basis function for transform MCMC samples to rvar object
#' @importFrom data.table :=
#' @importFrom data.table %like%

stan_output_to_rvar <- function(obj,...){

  #posterior
  stan_samps <- obj$lambda_draws
  stan_samps[, draw := 1:.N]
  I <- obj$contin_table_long

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

  out

}
