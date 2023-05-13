#' Basis function for computing pFDR
#' @importFrom data.table :=
#' @importFrom data.table %like%
#' @noRd
pFDR0 <- function(lambda_draws,
                   lambda_est,
                   optim = FALSE,
                   alpha = NULL,
                   k){

  BayesTIE <- lambda_draws %>%
    posterior::as_draws_rvars() %>%
    .$lambda %>%
    {posterior::Pr(. <= 1)} %>%
    {ifelse( lambda_est>k, ., 1)}

  pFDR <- BayesTIE %>%
    {.[.<1]} %>%
    mean() %>%
    {ifelse( is.na(.), 0, .)}

  return(
    list(
      k = k,
      optim = optim,
      pFDR = pFDR,
      lambda_est = lambda_est,
      BayesTIE = BayesTIE
    )
  )

}

utils::globalVariables(c("."))

