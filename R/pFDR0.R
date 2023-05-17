#' Basis function for computing pFDR
#' @importFrom data.table :=
#' @importFrom data.table %like%
#' @noRd
pFDR0 <- function(lambda_draws,
                  test_stat,
                  optim = FALSE,
                  alpha = NULL,
                  k){

  lambda_s1 <- lambda_draws  %>%
    {posterior::Pr(. <= 1)}

  BayesTIE <- lambda_s1 %>%
    {ifelse( test_stat > k, ., 1)}

  pFDR <- lambda_s1 %>%
    .[test_stat > k] %>%
    mean() %>%
    {ifelse( is.na(.), 0, .)}

  return(
    list(
      k = k,
      optim = optim,
      pFDR = pFDR,
      test_stat = test_stat,
      BayesTIE = BayesTIE
    )
  )

}

utils::globalVariables(c("."))

