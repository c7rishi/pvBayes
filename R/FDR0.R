#' Basis function for computing FDR0
#' @importFrom data.table :=
#' @importFrom data.table %like%
#' @noRd
FDR0 <- function(lambda_draws,
                 test_stat = function(x) quantile(x, 0.05),
                 optim = FALSE,
                 alpha = NULL,
                 k = 1.0){

  res_pFDR <- pFDR0(lambda_draws = lambda_draws,
                    test_stat = test_stat,
                    optim = optim,
                    k = k)

  v <- lambda_draws %>%
    {posterior::Pr(.>1)} %>%
    unname()

  d <- res_pFDR$discovery %>% unname()

  D <- sum(d)

  FDR <- ifelse(D != 0, sum( (1 - v)*d ) / D,  0)


  return(
    c(
      res_pFDR,
      list(FDR = FDR)
    )
  )

}

utils::globalVariables(c("."))

