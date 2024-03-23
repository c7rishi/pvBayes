#' Basis function for computing pFDR
#' @importFrom data.table :=
#' @importFrom data.table %like%
#' @noRd
pFDR0 <- function(lambda_draws,
                  test_stat = function(x) quantile(x, 0.05),
                  optim = FALSE,
                  alpha = NULL,
                  k = 1.0,
                  replace_prob_null_NA = FALSE) {
  # browser()
  posterior_prob_null <- mean(lambda_draws <= 1)
  test_stat_mat <- test_stat(lambda_draws)
  discovery_mat <- (1 * (test_stat_mat > k)) %>%
    ifelse(
      is.na(.),
      0,
      .
    )
  n_discovery <- sum(discovery_mat, na.rm = TRUE) %>%
    {
      if (is.na(.)) 0
      else if (is.infinite(.)) 0
      else .
    }
  pFDR <- if (n_discovery == 0) {
    0
  } else {
    sum(posterior_prob_null * discovery_mat)/n_discovery
  }

  return(
    list(
      k = k,
      optim = optim,
      pFDR = pFDR,
      test_stat = test_stat_mat,
      BayesTIE = posterior_prob_null,
      range_test_stat = range(test_stat_mat),
      discovery = 1 * discovery_mat
    )
  )
}


utils::globalVariables(c("."))

