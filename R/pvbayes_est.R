#' Bayes estimate using MCMC samples
#' @param lambda_draws MCMC samples obtianed from `pvbayes`
#' @param test_method The choice of test statistic based on posterior distribution
#' @param alpha Confidence level s.t. pFDR(k)=alpha
#' @returns
#' \itemize{
#'   \item k - Critical value for calculating pFDR
#'   \item optim - Indicate if the k value is optimized or specified
#'   \item pFDR - Positive false discovery rate value
#'   \item BayesTIE - Bayesian type-I error
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)

#'}
#' @export
pvbayes_est <- function(lambda_draws,
                        test_stat,
                        prob = NULL,
                        alpha = .05,
                        ...){

  res <-
    pFDR(lambda_draws = lambda_draws,
         test_stat =  test_stat,
         optim = TRUE,
         alpha = alpha)

  sig_naive <- lambda_draws %>%
    posterior::quantile2(0.05) %>%
    {ifelse( .> 1, 1, 0)}

  return(
    c(
      res,
      list(
        sig_naive = sig_naive
      )
    )

  )

}
