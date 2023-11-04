#' Signal detection function using MCMC samples with pDFR threshold
#' @param lambda_draws An rvar object of MCMC samples obtained from `pvbayes`.
#' @param test_stat A function of test statistic (e.g. mean, quantile).
#' @param alpha A real value between (0,1). Confidence level s.t. pFDR(k)=alpha.
#' @returns
#' \itemize{
#'   \item k - Critical value for calculating pFDR.
#'   \item optim - Indicate if the k value is optimized.
#'   \item pFDR - Positive false discovery rate.
#'   \item test_stat - A matrix of test statistics in each cell.
#'   \item BayesTIE - A matrix of Bayesian type-I error.
#'   \item range_test_stat - A vector of maximal and minimal test statistics.
#'   \item sig_pfdr - A matrix of discoveries (0 or 1) using pFDR as the threshold.
#'   \item sig_naive - A reference matrix of discoveries (0 or 1) using .05 quantile.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#'}
#' @export
pvbayes_est <- function(lambda_draws,
                        test_stat,
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
