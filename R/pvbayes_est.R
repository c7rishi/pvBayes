#' Bayes estimate using MCMC samples
#' @param lambda_draws MCMC samples obtianed from `pvbayes`
#' @param est_quantile Bayesian estimates using pth-quantile of the MCMC samples
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
#' mod <- pvbayes(contin_table = statin46, model = "horseshoe")
#'
#' #obtain the MCMC samples
#' mod$lambda_draws
#'
#' pvbayes_est(mod$lambda_draws, .5)
#'}
#' @export
pvbayes_est <- function(lambda_draws,
                         est_quantile,
                         alpha = .05){

  lambda_est <- lambda_draws %>%
    posterior::as_draws_rvars() %>%
    .$lambda %>%
    posterior::quantile2(est_quantile)

  res <-
    pFDR(lambda_draws = lambda_draws,
          lambda_est =  lambda_est,
          optim = TRUE,
          alpha = alpha)

  return(
    res
  )

}

