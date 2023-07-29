#' Compute pFDR
#' @param lambda_draws MCMC samples obtianed from `pvbayes`
#' @param test_stat Bayesian estimates
#' @param optim logical. Use specified critical value or optimize pFDR(k)=alpha. Default is TRUE, if FALSE then k must be specified
#' @param alpha Confidence level s.t. pFDR(k)=alpha
#' @param k Critical value
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
#' }
#' @export
pFDR <- function(lambda_draws,
                 test_stat,
                 optim = TRUE,
                 alpha = .05,
                 k = NULL){


  temp <- function(x){

    res <- pFDR0(lambda_draws = lambda_draws,
                 test_stat = test_stat,
                 k = x)$pFDR - alpha

    return(res)

  }

  range_test_stat <- pFDR0(lambda_draws = lambda_draws,
                           test_stat = test_stat,
                           k = 1)$range_test_stat

  if (optim == TRUE){

    k.optim <- tryCatch(
      {stats::uniroot(
        temp,
        interval = range_test_stat
      )$root},
      error = function(e){
        range_test_stat[2]
      }
      )

  } else {

    if (is.null(k)){ stop("k must be specified!") }

    k.optim <- k

  }

  return(
    pFDR0(lambda_draws = lambda_draws,
          test_stat = test_stat,
          optim = optim,
          k = k.optim)
  )

}

