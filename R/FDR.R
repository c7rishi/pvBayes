#' Optimization of FDR
#' @param lambda_draws An rvar object of MCMC samples obtained from `pvbayes`.
#' @param test_stat A function of test statistic (e.g. mean, quantile).
#' @param optim A logical value. Use specified critical value or optimize pFDR(k)=alpha. Default is TRUE, otherwise k must be specified.
#' @param alpha A real value between (0,1). Confidence level s.t. pFDR(k)=alpha.
#' @param k A real value. Specified critical value.
#' @param n_eval An integer. Amount of grids when doing optimization.
#' @param thresh description
#' @returns
#' \itemize{
#'   \item k - Critical value for calculating pFDR.
#'   \item optim - Indicate if the k value is optimized.
#'   \item pFDR - Positive false discovery rate.
#'   \item test_stat - A matrix of test statistics in each cell.
#'   \item BayesTIE - A matrix of Bayesian type-I error.
#'   \item range_test_stat - A vector of maximal and minimal test statistics.
#'   \item sig_pfdr - A matrix of discoveries (0 or 1) using pFDR as the threshold.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#' }
#' @export
FDR <- function(lambda_draws,
                 test_stat,
                 optim = TRUE,
                 alpha = .05,
                 k = NULL,
                 n_eval = 100,
                 thresh = 1.05){

  temp <- function(x){

    res <- FDR0(lambda_draws = lambda_draws,
                 test_stat = test_stat,
                 k = x)$FDR - alpha

    return(res)

  }


  range_test_stat <- FDR0(lambda_draws = lambda_draws,
                           test_stat = test_stat,
                           k = 1)$range_test_stat

  if (optim == TRUE){

    k_val <- seq(from = max(range_test_stat[1], thresh),
                 to = min(range_test_stat[2], thresh),
                 length.out = n_eval)

    temp_val <- sapply(k_val, temp)

    k.optim <- ifelse(sum(temp_val<=0)==0,
                      range_test_stat[2],
                      min(k_val[temp_val<=0])
    )


  } else {

    if (is.null(k)){ stop("k must be specified!") }

    k.optim <- k

  }

  return(
    FDR0(lambda_draws = lambda_draws,
          test_stat = test_stat,
          optim = optim,
          k = k.optim)
  )

}

