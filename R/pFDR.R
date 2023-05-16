#' Compute pFDR
#' @param lambda_draws MCMC samples obtianed from `pvbayes`
#' @param lambda_est Bayesian estimates
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
#' mod <- pvbayes(contin_table = statin46, model = "horseshoe")
#'
#' est <- lambda_draws %>% apply(2, median)
#'
#' # K is optimized
#' pFDR(par_draws = mod$lambda_draws, par_est = est)
#'
#' # K is specified
#' pFDR(par_draws = mod$lambda_draws, par_est = est, optim = FALSE, k = 1.1)
#' }
#' @export
pFDR <- function(lambda_draws,
                 lambda_est,
                 optim = TRUE,
                 alpha = .05,
                 k = NULL){


  temp <- function(x){

    res <- pFDR0(lambda_draws = lambda_draws,
                 lambda_est = lambda_est,
                 k = x)$pFDR - alpha

    return(res)

  }

  if (optim == TRUE){

    k.optim <- tryCatch(
      {stats::uniroot(temp,
                      interval = c(0,max(par_est))
      )$root},
      error = function(e){
        max(lambda_est)
      }
    )

    return(
      pFDR0(lambda_draws = lambda_draws,
            lambda_est = lambda_est,
            optim = TRUE,
            k = k.optim)
    )

  } else {

    if (is.null(k)){ stop("k must be specified!") }

    return(
      pFDR0(lambda_draws = lambda_draws,
            lambda_est = lambda_est,
            optim = FALSE,
            k = k)
    )

  }




}

