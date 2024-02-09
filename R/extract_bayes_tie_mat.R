#' Extract Bayesian Type I error matrix from pvBayes_est object
#' @param pvBayes_est_obj A pvBayes_bayes object with lamabda draws
#' @returns
#' \itemize{
#'   \item k - Critical value for calculating pFDR.
#'   \item optim - Indicate if the k value is optimized.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#' }
#' @export
extract_bayes_tie_mat <- function(pvbayes_est_obj,
                                  ...){

  if (!("pvbayes_est" %in% class(pvbayes_est_obj))) {
    stop("A 'pvbayes_est' object is required!")
  }

  return(pvbayes_est_obj$BayesTIE)
}

