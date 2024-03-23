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
pvbayes_est <- function(pvbayes_obj,
                        test_stat = function(x){quantile(x, 0.05)},
                        alpha = .05,
                        thresh = 1.1,
                        m = "FDR"){


  if (!("pvbayes" %in% class(pvbayes_obj))) {
    stop("A 'pvbayes' object is required!")
  }

  # if (m == "pFDR") {

  out <-
    pFDR(lambda_draws = pvbayes_obj$draws$lambda,
         test_stat =  test_stat,
         optim = TRUE,
         alpha = alpha,
         thresh = thresh)

  # } else if (m == "FDR") {
  #
  #   out <-
  #     FDR(lambda_draws = pvbayes_obj$draws$lambda,
  #         test_stat = test_stat,
  #         optim = TRUE,
  #         thresh = thresh,
  #         alpha = alpha)
  #
  #
  # }

  if (sum(out$discovery) > 0) {
    discovery_global <- (
      (pvbayes_obj$draws$lambda %>%
         posterior::draws_of() %>%
         apply(1, function(x)max(x[,-ncol(out$test_stat)]) ) %>%
         quantile(0.05) %>% unname()) > ifelse(is.na(out$k), max(pvbayes_obj$draws$lambda), out$k )
    ) * 1
  } else{
    discovery_global <- 0
  }



  # discovery_naive <- pvbayes_obj$draws$lambda %>%
  #   posterior::quantile2(0.05) %>%
  #   {ifelse( .> 1, 1, 0)}

  out <- c(
    out,
    list(
      discovery_global = discovery_global#,
      #discovery_naive = discovery_naive
    )
  )

  class(out) <- "pvbayes_est"

  return(
    out
  )

}
