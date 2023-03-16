#' Bayes estimate using MCMC samples
#' @param lambda_draws MCMC samples obtianed from `pvbayes`
#' @param est_quantile Bayesian estimates using pth-quantile of the MCMC samples
#' @param alpha Confidence level s.t. pFDR(k)=alpha
#' @export k Critical value for calculating pFDR
#' @export optim Indicate if the k value is optimized or specified
#' @export pFDR Positive false discovery rate value
#' @export BayesTIE Bayesian type-I error

pvbayes_est <- function(lambda_draws, 
                        est_quantile, 
                        alpha = .05){
  
  lambda_est <-
    lambda_draws %>%
    map_dbl(function(x)quantile(x,est_quantile))
  
  res <- 
    pFDR(par_draws = lambda_draws,
         par_est =  lambda_est,
         optim = TRUE,
         alpha = alpha)

  return(
    res
  )
  
}
