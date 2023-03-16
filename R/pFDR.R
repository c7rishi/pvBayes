#' Computing pFDR
#' @param par_draws MCMC samples obtianed from `pvbayes`
#' @param par_est Bayesian estimates
#' @param optim logical. Use specified critical value or optimize pFDR(k)=alpha. Default is TRUE, if FALSE then k must be specified
#' @param alpha Confidence level s.t. pFDR(k)=alpha
#' @param k Critical value
#' @export optim Indicate if the k value is optimized or specified
#' @export pFDR Positive false discovery rate value
#' @export BayesTIE Bayesian type-I error
#' @export k Critical value for calculating pFDR
#' 
pFDR <- function(par_draws, 
                 par_est, 
                 optim = TRUE, 
                 alpha = .05, 
                 k = NULL){
  
  temp <- function(x){
    
    res <- pFDR0(par_draws = par_draws, 
                 par_est = par_est, 
                 k = x)$pFDR - alpha
    
    return(res)
    
  }
  
  if (optim == TRUE){
    
    k.optim <-
      uniroot(temp, 
              interval = c(0,max(par_est))
      )$root
    
    return(
      pFDR0(par_draws = par_draws,
            par_est = par_est,
            optim = TRUE,
            k = k.optim)
    ) 
    
  } else {
    
    if (is.null(k)){ stop("k must be specified!") }
    
    return(
      pFDR0(par_draws = par_draws,
            par_est = par_est,
            optim = FALSE,
            k = k)
    ) 
    
  }
  
  
  
  
}

