#' Basis function for computing pFDR
#' @param par_draws MCMC samples obtianed from `pvbayes`
#' @param par_est Bayesian estimates
#' @param optim logical. Use specified critical value or optimize pFDR(k)=alpha. If FALSE, then k must be specified
#' @param alpha Confidence level s.t. pFDR(k)=alpha
#' @param k Critical value
#' @export k Critical value for calculating pFDR
#' @export optim Indicate if the k value is optimized or specified
#' @export pFDR Positive false discovery rate value
#' @export BayesTIE Bayesian type-I error

pFDR0 <- function(par_draws, 
                  par_est, 
                  optim = FALSE, 
                  alpha = NULL, 
                  k){
  
  par_name <- str_extract(names(par_est)[1], "^[^\\[]+") 
  
  Pr <-
    par_draws %>%
    map_dbl( function(x) mean( x<=1 ) ) %>%
    map2_dbl(par_est, function(x,y) ifelse( y>k, x, 1) ) %>% 
    as_tibble(
      rownames = par_name
    ) %>%
    rename(BayesTIE = value) %>% 
    mutate(Est = par_est)
  
  pFDR <- 
    Pr$BayesTIE %>% 
    {.[.<1]} %>% 
    mean() %>% 
    {ifelse( is.na(.), 0, .)}
    
  return(
    list(
      k = k,
      optim = optim,
      pFDR = pFDR, 
      BayesTIE = Pr
    )
  )
  
}

