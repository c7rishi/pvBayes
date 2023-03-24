#' Basis function for computing pFDR
pFDR0 <- function(par_draws,
                  par_est,
                  optim = FALSE,
                  alpha = NULL,
                  k){

  par_name <- stringr::str_extract(names(par_est)[1], "^[^\\[]+")

  Pr <-
    par_draws %>%
    purrr::map_dbl( function(x) mean( x<=1 ) ) %>%
    purrr::map2_dbl(par_est, function(x,y) ifelse( y>k, x, 1) ) %>% #mapply
    dplyr::as_tibble(
      rownames = par_name
    ) %>%
    dplyr::rename(BayesTIE = value) %>%
    dplyr::mutate(Est = par_est)

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

