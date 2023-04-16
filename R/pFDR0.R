#' Basis function for computing pFDR
#' @importFrom data.table :=
#' @importFrom data.table %like%
pFDR0 <- function(par_draws,
                  par_est,
                  optim = FALSE,
                  alpha = NULL,
                  k){

  par_name <- sub("\\[.*", "", names(par_est)[1])

  Pr <-
    par_draws %>%
    apply(2, function(x) mean( x<=1 ) ) %>%
    mapply( function(x, y) ifelse( y>k, x, 1), ., par_est) %>%
    data.table::as.data.table(keep.rownames=TRUE) %>%
    data.table::setnames(old = c("rn", "."),
                         new = c(par_name,"BayesTIE")) %>%
    {.[ , Est := par_est]}

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

