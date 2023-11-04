#' Generate data set
#' @param row_marginals,col_marginals A (possibly named) vector of row and column marginal totals. Must add up to the same total. If named, the names are passed on to the randomly generated matrix/matrices.
#' @param lambda_mat A numeric matrix of dimension length(row_marginals) x length(col_marginals). The (i, j)-th entry of lambda_mat determines the signal strength of the i-th adverse event and j-th drug pair.
#' @param signal_position A matrix of (0,1) to indicate the position of signals.
#' @param omega_vec A vector of zero inflation probabilities. Must be of the same length as col_marginals.
#' @returns
#' \itemize{
#'   \item k - Critical value for calculating pFDR
#'   \item optim - Indicate if the k value is optimized or specified
#' }
#' @examples
#' \dontrun{
#'
#'
#' }
#' @export
data_generation <- function(
    row_marginals = rowSums(pvLRT::statin46),
    col_marginals = colSums(pvLRT::statin46),
    lambda_mat,
    signal_position,
    omega_vec = rep(0, length(col_marginals))
){

  #browser()

  n_row <- length(row_marginals)
  n_col <- length(col_marginals)

  zi_protect <- which(signal_position == 1, arr.ind = TRUE) %>%
    { mapply(c, as.numeric(.[,"row"]), as.numeric(.[,"col"]), SIMPLIFY = FALSE) } %>%
    c(.,
      list(
        c(n_row, 0),
        c(0, n_col)
      )
    )

  data <- pvLRT::r_contin_table_zip(
    n = 1,
    row_marginals = row_marginals,
    col_marginals = col_marginals,
    signal_mat = lambda_mat,
    omega_vec = omega_vec,
    no_zi_idx = zi_protect
  )[[1]]


  return(
    list(
      contin_table = data,
      signal_true = signal_position
    )
  )

}



