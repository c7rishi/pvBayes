#' Generate signal matrix
#' @param I An integer of row number.
#' @param J An integer of column number.
#' @param sig_position A list of pairs (i, j) where signals are placed. To specify the entirety i-th row (or j-th column) use c(i, 0) (or c(0, j)).
#' @param sig_strength A list of signal strength corresponding to the position specified in `sig_position`.
#' @param non_sig_corr_mat A (J-1)x(J-1) matrix. Correlation matrix of multivariate-normal distribution for non-signal cells generation.
#' @param sig_mean_vec A positive numeric vector length by (J-1). Mean vector of multivariate-normal distribution for non-signal cells generation.
#' @param sig_corr_mat A (J-1)x(J-1) matrix. Correlation matrix of multivariate-normal distribution for signal cells generation.
#' @param var_vec A positive numeric vector length by (J-1). Variance vector of multivariate-normal for both signal and non-signal generation. Suggested to be relatively small value.
#' @returns
#' \itemize{
#'   \item lambda - A matrix of generated signals.
#'   \item signal_true - A matrix of (0,1) to indicate the position of signals.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#'}
#' @export
signal_generation <- function(
    I,
    J,
    sig_position,
    sig_strength,
    non_sig_corr_mat = matrix(0, nrow = J-1, ncol = J-1) %>% `diag<-`(1),
    # sig_mean_vec = NULL,
    sig_corr_mat = NULL,
    var_vec = rep(0.0001, J-1)
){


  #browser()

  signal_sub_mat <- matrix(0, nrow = I-1, ncol = J-1)
  strength_sub_mat <- matrix(0.99, nrow = I-1, ncol = J-1)

  for (k in 1:length(sig_position)) {
    if (sig_position[[k]][1] == 0 ){

      signal_sub_mat[,sig_position[[k]][2]] <- 1
      strength_sub_mat[,sig_position[[k]][2]] <- sig_strength[[k]]

    } else if ( sig_position[[k]][2] == 0 ) {

      signal_sub_mat[sig_position[[k]][1],] <- 1
      strength_sub_mat[sig_position[[k]][1],] <- sig_strength[[k]]

    } else {

      signal_sub_mat[sig_position[[k]][1], sig_position[[k]][2]] <- 1
      strength_sub_mat[sig_position[[k]][1], sig_position[[k]][2]] <- sig_strength[[k]]
    }

  }

  sig_lambda_sub_mat <- strength_sub_mat %>% apply(
    MARGIN = 1,
    FUN = function(mean_vec){
      mvtnorm::rmvnorm(
        n = 1,
        mean = log(mean_vec),
        sigma = diag(sqrt(var_vec)) %*% sig_corr_mat %*% diag(sqrt(var_vec))
      ) %>%
        exp()
    }
  ) %>%
    t()

  for (i in 1:(I-1)) {
    for (j in 1:(J-1)) {

      sig_lambda_sub_mat[i,j] <- ifelse(strength_sub_mat[i,j] <=1,
                                        min(sig_lambda_sub_mat[i,j], 0.99),
                                        max(sig_lambda_sub_mat[i,j], strength_sub_mat[i,j]))

    }
  }

  signal_sub_mat <- ifelse(sig_lambda_sub_mat<=1, 0, 1)

  signal_mat <- signal_sub_mat %>%
    rbind( rep(0, J-1) ) %>%
    cbind( rep(0, I) )

  # lambda_sub_mat <- ifelse(
  #   signal_sub_mat == 1,
  #   sig_lambda_sub_mat,
  #   non_sig_lambda_sub_mat
  # )

  lambda_mat <- sig_lambda_sub_mat %>%
    rbind( rep(0.99, J-1) ) %>%
    cbind( rep(0.99, I) )

  return(
    list(
      lambda = lambda_mat,
      signal_true = signal_mat
    )
  )

}


