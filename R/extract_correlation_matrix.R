#' Extract correlation matrix from pvBayes object
#' @param obj A pvBayes object with lamabda draws
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
extract_correlation_matrix <- function(obj,
                                       par,
                                       by_row = FALSE,
                                       method = "pearson",
                                       log_scale = FALSE,
                                       ...){
  # browser()
  draws_list <- extract_draws_list(obj, par)

  tmp <- draws_list %>%
    lapply(
      function(mat){
        if(log_scale){
          mat <- log(mat)
        }
        if (by_row) {
          out <- t(mat) %>%
            cor(method = method,...)
        } else {
          out <- mat %>%
            cor(method = method,...)
        }
        return(out)
      }
    )

  tmp_array <- array(
    NA,
    dim = c(length(tmp),
            dim(tmp[[1]])),
    dimnames = c(
      list(paste0("iter_", seq_along(tmp))),
      dimnames(tmp[[1]])
    )
  )

  for (ii in seq_along(tmp)) {
    tmp_array[ii, , ] <- tmp[[ii]]
  }


  return(posterior::rvar(tmp_array))
}

