#' Extract lambda draws from pvBayes object and transform to a list form
#' @param obj A pvBayes object with lambda draws
#' @param par Name of draws transform to list
#' @returns
#' \itemize{
#'   \item A list of lambda draws, each element is one MCMC sample.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#' }
#' @export
extract_draws_list <- function(pvbayes_obj, par = "lambda"){

  if (!("pvbayes" %in% class(pvbayes_obj))) {
    stop("A 'pvbayes' object is required!")
  }

  dt_draws <- pvbayes_obj$draws[[par]] %>%
    posterior::as_draws_df() %>%
    data.table::as.data.table() %>%
    data.table::melt(
      id.vars = c(".chain", ".iteration",".draw"),
      variable.name = "pair_name",
      value.name = "draws"
    )

  dt_draws[
    ,
    pair_name_split := pair_name %>%
      stringr::str_remove("x\\[") %>%
      stringr::str_remove("\\]$") %>%
      strsplit("\\,")
  ][,
    `:=`(
      AE = sapply(pair_name_split, function(x) x[1]),
      Drug = sapply(pair_name_split, function(x) x[2])
    )
  ][
    ,
    mcmc_sample_idx := paste(.chain, .iteration, .draw, sep = "_")
  ]

  draws_list <- dt_draws[
    ,
    .(draw_list =
        .(
          .SD[, .(AE, Drug, draws)] %>%
            data.table::dcast(
              AE ~ Drug,
              value.var = "draws",
              fill = NA
            ) %>%
            {
              tmp <- .
              tmp[
                ,
                setdiff(colnames(tmp), "AE"),
                with = FALSE
              ] %>%
                as.matrix() %>%
                magrittr::set_rownames(tmp$AE)
            }
        )
    ),
    by = mcmc_sample_idx
  ]

  return(draws_list$draw_list)

}

