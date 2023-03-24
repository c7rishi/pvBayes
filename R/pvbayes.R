#' Fitting Poisson model with Bayesian methods for pharmacovigilance
#' @param contin_table IxJ contingency table showing pairwise counts of adverse events for I AE (along the rows) and J Drugs (along the columns)
#' @param model Select Bayesian methods
#' \itemize{
#'   \item `horseshoe` - model 1
#'   \item `beta_prime` - model 2
#' }
#' @param stan_obj User defined model in stan with input: `n` counts in each cells, `E` expected value for each cells, and `N`total number of cells IxJ
#' @param stan_seed A seed for the (P)RNG to pass to CmdStan.
#' @param stan_chains The number of Markov chains to run in CmdStan. The default is 2.
#' @param stan_parallel_chains The maximum number of MCMC chains to run in parallel in CmdStan.
#' @returns
#' \itemize{
#'   \item lambda_draws - A tibble contain the MCMC samples for each poisson rate parameter.
#'   \item contin_table_long - The converted contigency table in long format tibble.
#' }
#' @examples
#' library(pvLRT)
#' data(statin64)
#' mod <- pvbayes(contin_table = statin64, model = "horseshoe)
#'
#' #obtain the MCMC samples
#' mod$lambda_draws
#'
#' @export
pvbayes <- function(contin_table,
                    model = "horseshoe",
                    stan_obj = NULL,
                    stan_seed = 123,
                    stan_chains = 2,
                    stan_parallel_chains = 2,
                    ...){

  if (is.null(names(.pvBayes$stanmodels))) {
    msg <- glue::glue(
      "Compiled stan models not found. Please run pvBayes_setup()"
    )
    stop(msg)
  }

  table_E <-
    contin_table %>%
    {tcrossprod(rowSums(.),colSums(.))  / sum(.)} %>%
    `dimnames<-`(
      value = list(rownames(contin_table), colnames(contin_table))
    ) %>%
    dplyr::as_tibble( #dplyr
      rownames = "AE"
    ) %>%
    tidyr::pivot_longer(
      cols = -AE,
      names_to = "Drug",
      values_to = "E"
    )

  suppressMessages(
    table_long <-
      contin_table %>%
      dplyr::as_tibble(
        rownames = "AE"
      ) %>%
      tidyr::pivot_longer(
        cols = -AE,
        names_to = "Drug",
        values_to = "Count"
      )  %>%
      dplyr::left_join(
        y = table_E
      ) %>%
      dplyr::mutate(
        lambda = stringr::str_c("lambda[", as.character(c(1:nrow(table_E))), "]")
      )
  )


  if (is.null(stan_obj)){

    mod <- .pvBayes$stanmodels[[glue::glue(model,"model",.sep = "_")]]

  } else {

    mod <- stan_obj

  }

  mod.fit <- mod$sample(
    data = list(
      N = nrow(table_long),
      n = table_long$Count,
      E = table_long$E
    ),
    seed = stan_seed,
    chains = stan_chains,
    parallel_chains = stan_parallel_chains
  )

  lambda_draws <-
    mod.fit$draws(format = "df") %>%
    dplyr::select(dplyr::starts_with("lambda"))

  return(
    list(
      lambda_draws = lambda_draws,
      contin_table_long = table_long
    )

  )


}

