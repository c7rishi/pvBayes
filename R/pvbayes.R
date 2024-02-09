#' Fitting Poisson model with Bayesian methods for pharmacovigilance
#' @importFrom data.table `:=`
#' @importFrom data.table `%like%`
#' @param contin_table IxJ contingency table showing pairwise counts of adverse events for I AE (along the rows) and J Drugs (along the columns). Colnames and rownames are needed.
#' @param model Options of model:
#' \itemize{
#'   \item `poisson` - poisson model with horseshoe prior
#'   \item `poisson_indep` - poisson model under independent assumption with horseshoe prior
#'   \item `zip` - zero-inflated poisson model with horseshoe prior
#'   \item `zip_indep` - zero-inflated poisson model under independent assumption with horseshoe prior
#'   \item `beta_prime` - beta prime model
#' }
#' @param starting option of HMC start point
#' \itemize{
#'   \item `random` - random starting point chosen by cmdstanr
#'   \item `LRT` - likelihood ratio estimator by pvLRT
#'   \item `MAP`- maximum posterior estimator
#' }
#' @param stan_obj User defined model in stan with input: `n` counts in each cells, `E` expected value for each cells, and `N` total number of cells IxJ
#' @param stan_seed A seed for the (P)RNG to pass to CmdStan.
#' @param stan_chains The number of Markov chains to run in CmdStan. The default is 2.
#' @param stan_iter_sampling An integer of total sampling number.
#' @param retry An integer for maximum repeated number of Stan sampling process.
#' @param ... additional parameters, currently is not used.
#' @returns
#' \itemize{
#'   \item E - Contingency table of expected counts.
#'   \item draws - A list contains the MCMC samples for each parameters.
#'   \item MCMC_convergence - Information of Stan sampler.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#' }
#' @export
pvbayes <- function(contin_table,
                    model = "zip_horseshoe",
                    stan_obj = NULL,
                    starting = "random",
                    stan_seed = 123,
                    stan_chains = 4,
                    stan_iter_sampling = 1000,
                    stan_cov_rate = 0.5,
                    retry = 10,
                    return_stan = FALSE,
                    ...){

  #find compiled model
  if ( is.null(names(.pvBayes$stanmodels)) ) {

    msg <- glue::glue(
      "Compiled stan models not found. Please run pvBayes_setup()"
    )
    stop(msg)

  }

  # allow input stan object
  if ( is.null(stan_obj) ) {

    #check model name exist
    if( !model %in% names(.pvBayes$stanmodels) ) {

      msg <- glue::glue(
        "Model name does not exist. Please check model_name() to fit avaible models."
      )
      stop(msg)

    }
    mod <- .pvBayes$stanmodels[[model]]

  } else {

    mod <- stan_obj

  }

  I <- nrow(contin_table)
  J <- ncol(contin_table)
  name.c <- colnames(contin_table)
  name.r <- rownames(contin_table)

  table_E <- contin_table %>%
    {tcrossprod(rowSums(.), colSums(.)) / sum(.)}

  stan_data <- list(
    I = I,
    J = J,
    n = contin_table,
    E = table_E
  )

  starting_list <- NULL

  if (starting != "random") {

    if (starting == "LRT") {

      lambda_start <- ifelse(contin_table/table_E > 1, contin_table/table_E, 1)
      omega_start <- pvLRT::pvlrt(contin_table, nsim = 2) %>%
        pvLRT::extract_zi_probability()

    } else if ( starting == "MAP" ) {

      temp_optim <-
        mod$optimize(
          data = stan_data,
          seed = stan_seed,
          refresh = 0
        )

      lambda_start <- temp_optim$mle("lambda") %>%
        matrix(
          nrow = I,
          ncol = J
        )
      omega_start <- temp_optim$mle("omega") %>% unname()

    }

    starting_list <- lapply(
      1:stan_chains,
      function(i)
        list(
          log_lambda = log(lambda_start),
          omega = omega_start
        )
    )


  }

  stan_retry <- TRUE
  n_retry <- 0

  while( stan_retry ) {

    all_stan_in <- list(
      data = stan_data,
      seed = stan_seed,
      chains = stan_chains,
      refresh = 500,
      init = starting_list,
      iter_sampling = stan_iter_sampling
    )

    mod.fit <- do.call(
      mod$sample,
      all_stan_in
    )

    div_rate <- mod.fit$diagnostic_summary(
      quiet = TRUE
    )$num_divergent %>%
      {sum(.)/(stan_iter_sampling * stan_chains)}

    stan_retry <- {div_rate > stan_cov_rate} & {n_retry <= retry}
    n_retry <- n_retry + 1
    stan_seed <- stan_seed + 1

  }

  convergence_msg <- ifelse(stan_retry, "Error", "Completed")

  # browser()

  par_vec <- c(
    "lambda",
    "omega",
    "kappa",
    "zi",
    "n_pred",
    "rho_Drug",
    "rho_AE"
  )

  draws_list <-  1:length(par_vec) %>%
    setNames(
      par_vec
    ) %>%
    lapply(
      FUN = function(k){
        temp <- tryCatch(
          {
            mod.fit$draws(format = "draws_matrix",
                          variables = par_vec[k]) %>%
              posterior::as_draws_rvars() %>%
              .[[par_vec[k]]]
          },
          error = function(e){

            NULL

          }
        )

        if ( par_vec[k] %in% c("lambda",
                               "lambda_indep",
                               "lambda_resid",
                               "n_pred")
             & !is.null(temp)) {

          dimnames(temp) <- list(name.r, name.c)

        }

        return(temp)
      }
    )

  out_list <- list(
    model = model,
    data = contin_table,
    E = table_E,
    draws = draws_list,
    MCMC_convergence = c(glue::glue("Divergence rate: {div_rate*100}%"),
                         convergence_msg)
  )

  class(out_list) <- "pvbayes"

  return(
    out_list
  )


}

