#' Fitting Poisson model with Bayesian methods for pharmacovigilance
#' @importFrom data.table `:=`
#' @importFrom data.table `%like%`
#' @param contin_table IxJ contingency table showing pairwise counts of adverse events for I AE (along the rows) and J Drugs (along the columns)
#' @param model Select Bayesian methods
#' \itemize{
#'   \item `horseshoe` - model 1
#'   \item `beta_prime` - model 2
#' }
#' @param stan_obj User defined model in stan with input: `n` counts in each cells, `E` expected value for each cells, and `N` total number of cells IxJ
#' @param stan_seed A seed for the (P)RNG to pass to CmdStan.
#' @param stan_chains The number of Markov chains to run in CmdStan. The default is 2.
#' @param stan_parallel_chains The maximum number of MCMC chains to run in parallel in CmdStan.
#' @param ... additional parameters
#' @returns
#' \itemize{
#'   \item lambda_draws - A tibble contain the MCMC samples for each poisson rate parameter.
#'   \item contin_table_long - The converted contingency table in long format tibble.
#' }
#' @examples
#' \dontrun{
#' library(pvLRT)
#' data(statin46)
#' mod <- pvbayes(contin_table = statin46, model = "horseshoe")
#'
#' #obtain the MCMC samples
#' mod$lambda_draws
#' }
#' @export
pvbayes <- function(contin_table,
                    model = "horseshoe",
                    stan_obj = NULL,
                    stan_seed = 123,
                    stan_chains = 4,
                    stan_parallel_chains = 4,
                    ...){

  #find compiled model
  if (is.null(names(.pvBayes$stanmodels))) {
    msg <- glue::glue(
      "Compiled stan models not found. Please run pvBayes_setup()"
    )
    stop(msg)
  }

  # allow input stan object
  if (is.null(stan_obj)){

    #check model name exist
    if( !model %in% names(.pvBayes$stanmodels) ){
      msg <- glue::glue(
        "Model name not exists. Please check ?pvbayes() to find avaible models."
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

  dots <- list(...)
  additional_var <- list(
    scale_global = 1,
    scale_local = 1,
    nu_global = 1,
    nu_local = 1,
    slab_scale = 1,
    slab_df = 1,
    c_alpha = 1,
    c_beta = 1,
    gamma = 10,
    b = 1
  )

  for (name in names(additional_var)) {
    if (!is.null(dots[[name]])) {
      additional_var[[name]] <- dots[[name]]
    }
  }

  table_E <- contin_table %>%
    {tcrossprod(rowSums(.), colSums(.)) / sum(.)}

  stan_data = list(
    I = I,
    J = J,
    n = contin_table,
    E = table_E
  )

  mod.fit <- mod$sample(
    data = c(stan_data, additional_var),
    seed = stan_seed,
    chains = stan_chains,
    parallel_chains = stan_parallel_chains,
    refresh = 500
  )

  par_vec <- c("lambda",
               "omega",
               "kappa",
               "zi",
               "n_pred")

  draws_list <- list()

  for (k in 1:length(par_vec)){

    temp <- tryCatch(
      {mod.fit$draws(format = "draws_matrix", variables = par_vec[k]) %>%
          posterior::as_draws_rvars() %>%
          .[[par_vec[k]]]
      },
      error = function(e){
        NULL
      }
    )

    if (is.null(temp)) {next}

    if ( length(dim(temp)) == 1) {
      names(temp) <-  name.c
    } else{
      dimnames(temp) <- list(name.r, name.c)
    }

    draws_list[[par_vec[k]]] <- temp

  }

  return(
    list(
      E = table_E,
      draws = draws_list
    )
  )


}

