#' Fitting Poisson model with Bayesian methods for pharmacovigilance
#' @importFrom data.table :=
#' @importFrom data.table %like%
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
#' data(statin46)
#' mod <- pvbayes(contin_table = statin46, model = "horseshoe")
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

  name.c <- colnames(contin_table)
  name.r <- rownames(contin_table)

  table_long_E <- contin_table %>%
    {tcrossprod(rowSums(.),colSums(.))  / sum(.)} %>%
    `colnames<-`(name.c) %>%
    data.table::as.data.table() %>%
    {.[ , AE :=  name.r]} %>%
    data.table::melt(id.vars = c("AE"),
         measure.vars = name.c[name.c != "AE"],
         variable.name = "Drug", value.name = "E") %>%
    data.table::setkey(Drug, AE)

  table_long <- contin_table %>%
    data.table::as.data.table() %>%
    {.[ , AE :=  name.r]}  %>%
    data.table::melt(id.vars = c("AE"),
         measure.vars = name.c[name.c != "AE"],
         variable.name = "Drug", value.name = "Count") %>%
    data.table::setkey(Drug, AE) %>%
    merge(table_long_E) %>%
    {.[ , lambda :=  paste0("lambda[", 1:nrow(.), "]")]}

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
    data.table::as.data.table() %>%
    {.[,.SD, .SDcols = .[, grepl("^lambda", colnames(.))]]}

  return(
    list(
      lambda_draws = lambda_draws,
      contin_table_long = table_long
    )

  )


}

