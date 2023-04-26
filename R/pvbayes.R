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
#' @returns
#' \itemize{
#'   \item lambda_draws - A tibble contain the MCMC samples for each poisson rate parameter.
#'   \item contin_table_long - The converted contingency table in long format tibble.
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

  if (is.null(stan_obj)){

    mod <- .pvBayes$stanmodels[[glue::glue(model,"model",.sep = "_")]]

  } else {

    mod <- stan_obj

  }

  stopifnot(model %in%
              c("zip", "horseshoe", "beta_prime", "regularized_zip", "regularized_horseshoe", "regularized_zip2"))

  I <- nrow(contin_table)
  J <- ncol(contin_table)
  name.c <- colnames(contin_table)
  name.r <- rownames(contin_table)

  dots <- list(...)


  scale_global <- dots$scale_global
  nu_global <- dots$nu_global
  nu_local <- dots$nu_local
  slab_scale <- dots$slab_scale
  slab_df <- dots$slab_df

  if (is.null(scale_global)) scale_global <- 1
  if (is.null(nu_global)) nu_global <- 1
  if (is.null(nu_local)) nu_local <- 1
  if (is.null(slab_scale)) slab_scale <- 1
  if (is.null(slab_df)) slab_df <- 1

  table_E <- contin_table %>%
    {tcrossprod(rowSums(.), colSums(.)) / sum(.)}

  lambda_txt <- matrix(NA, I, J, dimnames = list(name.r, name.c))

  for (i in 1:I) {
    for (j in 1:J) {
      lambda_txt[i, j] <- glue::glue("lambda[{i},{j}]")
    }
  }

  table_long_txt <- lambda_txt %>%
    `colnames<-`(name.c) %>%
    data.table::as.data.table() %>%
    {.[ , AE :=  name.r]} %>%
    data.table::melt(id.vars = c("AE"),
                     measure.vars = name.c[name.c != "AE"],
                     variable.name = "Drug", value.name = "lambda") %>%
    data.table::setkey(Drug, AE)

  table_long_E <- table_E %>%
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
    merge(table_long_txt)

  mod.fit <- mod$sample(
    data = list(
      I = I,
      J = J,
      n = contin_table,
      E = table_E,
      nu_global = nu_global,
      nu_local = nu_local,
      scale_global = scale_global,
      slab_scale = slab_scale,
      slab_df = slab_df
    ),
    seed = stan_seed,
    chains = stan_chains,
    parallel_chains = stan_parallel_chains
  )

  par_vec <- c("lambda", "omega", "kappa", "zi")

  for (k in 1:length(par_vec)){

    temp <- try(
      {mod.fit$draws(format = "df", variables = par_vec[k]) %>%
          data.table::as.data.table()%>%
          {.[,.SD, .SDcols = .[, grepl(paste0("^", par_vec[k]), colnames(.))]]}
      },
      silent=TRUE)

    if (inherits(temp, "try-error")) {
      assign(paste0( par_vec[k],"_draws"), NULL)
    } else{
      assign(paste0( par_vec[k],"_draws"), temp)
    }

  }

  return(
    list(
      lambda_draws = lambda_draws,
      omega_draws = omega_draws,
      kappa_draws = kappa_draws,
      zi_draws = zi_draws,
      contin_table_long = table_long
    )

  )


}

