#' Fitting Poisson model with Bayesian methods for pharmacovigilance
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
    {.[,.SD, .SDcols = .[, grepl("^lambda|^omega", colnames(.))]]}

  return(
    list(
      lambda_draws = lambda_draws,
      contin_table_long = table_long
    )

  )


}

