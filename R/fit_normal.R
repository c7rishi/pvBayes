#' Fit a N(mu, sigma) model to univariate data y
#' @param y input data
#' @param method method of posterior approximation. Defaults to MCMC sampling.
#' @param ... passed to appropriate stan method (sample, variational, or optimize)
#' @export
fit_normal <- function(y, method = "sample", ...) {

  if (is.null(names(.pvBayes$stanmodels))) {
    msg <- glue::glue(
      "Compiled stan models not found. Please run pvBayes_setup()"
    )
    stop(msg)
  }

  stan_data <- list(y = c(y), N = length(y))

  stan_mod <- .pvBayes$stanmodels$normal_sample

  stan_fn <-  if (method == "variational") {
    stan_mod$variational
  } else if (method == "optimize") {
    stan_mod$optimize
  } else {
    stan_mod$sample
  }

  fit <- tryCatch(
    stan_fn(data = stan_data, ...),
    error = function(e) e
  )

  if (is(fit, "error")) browser()

  fit
}
