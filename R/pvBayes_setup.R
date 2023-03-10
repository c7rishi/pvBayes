.pvBayes <- new.env(parent = emptyenv())

#' Set up stan compiled objects for pvBayes package
#' @param ... additional arguments passed to cmdstanr::cmdstan_model()
#' @export
pvBayes_setup <- function(...) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    msg <- glue::glue(
      "{cmdstanr} not found. pvBayes uses {cmdstanr} interface \\
      to stan for MCMC sampling from the Bayesian models. Please \\
      install {cmdstanr} before using {pvBayes}. For windows systems \\
      we recommend using the wsl version of cmdstan if available.",
      .open = "<<", .close = ">>"
    )
    stop(msg)
  }

  msg <- glue::glue(
    "\nInitializing pvBayes with cmdstanr::cmdstan_model().\\
     It may take a few minutes to compile the stan programs, please wait...
     \n \\
     The compiled programs will be stored in cache for faster \\
     compilation next time onwards...\n"
  )
  message(msg)

  .pvBayes$stanmodels <- new.env(parent = emptyenv())

  cache_loc <- tools::R_user_dir(package = "pvBayes", which = "cache")
  if (!dir.exists(cache_loc)) dir.create(cache_loc)

  stan_source_files <- c(
    list.files("./stan/", pattern = "\\.stan$", full.names = TRUE),
    list.files("./inst/stan/", pattern = "\\.stan$", full.names = TRUE)
  )
  stan_models <- stan_source_files %>%
    strsplit("\\/") %>%
    sapply(tail, 1) %>%
    gsub("\\.stan", "", .)

  for (ii in seq_along(stan_models)) {
    # browser()
    mod <- stan_models[ii]
    stan_file <- stan_source_files[ii]
    exe_file <- glue::glue("{cache_loc}/{mod}_stan.exe")
    key <- list(...)
    msg <- glue::glue("Compiling model '{mod}'...\n")
    message(msg)
    .pvBayes$stanmodels[[mod]] <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      exe_file = exe_file,
      compile = TRUE,
      ...
    )
  }

  msg <- glue::glue("\nDone!\n")
  message(msg)
}


.onLoad <- function(lib, pkg) {

  # user_permission <- utils::askYesNo("Set up stan compiled models with default options?")

  tmp <- tryCatch(
    pvBayes_setup(),
    error = function(e) e
  )
  if (is (tmp, "error")) {
    msg <- glue(
      "Running pvBayes_setup() with default setting failed...\n\\
      {pvBayes} uses the cmdstanr interface to stan to run \\
      the MCMC samplers. Please run pvBayes_setup() to \\
      activate the stan models. See ?pvBayes_setup() \\
      for more details.",
      .open = "<<", .close = ">>"
    )
    message(msg)
  }


}
