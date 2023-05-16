#' Obtain model names
#'
#' @export
model_name <- function(){

  stan_source_files <- system.file("stan/", package = "pvBayes") %>%
    list.files(pattern = "\\.stan$", full.names = TRUE)

  if (length(stan_source_files) == 0) {
    stan_source_files <- c(
      list.files("./inst/stan/", pattern = "\\.stan$", full.names = TRUE)
    )
  }

  model_names <- stan_source_files %>%
    basename() %>%
    {sub("\\.stan$", "", .)}

  return(model_names)

}
