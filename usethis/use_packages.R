# license
usethis::use_gpl3_license()

# create package level documentation
usethis::use_package_doc()
# use roxygen
usethis::use_roxygen_md()

# namespace
usethis::use_namespace()

# github remote cmdstanr
usethis::use_dev_package("cmdstanr", type= "Imports", remote = "stan-dev/cmdstanr")

# import
usethis::use_import_from("glue", "glue")

# %>% pipe
usethis::use_pipe()



# Rcpp, Rcpp armadillo

usethis::use_rcpp_armadillo()


# import packages
import_packages <- c(
  "magrittr",
  "glue",
  "stats", "methods", "utils"
)
for (pkg in import_packages) usethis::use_package(pkg)


usethis::use_import_from("cmdstanr", "cmdstan_model")
