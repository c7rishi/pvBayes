#' Bisection function for solving function
#' @importFrom data.table :=
#' @importFrom data.table %like%
#' @noRd
bisection <- function(f, a, b, n = 1000, tol = 1e-5, ...) {

  c <- (a + b) / 2
  if (
    n == 1 |
    ((abs(f(a) - f(c)) < tol) & (abs(f(b) - f(c)) < tol)) |
    ((b - a) < tol)
  )  return(c)

  if (f(a) * f(c) <= 0) {
    # cat("a=", a, "b=", c, "\n")
    bisection(f, a = a, b = c, n = n - 1)
  } else if (f(c) * f(b) <= 0) {
    # cat("a=", c, "b=", b, "\n")
    bisection(f, a = c, b = b, n = n - 1)
  } else {
    stop("interval does not have opposite signs")
  }

}

