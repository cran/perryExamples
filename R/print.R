# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## print method shared by lasso-type methods
print_lasso_type <- function(x, zeros = FALSE, ...) {
  # print coefficients
  cat("\nCoefficients:\n")
  print(coef(x, zeros=zeros), ...)
  # return object invisibly
  invisible(x)
}

#' @export
print.lasso <- print_lasso_type

#' @export
print.ladlasso <- print_lasso_type

#' @export
print.ridge <- function(x, ...) {
  # print coefficients
  cat("\nCoefficients:\n")
  print(coef(x), ...)
  # return object invisibly
  invisible(x)
}
