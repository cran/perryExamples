# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @export
coef.ladlasso <- function(object, zeros = TRUE, tol = .Machine$double.eps^0.5,
                       ...) {
  ## extract coefficients
  coef <- object$coefficients
  ## if requested, omit zero coefficients
  if(!isTRUE(zeros)) {
    if(is.null(dim(coef)))
      coef <- coef[abs(coef) > tol]
    else {
      keep <- apply(abs(coef) > tol, 1, any)
      coef <- coef[keep, , drop = FALSE]
    }
  }
  ## return coefficients
  coef
}

#' @export
coef.lasso <- function(object, zeros = TRUE, ...) {
  ## extract coefficients
  coef <- object$coefficients
  ## if requested, omit zero coefficients
  if(!isTRUE(zeros)) {
    if(is.null(dim(coef)))
      coef <- coef[coef != 0]
    else {
      keep <- apply(coef != 0, 1, any)
      coef <- coef[keep, , drop = FALSE]
    }
  }
  ## return coefficients
  coef
}
