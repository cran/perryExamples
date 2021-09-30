# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## add default column names to matrix
addColnames <- function(x) {
  # 'x' needs to be a matrix
  if(is.null(colnames(x))) colnames(x) <- paste("x", seq_len(ncol(x)), sep="")
  x
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
  if(!check || all(is.na(match(c("Intercept", "(Intercept)"), colnames(x))))) {
    cbind("(Intercept)"=rep.int(1, nrow(x)), x)
  } else x
}

## backtransform regression coefficients to original scale
backtransform <- function(beta, muX, sigmaX, muY, intercept = TRUE) {
  bt <- function(b) {
    b <- b / sigmaX
    a <- muY - sum(b * muX)  # intercept
    c("(Intercept)"=a, b)
  }
  if(is.null(dim(beta))) bt(beta) else apply(beta, 2, bt)
}

## remove intercept column from design matrix
removeIntercept <- function(x, pos) {
  if(missing(pos)) {
    pos <- match(c("Intercept", "(Intercept)"), colnames(x), nomatch = 0)
    if(any(pos > 0)) x[, -pos, drop=FALSE] else x
  } else x[, -pos, drop=FALSE]
}
