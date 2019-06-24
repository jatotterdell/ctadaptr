#' Probability that each cell is best
#'
#' @param mat A matrix of cell mean posterior draws.
#'
#' @return The probability that each column of \code{mat}
#'   is the maximum.
#'
#' @export
prob_best <- function(mat) {
  prop.table(table(factor(max.col(mat), levels = 1:ncol(mat))))
}

#' Probability that one cell is superior
#'
#' @param mat A matrix of cell mean posterior draws.
#' @param col The columnm being assessed for superiority.
#' @param eps The superiority margin.
#'
#' @return The probability that \code{col} is superior to
#'   to the other columns in \code{mat} by margin \code{eps}.
#'
#' @export
prob_superior <- function(mat, col, eps) {
  mean(apply(mat, 1, function(x) all(x[col] - x[-col] > eps)))
}

#' Probability that all cells are equivalent
#'
#' @param mat A matrix of cell mean posterior draws.
#' @param eps The equivalence margin.
#'
#' @return The probability that all cells in \code{mat}
#'   are equivalent within the bound \code{eps}.
prob_all_equivalent <- function(mat, eps) {
  a <- rep(1/ncol(mat), ncol(mat))
  mean(apply(sweep(mat, 1, drop(mat %*% a), `-`), 1,
             function(x) all(abs(x) <= eps)))
}

#' Probability that all cells are noninferior.
#'
#' @param mat A matrix of cell mean posterior draws.
#' @param col The reference column for noninferiority.
#' @param eps The noninferiority bound.
#'
#' @return The probability that all cells in \code{mat}
#'   are noninferior to the reference cell \code{col} by
#'   a non-inferiority margin of \code{eps}.
prob_all_noninferior <- function(mat, col, eps) {
  mean(apply(mat, 1, function(x) all(x[-col] - x[col] > -eps)))
}
