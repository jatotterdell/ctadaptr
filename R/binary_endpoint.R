#' Calculate density of beta-binomial distribution
#'
#' @param x The value at which to evaluate the density
#' @param n The sample size
#' @param a First parameter
#' @param b Second parameter
#'
#' @return Value of beta-binomial(n,a,b) evaluated at x
#'
#' @examples
#' dbetabinom(5, 10, 2, 3)
#'
#' @export
dbetabinom <- function(x, n, a = 1, b = 1){
  if(!(all(c(a, b) > 0))) stop("a and b must be > 0")
  if(any(n < 1)) stop("n must be > 0")
  if(any(x < 0)) stop("x must be >= 0")

  num <- lgamma(a + b) + lgamma(n + 1) + lgamma(x + a) + lgamma(n - x + b)
  den <- lgamma(a) + lgamma(b) + lgamma(x + 1) + lgamma(n - x + 1) + lgamma(n + a + b)
  prob <- exp(num - den)
  prob
}

#' Draw random variates from beta-binomial distribution
#'
#' @import stats
#'
#' @param n The number of random values to sample
#' @param m The sample size
#' @param a First parameter
#' @param b Second parameter
#'
#' @examples
#' rbetabinom(2, 10, 2, 3)
#'
#' @export
rbetabinom <- function(n, m, a = 1, b = 1) {
  if(!(all(c(a, b) > 0))) stop("a and b must be > 0")
  if(!(all(n > 0))) stop("n must be > 0")

  stats::rbinom(n, m, stats::rbeta(n, a, b))
}


#' Beta inequaility
#'
#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using numerical integration
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param ... other arguments passed to integrate/quadgk function
#'
#' @return The value of the integral
#'
#' @examples
#' beta_ineq(5, 5, 3, 7)
#'
#' @export
beta_ineq <- function(a, b, c, d, delta = 0, ...) {

  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")

  integrand <- function(x) { stats::dbeta(x, a, b)*stats::pbeta(x - delta, c, d) }
  tryCatch(
    integrate(integrand, delta, 1, ...)$value,
    error = function(err) NA)
}

#' Beta inequaility - normal approximation
#'
#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Normal approximation.
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#'
#' @return The value of the integral
#'
#' @examples
#' beta_ineq_approx(5, 5, 3, 7)
#'
#' @export
beta_ineq_approx <- function(a, b, c, d, delta = 0) {
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")

  m1 <- a / (a + b)
  v1 <- a*b / ( (a + b)^2 * (a + b + 1))
  m2 <- c / (c + d)
  v2 <- c*d / ( (c + d)^2 * (c + d + 1))
  z <- (m1 - m2 - delta) / sqrt(v1 + v2)
  return(stats::pnorm(z))
}

#' Beta inequality - Monte Carlo
#'
#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Monte Carlo method.
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param sims The number of Monte Carlo variates to generate for estimation
#'
#' @return The value of the integral
#'
#' @examples
#' beta_ineq_sim(5, 5, 3, 7)
#'
#' @export
beta_ineq_sim <- function(a, b, c, d, delta = 0, sims = 10000) {
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")

  X <- stats::rbeta(sims, a, b)
  Y <- stats::rbeta(sims, c, d)
  mean(X > Y + delta)
}

#' Calculate the predicted probability of success
#'
#' @param a First parameter of first beta random variable
#' @param b Second parameter of first beta random variable
#' @param c First paramter of second beta random variable
#' @param d Second parameter of second beta random variable
#' @param m1 Sample size to predict for first beta random variable
#' @param m2 Sample size to predict for second beta random variable
#' @param k_ppos The posterior probability cut-point to be assessed
#' @param post_method Method to use for calculating posterior probability.
#' One of `exact`, `approx` (default) or `sim`
#'
#' @return The predicted probability of success
#'
#' @import data.table
#' @export
calc_ppos <- function(a, b, c, d, m1, m2, k_ppos, post_method = "approx") {
  require(data.table)
  if(!(all(c(a, b, c, d, m1, m2) > 0))) stop("a, b, c, d, m1, m2 must be > 0")
  if(k_ppos < 0 | k_ppos > 1) stop("k_ppos must be in [0, 1]")

  calc_post <- switch(post_method,
                      "exact" = beta_ineq,
                      "approx" = beta_ineq_approx,
                      "sim" = beta_ineq_sim)

  y1pred <- rbetabinom(10000, m1, a, b)
  y2pred <- rbetabinom(10000, m2, c, d)
  ypred <- data.table(y1pred = y1pred, y2pred = y2pred)[, .N, keyby = list(y1pred, y2pred)]
  ypred[, `:=`(P = Vectorize(calc_post)(a + y1pred,
                                        b + m1 - y1pred,
                                        c + y2pred,
                                        d + m2 - y2pred))]
  ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
}
