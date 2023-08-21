#' @title Simulated Normality (SiNo) Test
#'
#' @description This generates the null distirbution for the SiNo Test.
#'
#' @param size The sample size
#' @param prec The number of equally spaced points at which the density
#' is to be estimated using \code{stats::density}
#' @param dlim The value of \code{cut} in \code{stats::density}, which
#' dictates how the density values at extremes are estimated (values dropped
#' to approximately 0 at the extremes)
#' @param sims The number of Monte Carlo simulations used to
#' generate the null distribution
#'
#' @return An object of class \code{density}
#'
#' @examples
#' sino_null(550)
#'
#' @seealso [sino_test()]
#'
#' @import stats
#' @importFrom cubature cubintegrate
#'
#' @export

sino_null <- function(size, prec = 1024, dlim = 10, sims = 1000) {

  ## Calculate test statistic
  sino_stat <- function(dat, prec, dlim) {
    dat <- (dat - mean(dat)) / sd(dat)
    kde <- density(dat, n = prec, cut = dlim)
    kde_pdf <- approxfun(kde, yleft = 0, yright = 0)
    abs_diff <- function(v) {
      return(abs(kde_pdf(v) - dnorm(v)))
    }
    return(cubintegrate(abs_diff, lower = -Inf, upper = Inf,
                        method = "pcubature")$integral)
  }

  ## Estimate the null distribution
  i <- 1
  test_stats <- rep(0, sims)
  while (i <= sims) {
    set.seed(i)
    rand_samp <- rnorm(size)
    test_stats[i] <- sino_stat(rand_samp, prec, dlim)
    i <- i + 1
  }
  kde <- density(test_stats, n = prec, cut = dlim)
  return(kde)

}
