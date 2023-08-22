#' @title Simulated Normality (SiNo) Test
#'
#' @description The Simulated Normality (or SiNo) Test is a Lilliefors-like
#' Kernel Density Estimate Test for Normality. It uses the integral of the
#' L1-distance over the reals as the test statistic. The null distribution
#' is estimated using Monte Carlo simulations. The realized sample
#' is first standardized (z-score) and then, the kernel density estimate
#' is obtained based on the standardized sample. This kernel density
#' estimate is then compared with the standard normal density.
#'
#' @param x A vector containing the data
#' @param simulations The number of Monte Carlo simulations used to
#' generate the null distribution when \code{ref_dist = NULL}
#' @param precision The number of equally spaced points at which the density
#' is to be estimated using \code{stats::density}
#' @param def_lim The value of \code{cut} in \code{stats::density}, which
#' dictates how the density values at extremes are estimated (values dropped
#' to approximately 0 at the extremes)
#' @param ref_dist The null distribution, either as an object of class
#' \code{density} (i.e., kernel density estimate) or a cumulative distribution
#' function (i.e., an object of class \code{function})
#' @param sig_level The level of significance
#'
#' @return A named list containing the following components: \cr
#' \itemize{
#'  \item{\code{data_name}}{ - The name of the dataset provided}
#'  \item{\code{data}}{ - The dataset placed under the test}
#'  \item{\code{statistic}}{ - The test statistic}
#'  \item{\code{p_value}}{ - The p-value}
#'  \item{\code{simulations}}{ - The value passed to \code{simulations} when
#'   the function was called}
#'  \item{\code{precision}}{ - The value passed to \code{precision} when
#'   the function was called}
#'  \item{\code{def_lim}}{ - The value passed to \code{def_lim} when
#'   the function was called}
#'  \item{\code{ref_dist}}{ - The null distribution used, either as a
#'  density object or a cumulative distribution function}
#' }
#'
#' @examples
#' sino_test(rnorm(1000))
#' sino_test(rexp(100), simulations = 500)
#'
#' @note If \code{ref_dist = NULL}, the results of the test will be
#' non-deterministic for sample sizes less than 5 or greater than 1000. This is
#' because the null distribution is generated every time using Monte Carlo
#' simulations for such cases. Pre-computed null distributions are available
#' only for samples with size between 5 and 1000 (inclusive). Thus, running the
#' test on the same dataset more than once with \code{ref_dist = NULL} may
#' yield different test statistics (but they are expected to be close).
#'
#' @keywords normality
#'
#' @seealso [sino_null()]
#'
#' @import stats
#' @importFrom cubature cubintegrate
#' @importFrom spatstat.explore CDF
#'
#' @export

sino_test <- function(x, simulations = 1000, precision = 1024, def_lim = 10,
                      ref_dist = NULL, sig_level = 0.05) {

  ## Function for calculating test statistic
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

  ## Function for estimating the null distribution
  sino_null <- function(size, prec, dlim, sims) {
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

  ## Perform the test
  if (is.null(ref_dist)) {
    if ((length(x) >= 5) && (length(x) <= 1000)) {
      null_dist <- SinoTest::null_dists[[length(x) - 4]]
    } else {
      null_dist <- sino_null(length(x), precision, def_lim, simulations)
    }
  } else {
    null_dist <- ref_dist
  }
  test_stat <- sino_stat(x, precision, def_lim)
  if (class(null_dist) == "density") {
    cdf_null <- CDF(null_dist)
  } else if (class(null_dist) == "function") {
    cdf_null <- null_dist
  }
  p_val <- 1 - cdf_null(test_stat)
  if (p_val <= sig_level) {
    interp <- "Reject the null hypothesis"
  } else {
    interp <- "Fail to reject the null hypothesis"
  }
  if (p_val == 0) {
    p_val <- "< 0.00001"
  }
  message(paste("\nSimulated Normality Test\n",
                "Null Hypothesis: The population distribution is a ",
                "normal distirbution.\n",
                "Data: ", deparse(substitute(x)), "\n",
                "Test Statistic: ", test_stat, "\n",
                "p-value: ", p_val, "\n",
                "Result: ", interp,
                " at ", sig_level, " level of significance\n", sep = ""))
  results <- list(data_name = deparse(substitute(x)), data = x,
                  statistic = test_stat, p_value = p_val,
                  simulations = simulations, precision = precision,
                  def_lim = def_lim,
                  ref_dist = ifelse(is.null(ref_dist), cdf_null,
                                    deparse(substitute(ref_dist))))

}
