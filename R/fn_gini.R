#' Gini Coefficient of a Count Distribution
#'
#' @description
#' Evenness of a count vector: `0` when every element is identical, approaching
#' `1` as the counts concentrate into a few elements. Standard summary for
#' library representation (reads or cells per guide, per construct, per clone).
#'
#' Note that sampling noise alone puts the Gini well above `0` at low coverage,
#' so a raw value is not interpretable on its own. Compare it against the value
#' a perfectly even library of the same size and depth would give, which
#' [fn_gini_poisson()] simulates.
#'
#' @param x Numeric. Counts, one per feature. Zeros are meaningful and are kept;
#'   `NA` values are removed.
#'
#' @return A single numeric between 0 and 1, or `NA_real_` if all counts are
#'   zero.
#'
#' @examples
#' fn_gini(rep(10, 100))          # 0, perfectly even
#' fn_gini(c(rep(0, 99), 1000))   # near 1, all reads on one feature
#'
#' \dontrun{
#' # a Gini is only interpretable against its sampling floor
#' counts <- rpois(1000, 4)
#' fn_gini(counts)
#' fn_gini_poisson(n = 1000, lambda = 4)
#' }
#' @export
fn_gini <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) stop("fn_gini needs at least one non-NA value.")
  if (any(x < 0)) stop("fn_gini is undefined for negative counts.")
  if (sum(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n^2 * mean(x))
}

#' Gini Coefficient Expected from Sampling Alone
#'
#' @description
#' Simulates the Gini coefficient a perfectly even library would show at a given
#' size and depth, i.e. the floor imposed by Poisson sampling. Subtract this
#' from an observed [fn_gini()] to get the unevenness that is actually in the
#' library rather than in the sequencing.
#'
#' @param n Integer. Number of features (guides, constructs, clones).
#' @param lambda Numeric. Mean counts per feature.
#' @param n.sims Integer. Simulation replicates. Default: `200`.
#'
#' @return A single numeric: the mean simulated Gini.
#'
#' @examples
#' set.seed(1)
#' fn_gini_poisson(n = 1000, lambda = 4)    # ~0.27 - the floor at 4x
#' fn_gini_poisson(n = 1000, lambda = 100)  # much lower at 100x
#' @export
fn_gini_poisson <- function(n, lambda, n.sims = 200) {
  if (length(n) != 1 || length(lambda) != 1)
    stop("fn_gini_poisson takes a single n and lambda.")
  mean(replicate(n.sims, fn_gini(stats::rpois(n, lambda))))
}
