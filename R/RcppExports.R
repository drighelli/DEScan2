# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' rcpparma_get_disjoint_max_win
#' Computes the disjoint max_win matrix.
#' @param z0 a matrix.
#' @param sigwin sigwin.
#' @param zthresh zthresh.
#' @param nmax nmax.
#' @param verbose verbose.
#' @return a matrix of three columns (bin_idx, win_idx, z_val) idxs in C style.
#' @keywords internal
rcpparma_get_disjoint_max_win <- function(z0, sigwin = 10L, zthresh = 10, nmax = 9999999L, verbose = TRUE) {
    .Call('_DEScan2_rcpparma_get_disjoint_max_win', PACKAGE = 'DEScan2', z0, sigwin, zthresh, nmax, verbose)
}

