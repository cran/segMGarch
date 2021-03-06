# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @keywords internal
stl_sort <- function(x) {
    .Call(`_segMGarch_stl_sort`, x)
}

#' @keywords internal
rcpp_rev <- function(x) {
    .Call(`_segMGarch_rcpp_rev`, x)
}

#' @keywords internal
func_coef <- function(z, scale) {
    .Call(`_segMGarch_func_coef`, z, scale)
}

#' @keywords internal
func_input <- function(coef, sgn, sq, diag) {
    .Call(`_segMGarch_func_input`, coef, sgn, sq, diag)
}

#' @keywords internal
func_dc <- function(z, gamma) {
    .Call(`_segMGarch_func_dc`, z, gamma)
}

#' @keywords internal
func_cusum <- function(z) {
    .Call(`_segMGarch_func_cusum`, z)
}

#' @keywords internal
func_cusum_vec <- function(z) {
    .Call(`_segMGarch_func_cusum_vec`, z)
}

#' @keywords internal
func_dc_by <- function(z, gamma, dmby, dtby) {
    .Call(`_segMGarch_func_dc_by`, z, gamma, dmby, dtby)
}

#' @keywords internal
func_mvt_ar <- function(ar, res) {
    .Call(`_segMGarch_func_mvt_ar`, ar, res)
}

func_density <- function(z, c) {
    .Call(`_segMGarch_func_density`, z, c)
}

func_density_by <- function(z, dmby, dtby, c) {
    .Call(`_segMGarch_func_density_by`, z, dmby, dtby, c)
}

func_input_off <- function(coef, sgn, log, sq) {
    .Call(`_segMGarch_func_input_off`, coef, sgn, log, sq)
}

func_input_on_boot <- function(coef, log, sq, avg) {
    .Call(`_segMGarch_func_input_on_boot`, coef, log, sq, avg)
}

func_input_off_boot <- function(coef, sgn, log, sq, avg) {
    .Call(`_segMGarch_func_input_off_boot`, coef, sgn, log, sq, avg)
}

