% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simMGarch.R
\docType{class}
\name{simMGarch-class}
\alias{simMGarch-class}
\title{An S4 class for a nonstationary CCC model.}
\description{
A specification class to create an object of a simulated piecewise constant conditional correlation (CCC) model
denoted by \eqn{r_t = (r_{1, t}, \ldots, r_{n, t})^T}, \eqn{t=1, \ldots, n} with
\eqn{r_{i, t}= \sqrt{h_{i, t}}\epsilon_{i, t}} where \eqn{h_{i, t}= \omega_i(t) + \sum_{j=1}^p \alpha_{i, j}(t)r_{i, t-j}^2 + \sum_{k=1}^q \beta_{i, k}(t)h_{i, t-k}}.
In this package, we assume a piecewise constant CCC with \eqn{p=q=1}.
}
\section{Slots}{

\describe{
\item{\code{y}}{The \eqn{n \times d} time series.}

\item{\code{cor_errors}}{The \eqn{n \times d} matrix of the errors.}

\item{\code{h}}{The \eqn{n \times d} matrix of the time-varying variances.}

\item{\code{n}}{Size of the time series.}

\item{\code{d}}{The number of variables (assets).}

\item{\code{r}}{A sparsity parameter to conrol the impact of changepoint across the series.}

\item{\code{multp}}{A parameter to control the covariance of errors.}

\item{\code{changepoints}}{The vector with the location of the changepoints.}

\item{\code{pw}}{A logical parameter to allow for changepoints in the error covariance matrix.}

\item{\code{a0}}{The vector of the parameters a0 in the individual GARCH processes denoted by \eqn{\omega_i(t)} in the above formula.}

\item{\code{a1}}{The vector of the parameters a1 in the individual GARCH processes denoted by \eqn{\alpha_i(t)} in the above formula.}

\item{\code{b1}}{The vector of the parameters b1 in the individual GARCH processes denoted by \eqn{\beta_i(t)} in the above formula.}

\item{\code{BurnIn}}{The size of the burn-in sample. Note that this only applies at the first simulated segment. Default is 50.}
}}

\examples{
pw.CCC.obj <- new("simMGarch")
pw.CCC.obj <- pc_cccsim(pw.CCC.obj)
par(mfrow=c(2,2))
ts.plot(pw.CCC.obj@y[1,]);ts.plot(pw.CCC.obj@y[2,])
ts.plot(pw.CCC.obj@h[1,]);ts.plot(pw.CCC.obj@h[1,])
}
\references{
Cho, H. and Korkas, K.K., 2022. High-dimensional GARCH process segmentation with an application to Value-at-Risk. Econometrics and Statistics, 23, pp.187-203.
}
