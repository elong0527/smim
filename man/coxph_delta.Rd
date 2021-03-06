% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxph_detail.R, R/coxph_detail_old.R
\name{coxph_delta}
\alias{coxph_delta}
\title{Delta adjusted Cox Imputation Model}
\usage{
coxph_delta(fit, pattern, delta)

coxph_delta(fit, pattern, delta)
}
\arguments{
\item{fit}{a Cox model object, i.e., the result of \code{coxph}.}

\item{pattern}{A integer vector of pattern indicator for each subject
* 1  = observed event time
* 2  = censoring at random
* 3+ = censoring not at random}

\item{delta}{A numeric vector of adjusted delta in the same length of pattern.
* delta = 1 is for observed event time and censoring at random.
* delta > 1 is for censoring time with hazard is larger than the hazard in CAR situation.
* delta < 1 is for censoring time with hazard is smaller than the hazard in CAR situation.}
}
\description{
Delta adjusted Cox Imputation Model

Delta adjusted Cox Imputation Model
}
