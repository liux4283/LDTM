% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_likelihood.R
\name{Loglik}
\alias{Loglik}
\title{Title: Likelihood function for Dirichlet and Multinomial models}
\usage{
Loglik(Y, X, b, model)
}
\arguments{
\item{Y}{\code{matrix} of count outcomes}

\item{X}{\code{matrix} of covariates}

\item{b}{\code{matrix} of dimensions p x q; where p = number of covariates
& q = number of branches.
Note: When q=1, we are conditioning the likelihood on that specific branch}

\item{model}{\code{character} type of model to use for the Log Likelihood. Options are
(Dirichlet Multinomial = "dirmult", Multinomial = "mult", or
Dirichlet = "dir")}
}
\value{
The value of the Log-likelihood for the desired model.
}
\description{
Title: Likelihood function for Dirichlet and Multinomial models
}
\details{
This is adapted from \link{Tao Wang and Hongyu Zhao (2017)}
This function is used to extract the log likelihood from Dirichlet & Multinomial models
The log-likelihood value of a regression model is a way to measure the goodness of fit for a model.
}
