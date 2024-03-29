% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_tree_structure.R
\name{Ytree}
\alias{Ytree}
\title{Title: Create the Outcomes for Each Interior Knots from the Tree}
\usage{
Ytree(Y, treeinf)
}
\arguments{
\item{Y}{The taxa outcome as input}

\item{treeinf}{The information extracted from "phylo" structured tree file}
}
\value{
A set of n * 2 matrices, each of which represent the an interior knot
and its children branches.
}
\description{
The Dirichlet Multinational Tree model requires the
outcomes at each interior nodes. The original taxa outcomes will not be applied directly
into the regression models. The regression model is built
above is for each interior knot and its children. The information from
\code{treeinfo} will be used to build the membership for each interior cluster.
}
\details{
This is adapted from \link{Tao Wang and Hongyu Zhao (2017)}
}
