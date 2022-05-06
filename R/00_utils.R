.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to my package; this is a package adapted from other people's papers")
}


# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.devtools <- list(
#     devtools.path = "~/R-dev",
#     devtools.install.args = "",
#     devtools.name = "Randy",
#     devtools.desc.author = "Randy Jin <xin.2.jin@cuanschutz.edu> [aut, cre]",
#     devtools.desc.license = "What license is it under?",
#     devtools.desc.suggests = NULL,
#     devtools.desc = list()
#   )
#   toset <- !(names(op.devtools) %in% names(op))
#   if(any(toset)) options(op.devtools[toset])
#
#   invisible()
# }



# 0.1 norm L2 {{{------------------
#' Title L2 or other norms
#'
#' @param v
#'
#' @return
#' @export
#'
norm2 <- function(v) {
  sqrt(sum(v^2))
}

## Add other type of norms
# }}}---------------



# 0.2 soft {{{-------------------
#' Title soft activation function
#' @details this is the equation
#'
#' @param M
#' @param v
#'
#' @return
#' @export
#'
soft_fun <- function(M, v) {
  soft <- M
  ## number of cols
  p <- dim(M)[2]
  for (j in 1 : p) {
    M.j <- M[, j]
    temp.j <- abs(M.j) - v
    temp.j[temp.j < 0] <- 0
    soft[, j] <- sign(M.j) * temp.j
  }
  return(soft)
}

# }}}-------------------



# 0.3
