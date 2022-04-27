
# 2.1 Loglikehood {{{------------
#' Title: Likelihood function for Dirichlet and Multinomial models
#'
#' @param Y
#' @param X
#' @param b
#' @param model
#'
#' @return
#' @export
#'
#'
Loglik <- function(Y, X, b, model) {
  # Compute the log likelihood, constant part discarded
  # The likelihood is scaled. Be careful when computing AIC BIC
  X <- as.matrix(X)
  g <- exp(X %*% t(b))	# n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  if (model == "dirmult") {
    res <- 	sum(lgamma(gs) - lgamma(ys + gs) +
                  rowSums(lgamma(Y + g) - lgamma(g)))
  }
  if (model == "mult") {
    res <- 	sum(rowSums(Y * log(g)) - ys * log(gs))
  }
  if (model == "dir") {
    res <- 	sum(lgamma(gs) + rowSums((g-1) * log(Y) - lgamma(g)))
  }

  return(res / nrow(X))
}
# }}}--------------



# 2.2 Score {{{--------------
#' Title Score function for Dirichlet and Multinomial models
#'
#' @param Y
#' @param X
#' @param b
#' @param model
#'
#' @return
#' @export
#'
Score <- function(Y, X, b, model) {
  # Compute the Score function at b
  X <- as.matrix(X)
  S <- 0
  g <- exp(X %*% t(b))	# n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)

  # ?digamma the first derivative of a Gamma function
  ## digamma(x) = d/dx{ ln \Gamma(x)} = \Gamma'(x) / \Gamma(x)

  if (model == "dirmult") {
    S <-  t((digamma(gs) - digamma(ys + gs) +
               digamma(Y + g) - digamma(g)) * g) %*% X
  }

  if (model == "mult") {S <- t((Y / g - ys / gs) * g) %*% X}

  if (model == "dir") {S <-  t((digamma(gs) - digamma(g) + log(Y)) * g) %*% X}

  return(S / nrow(X))
}
# }}}---------------






# 2.3 Hessian {{{-----------------
#' Title Hessian matrix for Dirichlet and Multinomial models
#'
#' @param Y
#' @param X
#' @param b
#' @param model
#'
#' @return
#' @export
#'
Hessian <- function(Y, X, b, model){
  #	Compute the diagonal of the hessian matrix at b

  X <- as.matrix(X)
  H <- 0
  g <- exp(X %*% t(b))	# n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  if (model == "dirmult") {
    H <- t((trigamma(gs) - trigamma(ys + gs) +
              trigamma(Y + g) - trigamma(g)) * g^2 +
             (digamma(gs) - digamma(ys+gs) +
                digamma(Y + g) - digamma(g)) * g) %*% X^2
  }
  if (model == "mult") {
    H <- t((-Y / g^2 + ys / gs^2) * g^2  +
             (Y / g - ys / gs) * g) %*% X^2
  }
  if (model == "dir") {
    H <- t((trigamma(gs) - trigamma(g)) * g^2  +
             (digamma(gs)  - digamma(g) + log(Y)) * g) %*% X^2
  }
  #	Divided by the sample size
  H / nrow(X)
}
# }}}---------------















