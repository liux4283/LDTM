
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


# 2.4 Loglik_Dirichlet_tree{{{----------------
#' Title Log likelihood for Dirichlet Tree
#'
#' @param Ytree
#' @param X
#' @param B
#' @param model
#'
#' @return
#' @export
#'
Loglik_Dirichlet_tree <- function(Ytree, X, B, model) {

  ## this is the original B not the tree-wised
  Btree <- B_tree(Ytree, B)
  N.nodes <- length(Ytree)
  res <- 0
  for (j in 1 : N.nodes) {
    res <- res + Loglik(Ytree[[j]], X, Btree[[j]], model)
  }
  res
}
## }}}-------------------




# 2.5 Score_Dirichlet_tree {{{------------------
#' Title
#'
#' @param Ytree
#' @param X
#' @param B
#' @param model
#'
#' @return
#' @export
#'
Score_Dirichlet_tree <- function(Ytree, X, B, model) {
  Btree <- B_tree(Ytree, B)
  N.nodes <- length(Ytree)
  S <- NULL
  for (j in 1 : N.nodes) {
    S <- rbind(S, Score(Ytree[[j]], X, Btree[[j]], model))
  }
  S
}

## }}}---------------------




# 2.6 Penalized_Loglike {{{----------------------
## thinking about change the name
## this is the penalized likelihood

#' Title Penalized log-likelihood for Dirichlet Tree Multinomial Regression
#'
#' @param Ytree
#' @param X
#' @param B
#' @param model
#' @param alpha
#' @param lambda
#'
#' @return
#' @export
#'
f_fun <- function(Ytree,
                  X, B,
                  model,
                  alpha,
                  lambda) {

  loglik <- Loglik_Dirichlet_tree(Ytree, X, B, model)

  B1 <- B[, -1]
  pty <- sum(abs(B1)) * lambda * (1 - alpha)

  ## alpha is the gamma in the paper
  if (alpha != 0) {
    pty <- pty + sum(apply(B1, 2, norm2)) *
      lambda * alpha
  }

  f <- -loglik + pty
  return(f)
}
# }}}-----------------




# 2.7 Prediction {{{------------------
#' Title Predicted values for Outcomes
#'
#' @param Y
#' @param X
#' @param B
#' @param treeinfo
#'
#' @return
#' @export
#'
pred_Y <- function(Y, X, B, treeinfo) {
  X <- as.matrix(X)
  Ytree <- Y_tree(Y, treeinfo)
  Btree <- B_tree(Ytree, B)

  n <- dim(Y)[1]
  N.tips <- treeinfo$N.tips
  tree.edge <- treeinfo$tree.edge
  inci <- treeinfo$inci

  predY <- Y
  for (j in 1 : N.tips) {
    p.j <- rep(1, n)
    index.j <- setdiff(which(inci[j, ]), j)
    index.j <- c(index.j, j)
    for (k in 1 : (length(index.j) - 1)) {
      index.j.k <- tree.edge[, 1] == index.j[k]
      children.j.k <- tree.edge[index.j.k, 2]

      p.j.k <- exp(X %*% t(Btree[[index.j[k] - N.tips]]))
      p.j.k <- p.j.k[, children.j.k == index.j[k + 1]] / rowSums(p.j.k)
      p.j <- p.j * p.j.k
    }
    predY[, j] <- p.j * rowSums(Y)
  }
  return(predY)
}

## }}}--------------------






