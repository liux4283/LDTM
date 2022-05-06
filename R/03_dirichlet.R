# 3.1 Gradient Descent Method for

## this is gradient descent methods for
## penalized likelihood
QL_fun <- function(Ytree,
                   X,
                   W, ## extra? the starting B
                   model,
                   B1,
                   grad,
                   alpha,
                   lambda,
                   L){

  W1 <- W[, -1]
  # Mon Apr 11 08:15:53 2022 ------------------------------
  ## what is grad???
  ## equation (5)
  QL1 <- -Loglik_Dirichlet_tree(Ytree, X, W, model) +
    sum(diag(t(B1 - W1) %*% grad))

  QL2 <- sum(abs(B1)) * lambda * (1 - alpha)

  if (alpha != 0) {
    QL2 <- QL2 + sum(apply(B1, 2, norm2)) * lambda * alpha
  }

  QL3 <- sum((B1 - W1)^2) * L / 2

  QL <- QL1 + QL2 + QL3
  return(QL)
}

