
## norm L2 ------------------
norm2 <- function(v) {
  sqrt(sum(v^2))
}


## penalty likelihood ???-------------------
## this is just lasso? 
# Thu Apr 21 12:36:47 2022 ------------------------------
  ## equation (4)
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

## basically this is the equation (8) in the paper
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

## getting betas--------------------
lambda_fun <- 
  function(grad, ## the score function
           L, ## C in equation 8
           alpha, ## alpha
           lambda){

  V <- -grad / L
  p <- dim(V)[2]

  W <- V
  W <- soft_fun(V, lambda * (1 - alpha) / L)
  
  
  if (alpha != 0) {
    for (j in 1 : p) {
      W.j <- W[, j]
      norm2.j <- norm2(W.j)
      thres.j <- norm2.j - lambda * alpha / L
      if (thres.j <= 0) {
        W[, j] <- 0
      } else {
        W[, j] <- thres.j / norm2.j * W.j
      }
    }
  }
  return(W)
}


## grid methods
lambda_raw_fun <- 
  function(grad, 
           L, ## likelihood
           alpha, 
           lambda.raw = 2, 
           fac1 = 1.1, 
           fac2 = 0.96) {
    
  doit1 <- T
  while(doit1) {
    lambda.raw <- lambda.raw * fac1
    W <- lambda_fun(grad, L, alpha, lambda.raw)
    nonzeros <- sum(W != 0)
    if (nonzeros == 0) doit1 <- F
  }
  
  lambda.max <- lambda.raw
  
  doit2 <- T
  while(doit2) {
    lambda.max <- lambda.max * fac2
    W <- lambda_fun(grad, L, alpha, lambda.max)
    nonzeros <- sum(W != 0)
    
    if (nonzeros > 0) doit2 <- F
  }

  lambda.max <- lambda.max / fac2
  return(lambda.max)
}
