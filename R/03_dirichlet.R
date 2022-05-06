# 3.1 QL_fun {{{------------------


#' Title  Gradient Descent Method for penalized likelihood
#' this is gradient descent methods for
#' penalized likelihood
#'
#' @param Ytree
#' @param X
#' @param W
#' @param model
#' @param B1
#' @param grad
#' @param alpha
#' @param lambda
#' @param L
#'
#' @return
#' @export
#'
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

## }}}---------------------




## 3.2 lambda_fun --------------------
#' Title Getting betas use the max(beta, 1)
#'
#' @param grad
#' @param L
#' @param alpha
#' @param lambda
#'
#' @return
#' @export
#'
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


## 3.3 lambda_raw_fun {{{---------------
#' Title Grid methods to find which lambda minimize the loss
#'
#' @param grad
#' @param L
#' @param alpha
#' @param lambda.raw
#' @param fac1
#' @param fac2
#'
#' @return
#' @export
#'
#' @examples
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
## }}}-------------------------


## 3.4 Tree path {{{--------------------------

#' Title to get the tree path and alpha for penalty likelihood
#'
#' @param Y
#' @param X
#' @param treeinfo
#' @param alpha
#' @param cutoff
#' @param model
#' @param err.conv
#' @param iter.max
#' @param L.init
#' @param lambda.max
#' @param lambda.min
#'
#' @return
#' @export
#'

T_group_path <- function(Y,
                         X,
                         treeinfo,
                         alpha = 0.5,
                         cutoff = 0.8,
                         model = "dirmult",
                         err.conv = 1e-3,
                         iter.max = 30,
                         L.init = NULL,
                         lambda.max = NULL,
                         lambda.min = NULL) {
  ## trying with training data
  # View(XY.train)
  # Y <- data1
  # X <- as.matrix(tree_X)
  # treeinfo <- info
  # alpha = 0.5
  # cutoff = 0.8
  # model = "dirmult"
  # err.conv = 1e-3
  # iter.max = 30
  # L.init = NULL
  # lambda.max = NULL
  # lambda.min = NULL

  # L.init = NULL
  # lambda.max = NULL
  # lambda.min = NULL
  X <- as.matrix(X)
  n <- dim(Y)[1] ## 75
  # dim(Y) ## 75 * 28
  # dim(X) ## 75 * 102

  p <- dim(X)[2] ## 28
  p1 <- p - 1
  Ytree <- Y_tree(Y, treeinfo)


  ## number of rows of Y_tree???
  ## J = 54
  J <- dim(treeinfo$tree.edge)[1]

  ## starting values for betas
  B0 <- matrix(0, J, p)

  ## ??????
  ## likelihood starting with certain value ---------
  if (is.null(L.init)) L.init <- p1 * J / n

  L <- L.init

  ## lambda starting with certain value ----------
  if (is.null(lambda.max)) {
    # model <- "dirmult"
    Sc <- -Score_Dirichlet_tree(Ytree, X, B0, model)

    # View(Score_Dirichlet_tree)
    # View(Score)
    # View(Sc)
    # dim(Sc) # 54 * 102

    grad <- Sc[, -1]
    # dim(grad) # remove the first col
    lambda.max <- lambda_raw_fun(grad, L, alpha)

    print(paste("lambda.max:", lambda.max))
  }

  lambda <- lambda.max

  iter <- 0
  nonzeroprop <- 0
  output.path <- list()

  ## path iteration -------------------
  ## trying the initate parameters
  # alpha = 0.5
  # cutoff = 0.8
  # err.conv = 1e-3
  # iter.max = 30
  while (nonzeroprop < cutoff) {

    ### initiation -------------------
    L <- L.init
    B.new <- B.old <-
      W.new <- W.old <- B0
    B1.old <- W1.old <- B0[, -1]
    # View(B0)
    loop <- 0
    err.val <- err.val1 <- err.val2 <- 1
    theta.old <- 1
    # View(B0)
    # View(B1.old)

    (try.err0 <-
        ## will produce error
        try(f.old <-
              f_fun(Ytree, X, B.old,
                    model, alpha, lambda)))
    # View(f_fun)

    ### whether there is trying error -----------
    if (class(try.err0) == "try-error" |
        is.nan(try.err0)) {

      print(paste("iter ", iter, ": Loglik (B.old) error!", sep = ""))

      ## redo betas
      (B <- B0 <- B.new)
      (iter <- iter + 1)

      output.path[[iter]] <-
        list(B = B,
             lambda = lambda,
             loop = loop,
             err.val = err.val)

      nonzeroprop <- sum(B[, -1] != 0) / p1 / J
      nonzeroprop

      (lambda <- lambda.max * 0.95^iter)

      #### whether lambda provided -------------
      if (!is.null(lambda.min)) {
        if (lambda < lambda.min) break
      }
    } else {
      ### if there is no error -------------
      while ((loop < iter.max) &
             (err.val > err.conv)) {
        Sc <- -Score_Dirichlet_tree(Ytree, X, W.old, model)
        grad <- Sc[, -1]
        flag <- T
        loop.inn <- 0

        while (flag & loop.inn < 50) {
          (V <- W1.old - grad / L)

          B1.new <- V
          B1.new <- soft_fun(V, lambda * (1 - alpha) / L)
          if (alpha != 0) {
            for (j in 1 : p1) {
              B.j <- B1.new[, j]
              norm2.j <- norm2(B.j)
              thres.j <- norm2.j - lambda * alpha / L
              if (thres.j <= 0) {
                B1.new[, j] <- 0
              } else {
                B1.new[, j] <- thres.j / norm2.j * B.j
              }
            }
          }
          B.new[, -1] <- B1.new

          (try.err1 <-
              try(f.new <- f_fun(Ytree, X, B.new, model, alpha, lambda)))
          try.err2 <- try.err1

          if (class(try.err1) == "try-error" | is.nan(try.err1)) break

          (try.err2 <-
             try(QL <- QL_fun(Ytree, X, W.old,
                              model, B1.new,
                              grad, alpha,
                              lambda, L)))

          # Thu Apr 14 13:54:29 2022 ------------------------------
          # cat("dim of W", dim(W.old), "\n",
          #     "dim of B", dim(B1.new), "\n",
          #     "dim of grad", dim(grad), "\n")

          if (class(try.err2) == "try-error" | is.nan(try.err2)) break

          if (QL >= f.new) {
            flag <- F
          } else {
            L <- L * 2
            loop.inn <- loop.inn + 1
          }
        }

        if (class(try.err1) == "try-error" |
            is.nan(try.err1) |
            class(try.err2) == "try-error" |
            is.nan(try.err2)) {
          print(paste("iter ", iter,
                      ": Loglik (B.new) error!", sep = ""))

          break

        } else {

          if (QL < f.new)
            print(paste("iter ", iter, ": descent_QL", sep = ""))
          if (f.new > f.old)
            print(paste("iter ", iter, ": ascent_f", sep = ""))
        }

        theta.new <- 2 / (loop + 3)
        c.temp <- (1 - theta.old) / theta.old * theta.new
        #theta.new <- (1 + sqrt(1 + 4 * theta.old^2)) / 2
        #c.temp <- (theta.old - 1) / theta.new

        if (f.new <= f.old) {
          W1.new <- B1.new + c.temp * (B1.new - B1.old)
          W.new[, -1] <- W1.new
        } else {
          B.new <- B.old
          W1.new <- B1.new + c.temp * (B1.old - B1.new)
          W.new[, -1] <- W1.new
        }

        #err.val1 <- max(abs(B.old - B.new))
        #err.val <- err.val1^2

        #err.val2 <- (f.old - f.new) / f.old
        err.val2 <- (f.old - f.new)
        err.val <- err.val2
        ##print(err.val)
        if (err.val > err.conv) {
          W.old <- W.new
          W1.old <- W.old[, -1]

          theta.old <- theta.new
          f.old <- f.new

          B.old <- B.new
          B1.old <- B.old[, -1]

          loop <- loop + 1
        }
      }

      (B <- B0 <- B.new)
      (iter <- iter + 1)
      output.path[[iter]] <- list(B = B,
                                  lambda = lambda,
                                  loop = loop,
                                  err.val = err.val)
      nonzeroprop <- sum(B[, -1] != 0) / p1 / J

      (lambda <- lambda.max * 0.95^iter)
      if (!is.null(lambda.min)) {
        if (lambda < lambda.min) break
      }
      #print(c(loop, err.val, iter))
    }
  }
  # View(output.path)
  return(list(output.path = output.path,
              alpha = alpha))
}

## }}}-------------------------



