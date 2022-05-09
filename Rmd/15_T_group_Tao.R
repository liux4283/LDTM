## Tree group ------------------------------------------------------------------
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


##
cv_T_group_fun <- 
  function(Y, X, 
           treeinfo, 
           alpha = 0.5, 
           cutoff = 0.8, 
           model = "dirmult",
           err.conv = 1e-3, 
           iter.max = 30,
           L.init = NULL, 
           lambda.max = NULL, 
           lambda.min = NULL,
           nfolds = 10, 
           CV = T, 
           cv.ind = NULL) {

  output <- T_group_path(Y, X, 
                         treeinfo, 
                         alpha, 
                         cutoff,
                         model,
                         err.conv, 
                         iter.max,
                         L.init, 
                         lambda.max, 
                         lambda.min)
  
  # cross validation starts here -------------
  
  if (CV == F) {
    return(output)
  } 
  
  ## cross validation
  else {
      output.path <- output$output.path
      path.length <- length(output.path)
      lambda.max <- output.path[[1]]$lambda
      lambda.min <- output.path[[path.length]]$lambda

    E <- matrix(NA, nrow = nfolds, ncol = path.length)
    if (is.null(cv.ind)) cv.ind <- kfold.dismo(1 : dim(X)[1], nfolds)

    for (i in 1 : nfolds) {

      cat("Starting CV fold #", i, sep = "", "\n")
      X1 <- X[cv.ind != i, , drop = FALSE]
      Y1 <- Y[cv.ind != i, , drop = FALSE]
      X2 <- X[cv.ind == i, , drop = FALSE]
      Y2 <- Y[cv.ind == i, , drop = FALSE]

      output.i <- T_group_path(Y1, X1, treeinfo,
                               alpha, cutoff = 1, 
                               model, err.conv, iter.max,
                               L.init, lambda.max, 
                               lambda.min)

      output.path.i <- output.i$output.path
      path.length.i <- length(output.path.i)
      for (j in 1 : path.length.i) {

        pred.Y2.i.j <- pred_Y(Y2, X2, output.path.i[[j]]$B, treeinfo)
        E[i, j] <- sum((Y2 - pred.Y2.i.j)^2 / pred.Y2.i.j / rowSums(Y2))
      }
    }

    cve <- apply(E, 2, mean, na.rm = T)
    cvse <- apply(E, 2, sd, na.rm = T) / sqrt(nfolds)
    
    min.i <- which.min(cve)[1]
    min.1se <- which(cve <= cve[min.i] + cvse[min.i])[1]

    output.cv <- list(cve = cve, 
                      cvse = cvse, 
                      min.i = min.i, 
                      min.1se = min.1se, 
                      E = E, 
                      output = output)
    return(output.cv)
  }
}


## Getting the information criteria --------------------------
IC_fun <- function(output, Y, X,
                   treeinfo, model = "dirmult") {

  print("lasso and group_lasso_only!")

  n <- dim(X)[1]
  Ytree <- Y_tree(Y, treeinfo)

  Aic.list <- Bic.list <-
    Df.list <- list()
  Aic <- Bic <- Df.Aic <-
    Df.Bic <- min.i.Aic <- min.i.Bic <- 
    rep(NA, length(output))

  for (j in 1 : length(output)) {

    output.path.j <- output[[j]]$output.path
    alpha.j <- output[[j]]$alpha
    Aic.j <- Bic.j <- Df.j <- NULL

    for (k in 1 : length(output.path.j)) {

      B.j.k <- output.path.j[[k]]$B
      loglik.j.k <- Loglik_Dirichlet_tree(Ytree, X, B.j.k, model = model)
      dev.j.k <- -2 * loglik.j.k * n

      if (alpha.j == 1) {
        Df.j.k <- sum(colSums(B.j.k[, -1] != 0) != 0) * dim(B.j.k[, -1])[1]
      } else {
        Df.j.k <- sum(B.j.k[, -1] != 0)
      }

      Aic.j.k <- dev.j.k + Df.j.k * 2
      Bic.j.k <- dev.j.k + Df.j.k * log(n)

      Df.j <- c(Df.j, Df.j.k)
      Aic.j <- c(Aic.j, Aic.j.k)
      Bic.j <- c(Bic.j, Bic.j.k)
    }

    Df.list[[j]] <- Df.j
    Aic.list[[j]] <- Aic.j
    Bic.list[[j]] <- Bic.j

    min.i.Aic[j] <- which.min(Aic.j)[1]
    Df.Aic[j] <- Df.j[min.i.Aic[j]]
    Aic[j] <- Aic.j[min.i.Aic[j]]

    min.i.Bic[j] <- which.min(Bic.j)[1]
    Df.Bic[j] <- Df.j[min.i.Bic[j]]
    Bic[j] <- Bic.j[min.i.Bic[j]]
  }

  return(list(Aic = Aic, 
              Bic = Bic, 
              Df.Aic = Df.Aic, 
              Df.Bic = Df.Bic,
              min.i.Aic = min.i.Aic, 
              min.i.Bic = min.i.Bic,
              Aic.list = Aic.list, 
              Bic.list = Bic.list, 
              Df.list = Df.list))
}
