
load("data/simu_data.RData")
source("12_Hongzhe_Dismo.R")
source("13_Dirichlet_tree_Tao.R")
source("14_functions_group_Tao.R")
source("15_T_group_Tao.R")


treeinfo1 <- tree_extract(tree)


treeinfo2 <- treeinfo1
treeinfo2$N.nodes <- 1
treeinfo2$tree.edge <- cbind(treeinfo1$N.tips + 1, 
                             1 : treeinfo1$N.tips)
treeinfo2$inci <- cbind(diag(1, treeinfo1$N.tips), 1) == 1
## information about the trees --------------------
## tips 28
## nodes 27
## edges 54
## edges + nodes 55
# View(treeinfo1$inci)

# ytree <- Y_tree(XY.train[[1]]$Y, treeinfo1)
# View(ytree)

set.seed(1)
type <- "group"
alpha <- c(0)
nfolds <- 5
CV <- F

## this is the tree group including  tips + nodes
if (type == "T_group") treeinfo <- treeinfo1
## this is the just the tips
if (type == "group") treeinfo <- treeinfo2

# repl.times <- 100
output <- est.Aic <- 
  est.Bic <- list()


min.i.Aic <- min.i.Bic <- 
  err.pred.Aic <- err.pred.Bic <-
  matrix(, repl.times, length(alpha))

# View(XY.train[[1]]$X)
# Thu Apr 14 10:17:16 2022 ------------------------------
## did not use cross validation
for (repl in seq_along(1:repl.times)) {
  
  XY.train.repl <- XY.train[[repl]]
  XY.test.repl <- XY.test[[repl]]
  output[[repl]] <- 
    est.Aic[[repl]] <- 
    est.Bic[[repl]] <- list()
  
  # View(XY.train.repl$X)
  # View(XY.train)
  
  for (j in seq_along(1:length(alpha))) {
    # repl <- 1
    ## cross-validation for tree group
    output[[repl]][[j]] <- 
      cv_T_group_fun(XY.train.repl$Y, 
                     XY.train.repl$X, 
                     treeinfo, 
                     alpha = alpha[j], 
                     nfolds = nfolds, 
                     CV = CV)
  }
  

  IC.repl <- IC_fun(output[[repl]], 
                    XY.train[[repl]]$Y, 
                    XY.train[[repl]]$X, 
                    treeinfo)
  min.i.Aic.repl <- IC.repl$min.i.Aic
  min.i.Bic.repl <- IC.repl$min.i.Bic
  # View(IC.repl)
  
  
  for (j in 1 : length(alpha)) {

    ## the output path for each repl
    output.path.repl.j <- output[[repl]][[j]]$output.path
    # View(alpha)
    # View(output.path.repl.j)

    est.Aic[[repl]][[j]] <- output.path.repl.j[[min.i.Aic.repl[j]]]$B
    est.Bic[[repl]][[j]] <- output.path.repl.j[[min.i.Bic.repl[j]]]$B

    pred.testY.Aic.repl.j <- pred_Y(XY.test[[repl]]$Y, 
                                    XY.test[[repl]]$X, 
                                    est.Aic[[repl]][[j]], 
                                    treeinfo)
    err.pred.Aic[repl, j] <- sum((XY.test.repl$Y - pred.testY.Aic.repl.j)^2 / 
                                   pred.testY.Aic.repl.j / 
                                   rowSums(XY.test.repl$Y))

    pred.testY.Bic.repl.j <- pred_Y(XY.test[[repl]]$Y,
                                    XY.test[[repl]]$X, 
                                    est.Bic[[repl]][[j]],
                                    treeinfo)
    err.pred.Bic[repl, j] <- sum((XY.test.repl$Y - pred.testY.Bic.repl.j)^2 /
                                   pred.testY.Bic.repl.j / 
                                   rowSums(XY.test.repl$Y))
  }
}


View(est.Aic)
View(output)
