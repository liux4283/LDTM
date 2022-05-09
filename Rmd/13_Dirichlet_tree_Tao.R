library(caper)
?clade.members.list

## tree_extract ---------------------------------------------
tree_extract <- function(tree) {

  N.nodes <- tree$Nnode
  N.tips <- length(tree$tip.label)
  tree.edge <- tree$edge

  node.tips1 <- 
    clade.members.list(tree, 
                       tips = T, 
                       tip.labels = F, 
                       include.nodes = F)
  inci <- matrix(, N.tips, length(node.tips1))
  
  for (j in 1 : N.tips) {
    inci[j, ] <- sapply(node.tips1, 
                        function(set) is.element(j, set))
  }

  treeinfo <- list()
  treeinfo$N.nodes <- N.nodes
  treeinfo$N.tips <- N.tips
  treeinfo$tree.edge <- tree.edge
  treeinfo$inci <- inci
  return(treeinfo)
}

## the tree matrix dimension
# dim(treeinfo$inci) ## 28 * 29
# View(treeinfo1$inci) ## 28 * 29
# View(treeinfo2$inci) ## 28 * 55

## Y_tree -----------------------------------------------------
## just one person's
## just one time
## with all the species and knobs


# Ynew = 
# # Y R R R 
# # R Y R R
# # R R Y R
# # R R R Y

Y_tree <- function(Y, treeinfo) {
  ## Y is taxa outcome
  # Y <- data1
  # treeinfo = info
  n <- dim(Y)[1]
  N.nodes <- treeinfo$N.nodes
  N.tips <- treeinfo$N.tips
  tree.edge <- treeinfo$tree.edge
  inci <- treeinfo$inci

  Ytree <- list()
  for (j in 1 : N.nodes) {
    index.j <- tree.edge[, 1] == (N.tips + j)
    children.j <- tree.edge[index.j, 2]
    Ytree[[j]] <- matrix(, n, length(children.j))
    
    for (k in 1 : length(children.j)) {
      Ytree.j.k <- Y[, inci[, children.j[k]], drop = FALSE]
      Ytree[[j]][, k] <- rowSums(Ytree.j.k)
    }
  }
  # Ytree
  return(Ytree)
}

## B_tree ----------------------------------------------------------------------
## this the final betas, but in the position corresponding 
## to the 
B_tree <- function(Ytree, B) {

  Btree <- list()
  a <- 1
  for (j in 1 : length(Ytree)) {
    a.j <- dim(Ytree[[j]])[2]
    Btree[[j]] <- B[a : (a + a.j - 1), ]
    a <- a + a.j
  }
  return(Btree)
}
# View(B_tree)


Loglik_Dirichlet_tree <- function(Ytree, X, B, model) {

  # Thu Apr 14 12:00:45 2022 ------------------------------
  ## this is the original B not the tree-wised
  Btree <- B_tree(Ytree, B)
  N.nodes <- length(Ytree)
  res <- 0
  for (j in 1 : N.nodes) {
    res <- res + Loglik(Ytree[[j]], X, Btree[[j]], model)
  }
  res
}
# View(Loglik)
# View(lgamma)


Score_Dirichlet_tree <- function(Ytree, X, B, model) {

  Btree <- B_tree(Ytree, B)
  N.nodes <- length(Ytree)
  S <- NULL
  for (j in 1 : N.nodes) {
    S <- rbind(S, Score(Ytree[[j]], X, Btree[[j]], model))
  }
  S
}
# View(Score)

# (y - X %*% bhat)^2 - plambda

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
