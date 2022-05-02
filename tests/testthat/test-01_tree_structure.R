test_that("belong function", {
  B <- c(2, 3, 4)
  A <- 3:10

  expect_equal(belong(setA = A, subB = B), c(FALSE, TRUE, TRUE))
})


## Ytree ----------------
test_that("Ytree for outcomes", {
  library(ape)
  library(caper)
  library(tidyverse)
  library(LDTM)

  tree1 <- ape::read.tree(text = "((A, B), (C, D));")
  plot(tree1, edge.width = 2)

  Trinfo <- treeinfo(tree1)
  str(Trinfo)

  Ymat <- Ytree(Y = toy_data1,
                treeinf = Trinfo)
  return(Ymat)
})



