test_that("belong function", {
  B <- c(2, 3, 4)
  A <- 3:10

  expect_equal(belong(setA = A, subB = B), c(FALSE, TRUE, TRUE))
})



test_that("Ytree for outcomes", {
  library(ape)
  library(caper)
  library(tidyverse)

  tree1 <- ape::read.tree(text = "((A, B), (C, D));")
  plot(tree1, edge.width = 2)

  ## import dataset
  # toy_data0 <- here::here("data", "toy_data.csv") %>%
  #   read.csv(row.names = 1) %>%
  #   mutate(ID = as.factor(ID)) %>%
  #   dplyr::select(-intercept, -eta, -random, -fixed)
  #
  # toy_data1 <- toy_data0 %>%
  #   dplyr::select(ID, visit, species, y) %>%
  #   pivot_wider(names_from = species,
  #               values_from = y) %>%
  #   unite("id_time", c(ID, visit), sep = "_") %>%
  #   column_to_rownames("id_time")

  Trinfo <- treeinfo(tree1)
  str(Trinfo)

  ## change the data into cluster
  ## each matrix contain the interior knot
  ## as well as its children
  Ymat <- Ytree(Y = data1,
                treeinf = Trinfo)
  Ymat
})



