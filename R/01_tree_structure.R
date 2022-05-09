
# 1.1 belong {{{---------
#' Title: Provide each elements from set B belongs to set membership of set A
#' @description `belong` provides whether each element of set B belong to set A
#'
#' @details
#' (TBF)
#'
#' @param setA whether the elements in this setA
#' @param subB the elements takes from subsetB
#'
#' @return A list of logic values indicates the membership;
#' the length of this list is the same as subB
#'
#' @examples
#' library(tidyverse)
#' B <- c(2, 3, 4)
#' A <- 3:10
#' belong(setA = A, subB = B)
#' @export

belong <- function(setA, subB) {

  logic <- unlist(purrr::map(subB, ~is.element(., setA)))

  return(logic = logic)
}
#}}}--------



# 1.2 treeinfo {{{---------
#' Title: Extract Structural Information from a Phylo Tree
#'
#' @description `treeinfo()` provides the basic information from a phylogenic tree
#' without the distances.
#'
#' @details
#' This is adapted from [Tao Wang and Hongyu Zhao (2017)]
#'
#' @param tree A "phylo" structured tree file
#'
#' @return A list object contains necessary information from the `tree`.
#' The output has the following properties:
#'
#' * Nnodes: `integer` the number of nodes in the given "phylo" tree
#' * Ntips: `integer` the number of tips/leaves in the given "phylo" tree
#' * Tredge: `matrix` the matrix for the edges information;
#' each row represents the two ends of the edges
#' * TreeMat: `matrix` the matrix consists of logical values;
#' the binary value of each taxa (by row) belongs to
#' a given leaf or inner node (by col).
#'
#' @export

treeinfo <- function(tree) {
  Nnodes <- tree$Nnode
  Ntips <- length(tree$tip.label)
  Tredge <- tree$edge
  Nodetips <- caper::clade.members.list(tree,
                                        tips = T,
                                        tip.labels = F,
                                        include.nodes = F)
  ## if the Ntips belongs to Nodetips
  ## return TRUE, otherwise FALSE;
  ## then locate at correct position in the matrix
  TreeMat <- purrr::map_dfr(Nodetips, belong, 1:Ntips)

  return(treeinfo = list(Nnodes = Nnodes,
                         Ntips = Ntips,
                         Tredge = Tredge,
                         TreeMat = TreeMat))
}

#}}} --------------




# 1.3 Ytree {{{-------------
#' Title: Create the Outcomes for Each Interior Knots from the Tree
#'
#' @description The Dirichlet Multinational Tree model requires the
#' outcomes at each interior nodes. The original taxa outcomes will not be applied directly
#' into the regression models. The regression model is built
#' above is for each interior knot and its children. The information from
#' `treeinfo` will be used to build the membership for each interior cluster.
#'
#' @details
#' This is adapted from [Tao Wang and Hongyu Zhao (2017)]
#'
#' @param Y The taxa outcome as input
#' @param treeinf The information extracted from "phylo" structured tree file
#'
#' @return A set of n * 2 matrices, each of which represent the an interior knot
#' and its children branches.
#'
#'
#' @export

Ytree <- function(Y, treeinf) {
  Ytree <- list()
  for (j in seq_along(1:treeinf$Nnodes)) {
    index.j <- treeinf$Tredge[, 1] == (treeinf$Ntips + j)
    children.j <- treeinf$Tredge[index.j, 2]
    Ytree[[j]] <- matrix(, nrow(Y), length(children.j))

    for (k in seq_along(1:length(children.j))) {
      Ytree.j.k <- Y[, as.matrix(treeinf$TreeMat)[, children.j[k]],
                     drop = FALSE]
      Ytree[[j]][, k] <- rowSums(Ytree.j.k)
    }
  }
  return(Ytree = Ytree)
}

# }}}-------------



# 1.4 Btree {{{-------------
#' Title: Create the Coefficients for Each Interior Knots Matrix
#'
#' @description The Dirichlet Multinational Tree model requires the
#' outcomes at each interior knots. We need to change a list of coefficients
#' into high dimensional matrix format corresponding to each matrix in `Ytree` outcomes.
#' The length for `Btree` outcomes is the same as `Ytree`; TBF
#'
#' @details
#' This is adapted from [Tao Wang and Hongyu Zhao (2017)]
#'
#' @param Ytree The results from `Ytree()`
#' @param B a list of covariate coefficients
#'
#' @return A set of n * 2 matrices, each of which represent the covariate coefficients
#' at the an interior knot and its children branches.
#'
#' @export
Btree <- function(Ytree, B) {
  Btree <- list()
  a <- 1
  for (j in 1 : length(Ytree)) {
    a.j <- ncol(Ytree[[j]])
    Btree[[j]] <- B[a : (a + a.j - 1), ]
    a <- a + a.j
  }
  return(Btree)
}
# }}}-----------


