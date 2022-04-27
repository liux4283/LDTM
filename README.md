
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDTM

<!-- badges: start -->
<!-- badges: end -->

The goal of LDTM is to …

## Installation

You can install the development version of LDTM from
[GitHub](https://github.com/Goodgolden/LDTM) with:

``` r
# install.packages("devtools")
devtools::install_github("Goodgolden/LDTM")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(LDTM)
#> Loading required package: tidyverse
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.6     ✓ dplyr   1.0.8
#> ✓ tidyr   1.2.0     ✓ stringr 1.4.0
#> ✓ readr   2.1.2     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
## basic example code
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->

## Introduction

### OTU operational taxonomic units

microbial genomic sequencess clustered by sequence similarity

partition sequences into discrete groups instead of traditional
taxonomic units.

the most abundant sequence in an OTU is the representative sequence

representative sequences from all the OTUs are used to construct a
phylogenetic tree among all the OTUs

Microbial community information == OTUs + counts + phylogenetic
relationship + taxonomy

### Human gut microbiome study in University of Pennsylvania

the effect of diet on gut microbiome composition

-   gut microbiome data
-   nutrient intake data

### distance based analysis

-   identify a few gut microbiome associated nutrients

-   unable to provide information on how dietary nutrient affect
    bacterial taxa

## Goal

identify both the key nutrients as well as the taxa the nutrient affect

Chen and Li (2013) adopted a regression-based approach

OTU abundance data as multivariate count responses, and nutrient as
covaraite

## Model

### Multinomial logistic regression

![
\\begin{split}
& Y = (Y_1, ..., Y_K)^T\\\\
& y = (y_1, ..., y_K)^T\\\\
& \\sum\_{k=1}^K Y_k = \\sum\_{k=1}^K y_k\\\\
& p = (p_1, ..., p_K)^T,\\ \\sum\_{k=1} ^K p_k =1\\\\
& f_M(y;\\ p) = \\frac {\\Gamma ({\\sum\_{k=1}^K y_k + 1})} 
{\\prod\_{k=1}^K \\Gamma ({y_k + 1})} \\prod\_{k=1}^Kp_k^{y_k}\\\\
\\end{split}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bsplit%7D%0A%26%20Y%20%3D%20%28Y_1%2C%20...%2C%20Y_K%29%5ET%5C%5C%0A%26%20y%20%3D%20%28y_1%2C%20...%2C%20y_K%29%5ET%5C%5C%0A%26%20%5Csum_%7Bk%3D1%7D%5EK%20Y_k%20%3D%20%5Csum_%7Bk%3D1%7D%5EK%20y_k%5C%5C%0A%26%20p%20%3D%20%28p_1%2C%20...%2C%20p_K%29%5ET%2C%5C%20%5Csum_%7Bk%3D1%7D%20%5EK%20p_k%20%3D1%5C%5C%0A%26%20f_M%28y%3B%5C%20p%29%20%3D%20%5Cfrac%20%7B%5CGamma%20%28%7B%5Csum_%7Bk%3D1%7D%5EK%20y_k%20%2B%201%7D%29%7D%20%0A%7B%5Cprod_%7Bk%3D1%7D%5EK%20%5CGamma%20%28%7By_k%20%2B%201%7D%29%7D%20%5Cprod_%7Bk%3D1%7D%5EKp_k%5E%7By_k%7D%5C%5C%0A%5Cend%7Bsplit%7D%0A "
\begin{split}
& Y = (Y_1, ..., Y_K)^T\\
& y = (y_1, ..., y_K)^T\\
& \sum_{k=1}^K Y_k = \sum_{k=1}^K y_k\\
& p = (p_1, ..., p_K)^T,\ \sum_{k=1} ^K p_k =1\\
& f_M(y;\ p) = \frac {\Gamma ({\sum_{k=1}^K y_k + 1})} 
{\prod_{k=1}^K \Gamma ({y_k + 1})} \prod_{k=1}^Kp_k^{y_k}\\
\end{split}
")

The link function is a multinomial-Poisson transformation

**might need to use Poissonization to simulate the data**

#### Overdisperson

-   the multinomial distribution is not appropriate

-   ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
    is a random variable with some prior distribution

-   ![\\Phi^d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi%5Ed "\Phi^d")
    is a (d-1) dimensional simplex,

the support of p is
![\\Phi^{K-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi%5E%7BK-1%7D "\Phi^{K-1}")

-   Dirichlet distribution density at
    ![u = (u_1, ..., u_K) \\in \\Phi^{K-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u%20%3D%20%28u_1%2C%20...%2C%20u_K%29%20%5Cin%20%5CPhi%5E%7BK-1%7D "u = (u_1, ..., u_K) \in \Phi^{K-1}")

-   Dirichlet is a conjugate prior to multinomial distribution

-   posterior is Dirichlet multionmial distribution, aka the Dirichlet
    compound multinomial distribution.

![
\\begin{split}
& a = (a_1, ..., a_K)^T,\\ a_i \> 0\\\\
& f_D (u;a) = \\frac {\\Gamma ({\\sum\_{k=1}^K \\alpha_k})} 
{\\prod\_{k=1}^K \\Gamma ({\\sum\_{k=1}^K \\alpha_k})} \\prod\_{k=1}^K u_k^{\\alpha_k - 1}\\\\
& f\_{DM}(y; \\alpha) = \\
\\int\_{u \\in \\Phi^{K-1}} f_M(y; u) f_D(u; \\alpha) \\\\
& = \\frac {\\Gamma ({\\sum\_{k=1}^K y_k + 1}) \\Gamma ({\\sum\_{k=1}^K \\alpha_k})} 
{\\Gamma ({\\sum\_{k=1}^K y_k} + {\\sum\_{k=1}^K \\alpha_k})} 
\\prod\_{k=1}^K \\frac {\\Gamma (y_k + \\alpha_k)} 
{\\Gamma ({y_k + 1}) \\Gamma ({\\alpha_k})}
\\end{split} 
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bsplit%7D%0A%26%20a%20%3D%20%28a_1%2C%20...%2C%20a_K%29%5ET%2C%5C%20a_i%20%3E%200%5C%5C%0A%26%20f_D%20%28u%3Ba%29%20%3D%20%5Cfrac%20%7B%5CGamma%20%28%7B%5Csum_%7Bk%3D1%7D%5EK%20%5Calpha_k%7D%29%7D%20%0A%7B%5Cprod_%7Bk%3D1%7D%5EK%20%5CGamma%20%28%7B%5Csum_%7Bk%3D1%7D%5EK%20%5Calpha_k%7D%29%7D%20%5Cprod_%7Bk%3D1%7D%5EK%20u_k%5E%7B%5Calpha_k%20-%201%7D%5C%5C%0A%26%20f_%7BDM%7D%28y%3B%20%5Calpha%29%20%3D%20%5C%0A%5Cint_%7Bu%20%5Cin%20%5CPhi%5E%7BK-1%7D%7D%20f_M%28y%3B%20u%29%20f_D%28u%3B%20%5Calpha%29%20%5C%5C%0A%26%20%3D%20%5Cfrac%20%7B%5CGamma%20%28%7B%5Csum_%7Bk%3D1%7D%5EK%20y_k%20%2B%201%7D%29%20%5CGamma%20%28%7B%5Csum_%7Bk%3D1%7D%5EK%20%5Calpha_k%7D%29%7D%20%0A%7B%5CGamma%20%28%7B%5Csum_%7Bk%3D1%7D%5EK%20y_k%7D%20%2B%20%7B%5Csum_%7Bk%3D1%7D%5EK%20%5Calpha_k%7D%29%7D%20%0A%5Cprod_%7Bk%3D1%7D%5EK%20%5Cfrac%20%7B%5CGamma%20%28y_k%20%2B%20%5Calpha_k%29%7D%20%0A%7B%5CGamma%20%28%7By_k%20%2B%201%7D%29%20%5CGamma%20%28%7B%5Calpha_k%7D%29%7D%0A%5Cend%7Bsplit%7D%20%0A "
\begin{split}
& a = (a_1, ..., a_K)^T,\ a_i > 0\\
& f_D (u;a) = \frac {\Gamma ({\sum_{k=1}^K \alpha_k})} 
{\prod_{k=1}^K \Gamma ({\sum_{k=1}^K \alpha_k})} \prod_{k=1}^K u_k^{\alpha_k - 1}\\
& f_{DM}(y; \alpha) = \
\int_{u \in \Phi^{K-1}} f_M(y; u) f_D(u; \alpha) \\
& = \frac {\Gamma ({\sum_{k=1}^K y_k + 1}) \Gamma ({\sum_{k=1}^K \alpha_k})} 
{\Gamma ({\sum_{k=1}^K y_k} + {\sum_{k=1}^K \alpha_k})} 
\prod_{k=1}^K \frac {\Gamma (y_k + \alpha_k)} 
{\Gamma ({y_k + 1}) \Gamma ({\alpha_k})}
\end{split} 
")

#### Limitation

-   all componnets must share a common variance parameter

-   components are mutually independent, up to the constraint that must
    sum up to 1

-   distribution fails to take into account the special and inherent
    property of microbiome count data (evolutionary relationships in the
    phylogenetic tree)

#### Generalization of Dirichlet Multinomial distribution

-   the relationships among the components of the count vector can be
    represented as a tree, node-by-node.

    -   each component has a independent variance

    -   components are correlated at subtree levels

-   a regression model with the effects of covariates

-   a regularized methods for selecting covarites (nutrients) that are
    associated with the count responses (OTUs)

#### Logistic Normal Multinomial distribution

-   Billheimer with Aitchison’s logistic normal distribution instead of
    Dirichlet

    -   covaraites to the count vector (Billheimer et al. 2001)
    -   link dietary nutrients with bacteria counts (Xia et al. 2013)

cannot exploit the tree structure information, logistic normal
multinomial does not have a closed-form expression.

### Dirichlet-Tree Multinomial Distributions

total number of counts, determined by the sequence depth, as an
ancillary statistics

-   analyssi conditioning on this number

A Tree
![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")
representing the hierarchical structure over the count responses

-   we will incorperate the tree structural information into the
    modeling

-   ![\\mathcal L](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%20L "\mathcal L")
    the set of leaf nodes

-   ![\\mathcal V](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%20V "\mathcal V")
    the set of internal nodes of
    ![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")

-   for each
    ![v \\in \\mathcal V](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v%20%5Cin%20%5Cmathcal%20V "v \in \mathcal V"),
    ![\\mathcal C_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%20C_v "\mathcal C_v")
    the child nodes of
    ![v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v "v")

-   for each
    ![v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v "v")
    and
    ![c \\in \\mathcal C_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%5Cin%20%5Cmathcal%20C_v "c \in \mathcal C_v"),
    ![\\delta\_{vc} (l) = 1(v\\rightarrow c, l \\in \\mathcal L)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_%7Bvc%7D%20%28l%29%20%3D%201%28v%5Crightarrow%20c%2C%20l%20%5Cin%20%5Cmathcal%20L%29 "\delta_{vc} (l) = 1(v\rightarrow c, l \in \mathcal L)")

-   ![\\mathcal L = \\{ 1, ..., K\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%20L%20%3D%20%5C%7B%201%2C%20...%2C%20K%5C%7D "\mathcal L = \{ 1, ..., K\}")

-   ![y\_{vc} = \\sum\_{l\\in\\mathcal L} \\delta\_{vc} (l) \\ y_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bvc%7D%20%3D%20%5Csum_%7Bl%5Cin%5Cmathcal%20L%7D%20%5Cdelta_%7Bvc%7D%20%28l%29%20%5C%20y_l "y_{vc} = \sum_{l\in\mathcal L} \delta_{vc} (l) \ y_l")

-   ![p\_{vc} = \\sum\_{l\\in \\mathcal L} \\delta\_{vc} (l) \\ p_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_%7Bvc%7D%20%3D%20%5Csum_%7Bl%5Cin%20%5Cmathcal%20L%7D%20%5Cdelta_%7Bvc%7D%20%28l%29%20%5C%20p_l "p_{vc} = \sum_{l\in \mathcal L} \delta_{vc} (l) \ p_l")

-   ![c \\in \\mathcal C_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%5Cin%20%5Cmathcal%20C_v "c \in \mathcal C_v")
    is the index for the subtree counts and probabilities

## 

-   Let
    ![v_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_0 "v_0")
    to be the root node of
    ![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T").

-   each leaf node is
    ![(v_0 \\equiv v_0^l )\\rightarrow v_1^l \\rightarrow ... \\rightarrow v\_{D_l -1}^l \\rightarrow (v\_{D_l}^l \\equiv l)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28v_0%20%5Cequiv%20v_0%5El%20%29%5Crightarrow%20v_1%5El%20%5Crightarrow%20...%20%5Crightarrow%20v_%7BD_l%20-1%7D%5El%20%5Crightarrow%20%28v_%7BD_l%7D%5El%20%5Cequiv%20l%29 "(v_0 \equiv v_0^l )\rightarrow v_1^l \rightarrow ... \rightarrow v_{D_l -1}^l \rightarrow (v_{D_l}^l \equiv l)")

    -   the path from
        ![v_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_0 "v_0")
        to
        ![l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;l "l"),
        where
        ![v_d^l \\in \\mathcal V, \\forall d = 0, ..., D\_{l-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_d%5El%20%5Cin%20%5Cmathcal%20V%2C%20%5Cforall%20d%20%3D%200%2C%20...%2C%20D_%7Bl-1%7D "v_d^l \in \mathcal V, \forall d = 0, ..., D_{l-1}")
    -   ![D_l \\geq 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D_l%20%5Cgeq%201 "D_l \geq 1")
        is the inner branches in this path
    -   ![b\_{vc} = \\frac {p\_{vc}} {\\sum\_{S \\in \\mathcal C_v} p\_{vS}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_%7Bvc%7D%20%3D%20%5Cfrac%20%7Bp_%7Bvc%7D%7D%20%7B%5Csum_%7BS%20%5Cin%20%5Cmathcal%20C_v%7D%20p_%7BvS%7D%7D "b_{vc} = \frac {p_{vc}} {\sum_{S \in \mathcal C_v} p_{vS}}"),
        with the
        ![\\sum\_{c \\in \\mathcal C_v} b\_{vc} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20b_%7Bvc%7D%20%3D%201 "\sum_{c \in \mathcal C_v} b_{vc} = 1")
        that each branch as a probability (**a conditional probability
        given we know where the branch started**).

![
p\_{l} = b\_{v_0 v_1^l} \\times b\_{v_1^l v_2^l} \\times ... \\times 
b\_{v\_{D_l - 1}^l} = \\prod\_{v\\in \\mathcal V} \\prod\_{c \\in \\mathcal C} b\_{vc}^{\\delta^{vc}(l)}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ap_%7Bl%7D%20%3D%20b_%7Bv_0%20v_1%5El%7D%20%5Ctimes%20b_%7Bv_1%5El%20v_2%5El%7D%20%5Ctimes%20...%20%5Ctimes%20%0Ab_%7Bv_%7BD_l%20-%201%7D%5El%7D%20%3D%20%5Cprod_%7Bv%5Cin%20%5Cmathcal%20V%7D%20%5Cprod_%7Bc%20%5Cin%20%5Cmathcal%20C%7D%20b_%7Bvc%7D%5E%7B%5Cdelta%5E%7Bvc%7D%28l%29%7D%0A "
p_{l} = b_{v_0 v_1^l} \times b_{v_1^l v_2^l} \times ... \times 
b_{v_{D_l - 1}^l} = \prod_{v\in \mathcal V} \prod_{c \in \mathcal C} b_{vc}^{\delta^{vc}(l)}
")

-   for
    ![p_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_l "p_l")
    is the product of branch probabilities from root
    ![v_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_0 "v_0")
    to final leaf
    ![l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;l "l")

**this is for the leaf level**

![
f_M(y; p) = f_M(y; b_v, v\\in \\mathcal V) = \\prod\_{v \\in \\mathcal V}
\\frac {\\Gamma (\\sum\_{c \\in \\mathcal C_v} y\_{vc} + 1)} 
{\\prod\_{c \\in \\mathcal C_v} \\Gamma ({y\_{vc} + 1})} \\prod\_{c \\in \\mathcal C_v} b\_{vc}^{y_vc}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Af_M%28y%3B%20p%29%20%3D%20f_M%28y%3B%20b_v%2C%20v%5Cin%20%5Cmathcal%20V%29%20%3D%20%5Cprod_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%0A%5Cfrac%20%7B%5CGamma%20%28%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20y_%7Bvc%7D%20%2B%201%29%7D%20%0A%7B%5Cprod_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20%5CGamma%20%28%7By_%7Bvc%7D%20%2B%201%7D%29%7D%20%5Cprod_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20b_%7Bvc%7D%5E%7By_vc%7D%0A "
f_M(y; p) = f_M(y; b_v, v\in \mathcal V) = \prod_{v \in \mathcal V}
\frac {\Gamma (\sum_{c \in \mathcal C_v} y_{vc} + 1)} 
{\prod_{c \in \mathcal C_v} \Gamma ({y_{vc} + 1})} \prod_{c \in \mathcal C_v} b_{vc}^{y_vc}
")

![b_v = \\{b\_{vc}, c\\in \\mathcal C_v\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_v%20%3D%20%5C%7Bb_%7Bvc%7D%2C%20c%5Cin%20%5Cmathcal%20C_v%5C%7D "b_v = \{b_{vc}, c\in \mathcal C_v\}")

**this is for the internal nodes**
![b_v, v \\in \\mathcal V](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_v%2C%20v%20%5Cin%20%5Cmathcal%20V "b_v, v \in \mathcal V")
the joint density of
![(b_v, v \\in \\mathcal V)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28b_v%2C%20v%20%5Cin%20%5Cmathcal%20V%29 "(b_v, v \in \mathcal V)")
has the form of:

![
f\_{DT} (u_v; \\alpha_v) = \\prod\_{v\\in \\mathcal V} f_D(u_v; \\alpha_v) = 
\\prod\_{v \\in \\mathcal V}
\\frac {\\Gamma (\\sum\_{c \\in \\mathcal C_v} \\alpha\_{vc})} 
{\\prod\_{c \\in \\mathcal C_v} \\Gamma ({\\alpha\_{vc}})} \\prod\_{c \\in \\mathcal C_v} u\_{vc}^{\\alpha_vc - 1}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Af_%7BDT%7D%20%28u_v%3B%20%5Calpha_v%29%20%3D%20%5Cprod_%7Bv%5Cin%20%5Cmathcal%20V%7D%20f_D%28u_v%3B%20%5Calpha_v%29%20%3D%20%0A%5Cprod_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%0A%5Cfrac%20%7B%5CGamma%20%28%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20%5Calpha_%7Bvc%7D%29%7D%20%0A%7B%5Cprod_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20%5CGamma%20%28%7B%5Calpha_%7Bvc%7D%7D%29%7D%20%5Cprod_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20u_%7Bvc%7D%5E%7B%5Calpha_vc%20-%201%7D%0A "
f_{DT} (u_v; \alpha_v) = \prod_{v\in \mathcal V} f_D(u_v; \alpha_v) = 
\prod_{v \in \mathcal V}
\frac {\Gamma (\sum_{c \in \mathcal C_v} \alpha_{vc})} 
{\prod_{c \in \mathcal C_v} \Gamma ({\alpha_{vc}})} \prod_{c \in \mathcal C_v} u_{vc}^{\alpha_vc - 1}
")

-   ![u_v = (u\_{vc}, c\\in \\mathcal C_v) \\in \\Phi^{K_v -1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_v%20%3D%20%28u_%7Bvc%7D%2C%20c%5Cin%20%5Cmathcal%20C_v%29%20%5Cin%20%5CPhi%5E%7BK_v%20-1%7D "u_v = (u_{vc}, c\in \mathcal C_v) \in \Phi^{K_v -1}")

-   ![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
    is the number of children of
    ![v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v "v")

-   ![\\alpha_v = (\\alpha_vc \> 0, c \\in \\mathcal C_v)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_v%20%3D%20%28%5Calpha_vc%20%3E%200%2C%20c%20%5Cin%20%5Cmathcal%20C_v%29 "\alpha_v = (\alpha_vc > 0, c \in \mathcal C_v)")

![
\\begin{split}
f\_{DTM} (y; \\alpha_v, v\\in\\mathcal V) 
& = \\prod\_{v\\in\\mathcal V} \\int\_{u_v \\in \\Phi^{K\_{v}-1}} f_M(y; u_v) f\_{D}(u_v; \\alpha_v) du_v \\\\
& = \\prod\_{v\\in\\mathcal V} \\frac {\\Gamma ({\\sum\_{c\\in \\mathcal C_v} y\_{vc} + 1}) \\Gamma ({\\sum\_{c\\in \\mathcal C_v} \\alpha\_{vc}})} 
{\\Gamma ({\\sum\_{c \\in \\mathcal C_v} y\_{vc} } + {\\sum\_{c\\in \\mathcal C_v}\\alpha\_{vc}})} 
\\prod\_{c\\in \\mathcal C_v} \\frac {\\Gamma (y\_{vc} + \\alpha\_{vc})} 
{\\Gamma ({y\_{vc} + 1}) \\Gamma ({\\alpha\_{vc}})}
\\end{split}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bsplit%7D%0Af_%7BDTM%7D%20%28y%3B%20%5Calpha_v%2C%20v%5Cin%5Cmathcal%20V%29%20%0A%26%20%3D%20%5Cprod_%7Bv%5Cin%5Cmathcal%20V%7D%20%5Cint_%7Bu_v%20%5Cin%20%5CPhi%5E%7BK_%7Bv%7D-1%7D%7D%20f_M%28y%3B%20u_v%29%20f_%7BD%7D%28u_v%3B%20%5Calpha_v%29%20du_v%20%5C%5C%0A%26%20%3D%20%5Cprod_%7Bv%5Cin%5Cmathcal%20V%7D%20%5Cfrac%20%7B%5CGamma%20%28%7B%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20y_%7Bvc%7D%20%2B%201%7D%29%20%5CGamma%20%28%7B%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Calpha_%7Bvc%7D%7D%29%7D%20%0A%7B%5CGamma%20%28%7B%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20y_%7Bvc%7D%20%7D%20%2B%20%7B%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%5Calpha_%7Bvc%7D%7D%29%7D%20%0A%5Cprod_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Cfrac%20%7B%5CGamma%20%28y_%7Bvc%7D%20%2B%20%5Calpha_%7Bvc%7D%29%7D%20%0A%7B%5CGamma%20%28%7By_%7Bvc%7D%20%2B%201%7D%29%20%5CGamma%20%28%7B%5Calpha_%7Bvc%7D%7D%29%7D%0A%5Cend%7Bsplit%7D%0A "
\begin{split}
f_{DTM} (y; \alpha_v, v\in\mathcal V) 
& = \prod_{v\in\mathcal V} \int_{u_v \in \Phi^{K_{v}-1}} f_M(y; u_v) f_{D}(u_v; \alpha_v) du_v \\
& = \prod_{v\in\mathcal V} \frac {\Gamma ({\sum_{c\in \mathcal C_v} y_{vc} + 1}) \Gamma ({\sum_{c\in \mathcal C_v} \alpha_{vc}})} 
{\Gamma ({\sum_{c \in \mathcal C_v} y_{vc} } + {\sum_{c\in \mathcal C_v}\alpha_{vc}})} 
\prod_{c\in \mathcal C_v} \frac {\Gamma (y_{vc} + \alpha_{vc})} 
{\Gamma ({y_{vc} + 1}) \Gamma ({\alpha_{vc}})}
\end{split}
")

that each component in the product
![\\prod\_{v \\in \\mathcal V}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cprod_%7Bv%20%5Cin%20%5Cmathcal%20V%7D "\prod_{v \in \mathcal V}")
corresponds to a interior node in the tree, which is a Dirichlet
multinomial distribution based on the accumulated counts along branches
of given node.

![
E\[p_l\] = \\prod\_{v\\in\\mathcal V} \\prod\_{c\\in\\mathcal C_v}
\\bigg(E\[b\_{vc}\] \\bigg)^{\\delta\_{vc}(l)} = 
\\prod\_{v\\in\\mathcal V} 
\\prod\_{c\\in\\mathcal C_v}\\bigg( \\frac {\\alpha\_{vc}} {\\sum\_{c\\in \\mathcal C_v} \\alpha\_{vc}} \\bigg)^{\\delta\_{vc}(l)}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AE%5Bp_l%5D%20%3D%20%5Cprod_%7Bv%5Cin%5Cmathcal%20V%7D%20%5Cprod_%7Bc%5Cin%5Cmathcal%20C_v%7D%0A%5Cbigg%28E%5Bb_%7Bvc%7D%5D%20%5Cbigg%29%5E%7B%5Cdelta_%7Bvc%7D%28l%29%7D%20%3D%20%0A%5Cprod_%7Bv%5Cin%5Cmathcal%20V%7D%20%0A%5Cprod_%7Bc%5Cin%5Cmathcal%20C_v%7D%5Cbigg%28%20%5Cfrac%20%7B%5Calpha_%7Bvc%7D%7D%20%7B%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Calpha_%7Bvc%7D%7D%20%5Cbigg%29%5E%7B%5Cdelta_%7Bvc%7D%28l%29%7D%0A "
E[p_l] = \prod_{v\in\mathcal V} \prod_{c\in\mathcal C_v}
\bigg(E[b_{vc}] \bigg)^{\delta_{vc}(l)} = 
\prod_{v\in\mathcal V} 
\prod_{c\in\mathcal C_v}\bigg( \frac {\alpha_{vc}} {\sum_{c\in \mathcal C_v} \alpha_{vc}} \bigg)^{\delta_{vc}(l)}
")

given
![\\sum\_{k=1}^K Y_k = \\sum\_{k=1}^K y_k; \\ E\[Y_l\] = E\[p_l\] \\sum\_{k=1} ^Ky_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bk%3D1%7D%5EK%20Y_k%20%3D%20%5Csum_%7Bk%3D1%7D%5EK%20y_k%3B%20%5C%20E%5BY_l%5D%20%3D%20E%5Bp_l%5D%20%5Csum_%7Bk%3D1%7D%20%5EKy_k "\sum_{k=1}^K Y_k = \sum_{k=1}^K y_k; \ E[Y_l] = E[p_l] \sum_{k=1} ^Ky_k")

Tips: 1. when
![\\mathcal V \\neq \\{v_0\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%20V%20%5Cneq%20%5C%7Bv_0%5C%7D "\mathcal V \neq \{v_0\}"),
each
![p_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_l "p_l")
has an independet variance 2. two components $p\_{v_i}, p\_{v_j} $ with
in the same subtree indexed by
![v \\in \\mathcal V](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v%20%5Cin%20%5Cmathcal%20V "v \in \mathcal V")
are correlated, due to the same ancestor of
![v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v "v").

### Dirichlet Tree Multinomial Regression Model

a regression of
![Y = (y_1, ..., y_K)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%20%3D%20%28y_1%2C%20...%2C%20y_K%29%5ET "Y = (y_1, ..., y_K)^T")
on
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
covariates
![X_1, ... X_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_1%2C%20...%20X_p "X_1, ... X_p")

for each
![v \\in \\mathcal V](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v%20%5Cin%20%5Cmathcal%20V "v \in \mathcal V")
and
![c \\in \\mathcal C_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%5Cin%20%5Cmathcal%20C_v "c \in \mathcal C_v")
we express
![\\log(\\alpha\_{vc})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clog%28%5Calpha_%7Bvc%7D%29 "\log(\alpha_{vc})")
as a linear combination of the covariates.

![
\\begin{split}
& \\log(\\alpha\_{vc}) = X^T \\beta\_{vc} \\ \\ \\ \\ (2.a)\\\\
& \\log(\\alpha\_{(i)vc}) = X_i^T \\beta\_{(i)vc} + Z_i^T b\_{(i){vc}}\\ \\ \\ \\ (2.b)\\\\
& X = (X_1, ..., X_p)^T\\\\
& \\beta\_{vc} = (\\beta\_{vc1}, ..., \\beta\_{vcp})^T\\\\
& y_i = (y\_{i1}, ..., y\_{iK})^T\\\\
& x_i = (x\_{i1}, ..., x\_{ip})^T\\\\
& b\_{vc} = (b\_{vc1}, ..., b\_{vcq})^T\\\\
& z_i = (z\_{i1}, ..., z\_{iq})^T\\\\
\\\\
& \\beta_v = (\\beta\_{vc}, c \\in \\mathcal C_v)\\\\
& \\beta = (\\beta_v, v \\in \\mathcal V) \\\\
\\\\
& \\beta\_{vi} \| b\_{vi} = (\\beta\_{(i)vc} \| b\_{(i)vc}, c \\in \\mathcal C_v)\\\\
& \\beta_i \| b\_{i} = (\\beta\_{(i)v}\| b\_{(i)v}, v \\in \\mathcal V)
\\end{split}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bsplit%7D%0A%26%20%5Clog%28%5Calpha_%7Bvc%7D%29%20%3D%20X%5ET%20%5Cbeta_%7Bvc%7D%20%5C%20%5C%20%5C%20%5C%20%282.a%29%5C%5C%0A%26%20%5Clog%28%5Calpha_%7B%28i%29vc%7D%29%20%3D%20X_i%5ET%20%5Cbeta_%7B%28i%29vc%7D%20%2B%20Z_i%5ET%20b_%7B%28i%29%7Bvc%7D%7D%5C%20%5C%20%5C%20%5C%20%282.b%29%5C%5C%0A%26%20X%20%3D%20%28X_1%2C%20...%2C%20X_p%29%5ET%5C%5C%0A%26%20%5Cbeta_%7Bvc%7D%20%3D%20%28%5Cbeta_%7Bvc1%7D%2C%20...%2C%20%5Cbeta_%7Bvcp%7D%29%5ET%5C%5C%0A%26%20y_i%20%3D%20%28y_%7Bi1%7D%2C%20...%2C%20y_%7BiK%7D%29%5ET%5C%5C%0A%26%20x_i%20%3D%20%28x_%7Bi1%7D%2C%20...%2C%20x_%7Bip%7D%29%5ET%5C%5C%0A%26%20b_%7Bvc%7D%20%3D%20%28b_%7Bvc1%7D%2C%20...%2C%20b_%7Bvcq%7D%29%5ET%5C%5C%0A%26%20z_i%20%3D%20%28z_%7Bi1%7D%2C%20...%2C%20z_%7Biq%7D%29%5ET%5C%5C%0A%5C%5C%0A%26%20%5Cbeta_v%20%3D%20%28%5Cbeta_%7Bvc%7D%2C%20c%20%5Cin%20%5Cmathcal%20C_v%29%5C%5C%0A%26%20%5Cbeta%20%3D%20%28%5Cbeta_v%2C%20v%20%5Cin%20%5Cmathcal%20V%29%20%5C%5C%0A%5C%5C%0A%26%20%5Cbeta_%7Bvi%7D%20%7C%20b_%7Bvi%7D%20%3D%20%28%5Cbeta_%7B%28i%29vc%7D%20%7C%20b_%7B%28i%29vc%7D%2C%20c%20%5Cin%20%5Cmathcal%20C_v%29%5C%5C%0A%26%20%5Cbeta_i%20%7C%20b_%7Bi%7D%20%3D%20%28%5Cbeta_%7B%28i%29v%7D%7C%20b_%7B%28i%29v%7D%2C%20v%20%5Cin%20%5Cmathcal%20V%29%0A%5Cend%7Bsplit%7D%0A "
\begin{split}
& \log(\alpha_{vc}) = X^T \beta_{vc} \ \ \ \ (2.a)\\
& \log(\alpha_{(i)vc}) = X_i^T \beta_{(i)vc} + Z_i^T b_{(i){vc}}\ \ \ \ (2.b)\\
& X = (X_1, ..., X_p)^T\\
& \beta_{vc} = (\beta_{vc1}, ..., \beta_{vcp})^T\\
& y_i = (y_{i1}, ..., y_{iK})^T\\
& x_i = (x_{i1}, ..., x_{ip})^T\\
& b_{vc} = (b_{vc1}, ..., b_{vcq})^T\\
& z_i = (z_{i1}, ..., z_{iq})^T\\
\\
& \beta_v = (\beta_{vc}, c \in \mathcal C_v)\\
& \beta = (\beta_v, v \in \mathcal V) \\
\\
& \beta_{vi} | b_{vi} = (\beta_{(i)vc} | b_{(i)vc}, c \in \mathcal C_v)\\
& \beta_i | b_{i} = (\beta_{(i)v}| b_{(i)v}, v \in \mathcal V)
\end{split}
")

![
l\_{DTM}(\\beta) = \\log L\_{DTM}(\\beta) = \\sum\_{v \\in \\mathcal V} \\log L_v(\\beta_v) 
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Al_%7BDTM%7D%28%5Cbeta%29%20%3D%20%5Clog%20L_%7BDTM%7D%28%5Cbeta%29%20%3D%20%5Csum_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%20%5Clog%20L_v%28%5Cbeta_v%29%20%0A "
l_{DTM}(\beta) = \log L_{DTM}(\beta) = \sum_{v \in \mathcal V} \log L_v(\beta_v) 
")

![
l_v(\\beta_v) = log L_v(\\beta_v) = 
\\sum\_{i=1}^n \\Bigg\[
  \\tilde \\Gamma \\Big(\\sum\_{c\\in \\mathcal C_v} \\alpha\_{ivc}\\Big)
  - \\tilde \\Gamma \\Big(\\sum\_{c\\in \\mathcal C_v} y\_{ivc} + 
    \\sum\_{c\\in \\mathcal C_v} \\alpha\_{ivc}\\Big) +
\\sum\_{c \\in \\mathcal C_v} 
  \\bigg\\{ \\tilde \\Gamma \\Big(y\_{ivc} +
    \\alpha\_{ivc}\\Big) -
\\tilde \\Gamma \\Big(\\alpha\_{ivc}\\Big) \\bigg\\}
\\Bigg\]
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Al_v%28%5Cbeta_v%29%20%3D%20log%20L_v%28%5Cbeta_v%29%20%3D%20%0A%5Csum_%7Bi%3D1%7D%5En%20%5CBigg%5B%0A%20%20%5Ctilde%20%5CGamma%20%5CBig%28%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Calpha_%7Bivc%7D%5CBig%29%0A%20%20-%20%5Ctilde%20%5CGamma%20%5CBig%28%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20y_%7Bivc%7D%20%2B%20%0A%20%20%20%20%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Calpha_%7Bivc%7D%5CBig%29%20%2B%0A%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20%0A%20%20%5Cbigg%5C%7B%20%5Ctilde%20%5CGamma%20%5CBig%28y_%7Bivc%7D%20%2B%0A%20%20%20%20%5Calpha_%7Bivc%7D%5CBig%29%20-%0A%5Ctilde%20%5CGamma%20%5CBig%28%5Calpha_%7Bivc%7D%5CBig%29%20%5Cbigg%5C%7D%0A%5CBigg%5D%0A "
l_v(\beta_v) = log L_v(\beta_v) = 
\sum_{i=1}^n \Bigg[
  \tilde \Gamma \Big(\sum_{c\in \mathcal C_v} \alpha_{ivc}\Big)
  - \tilde \Gamma \Big(\sum_{c\in \mathcal C_v} y_{ivc} + 
    \sum_{c\in \mathcal C_v} \alpha_{ivc}\Big) +
\sum_{c \in \mathcal C_v} 
  \bigg\{ \tilde \Gamma \Big(y_{ivc} +
    \alpha_{ivc}\Big) -
\tilde \Gamma \Big(\alpha_{ivc}\Big) \bigg\}
\Bigg]
")

![\\tilde \\Gamma (.) = \\log \\big(\\Gamma (.) \\big)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctilde%20%5CGamma%20%28.%29%20%3D%20%5Clog%20%5Cbig%28%5CGamma%20%28.%29%20%5Cbig%29 "\tilde \Gamma (.) = \log \big(\Gamma (.) \big)")

![\\alpha\_{ivc} = \\exp(x_i^T\\beta\_{vc})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%7Bivc%7D%20%3D%20%5Cexp%28x_i%5ET%5Cbeta_%7Bvc%7D%29 "\alpha_{ivc} = \exp(x_i^T\beta_{vc})")

![y\_{ivc} = \\sum\_{l \\in \\mathcal L} \\delta\_{vc} (l) y\_{il}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bivc%7D%20%3D%20%5Csum_%7Bl%20%5Cin%20%5Cmathcal%20L%7D%20%5Cdelta_%7Bvc%7D%20%28l%29%20y_%7Bil%7D "y_{ivc} = \sum_{l \in \mathcal L} \delta_{vc} (l) y_{il}")

![
l_v(\\beta_v) = \\log L_v(\\beta_v) = 
\\sum\_{i=1}^n \\Bigg\[
  \\log \\bigg(\\Gamma \\Big(\\sum\_{c\\in \\mathcal C_v} \\exp(x_i^T\\beta\_{vc})\\Big) \\bigg) 
  - \\log \\bigg(\\Gamma \\Big(\\sum\_{c\\in \\mathcal C_v} y\_{ivc} + 
    \\sum\_{c\\in \\mathcal C_v} \\exp(x_i^T\\beta\_{vc})\\Big) \\bigg) +
\\sum\_{c \\in \\mathcal C_v} 
  \\bigg\\{ \\log \\bigg(\\Gamma \\Big(y\_{ivc} +
    \\exp(x_i^T\\beta\_{vc})\\Big)\\bigg) -
\\log \\bigg(\\Gamma \\Big(\\exp(x_i^T\\beta\_{vc})\\Big) \\bigg)\\bigg\\}
\\Bigg\]
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Al_v%28%5Cbeta_v%29%20%3D%20%5Clog%20L_v%28%5Cbeta_v%29%20%3D%20%0A%5Csum_%7Bi%3D1%7D%5En%20%5CBigg%5B%0A%20%20%5Clog%20%5Cbigg%28%5CGamma%20%5CBig%28%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Cexp%28x_i%5ET%5Cbeta_%7Bvc%7D%29%5CBig%29%20%5Cbigg%29%20%0A%20%20-%20%5Clog%20%5Cbigg%28%5CGamma%20%5CBig%28%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20y_%7Bivc%7D%20%2B%20%0A%20%20%20%20%5Csum_%7Bc%5Cin%20%5Cmathcal%20C_v%7D%20%5Cexp%28x_i%5ET%5Cbeta_%7Bvc%7D%29%5CBig%29%20%5Cbigg%29%20%2B%0A%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%20%0A%20%20%5Cbigg%5C%7B%20%5Clog%20%5Cbigg%28%5CGamma%20%5CBig%28y_%7Bivc%7D%20%2B%0A%20%20%20%20%5Cexp%28x_i%5ET%5Cbeta_%7Bvc%7D%29%5CBig%29%5Cbigg%29%20-%0A%5Clog%20%5Cbigg%28%5CGamma%20%5CBig%28%5Cexp%28x_i%5ET%5Cbeta_%7Bvc%7D%29%5CBig%29%20%5Cbigg%29%5Cbigg%5C%7D%0A%5CBigg%5D%0A "
l_v(\beta_v) = \log L_v(\beta_v) = 
\sum_{i=1}^n \Bigg[
  \log \bigg(\Gamma \Big(\sum_{c\in \mathcal C_v} \exp(x_i^T\beta_{vc})\Big) \bigg) 
  - \log \bigg(\Gamma \Big(\sum_{c\in \mathcal C_v} y_{ivc} + 
    \sum_{c\in \mathcal C_v} \exp(x_i^T\beta_{vc})\Big) \bigg) +
\sum_{c \in \mathcal C_v} 
  \bigg\{ \log \bigg(\Gamma \Big(y_{ivc} +
    \exp(x_i^T\beta_{vc})\Big)\bigg) -
\log \bigg(\Gamma \Big(\exp(x_i^T\beta_{vc})\Big) \bigg)\bigg\}
\Bigg]
")

### Regularized Likelihood Estimation

parametrization for Dirichlet Tree Multinomial Regression Model

number of
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
time the number of tree branches
![\\sum\_{v\\in \\mathcal V} K_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bv%5Cin%20%5Cmathcal%20V%7D%20K_v "\sum_{v\in \mathcal V} K_v")

![
pnl\_{DTM}(\\beta; \\lambda, \\gamma) = -l\_{DTM}(\\beta) 
 + \\lambda \\bigg\\{ (1- \\gamma)\\sum\_{v \\in \\mathcal V} \\sum\_{c \\in \\mathcal C_v}\\\|\\beta\_{cv}\\\|\_{L1}
 + \\gamma \\sum\_{v \\in \\mathcal V} \\sum\_{c \\in \\mathcal C_v}\\\|\\beta\_{cv}\\\|\_{L2} \\bigg\\}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Apnl_%7BDTM%7D%28%5Cbeta%3B%20%5Clambda%2C%20%5Cgamma%29%20%3D%20-l_%7BDTM%7D%28%5Cbeta%29%20%0A%20%2B%20%5Clambda%20%5Cbigg%5C%7B%20%281-%20%5Cgamma%29%5Csum_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%20%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%5C%7C%5Cbeta_%7Bcv%7D%5C%7C_%7BL1%7D%0A%20%2B%20%5Cgamma%20%5Csum_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%20%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%5C%7C%5Cbeta_%7Bcv%7D%5C%7C_%7BL2%7D%20%5Cbigg%5C%7D%0A "
pnl_{DTM}(\beta; \lambda, \gamma) = -l_{DTM}(\beta) 
 + \lambda \bigg\{ (1- \gamma)\sum_{v \in \mathcal V} \sum_{c \in \mathcal C_v}\|\beta_{cv}\|_{L1}
 + \gamma \sum_{v \in \mathcal V} \sum_{c \in \mathcal C_v}\|\beta_{cv}\|_{L2} \bigg\}
")

### Algorithm of accelerated proximal gradient method

approximate
![l\_{DTM} (\\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;l_%7BDTM%7D%20%28%5Cbeta%29 "l_{DTM} (\beta)")
at
![\\beta^{(t)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%5E%7B%28t%29%7D "\beta^{(t)}")

![
\\begin{split}
& - \\tilde l\_{DTM}(\\beta; \\eta^{(t)}) = -l\_{DTM}(\\eta^{(t)}) - \\langle\\beta - \\eta^{(t)}, \\nabla l\_{DTM} (\\eta^{(t)})\\rangle + \\frac {C} {2} \\\|\\beta - \\eta^{(t)}\\\|^2_2 \\\\
& pn \\tilde{l}\_{DMT}(\\beta) = - \\tilde l\_{DTM} (\\beta; \\eta^{(t)}) + \\lambda_1
\\sum\_{v \\in \\mathcal V} \\sum\_{c \\in \\mathcal C_v}\\\|\\beta\_{cv}\\\|\_{1}
 + \\lambda_2\\sum\_{v \\in \\mathcal V} \\sum\_{c \\in \\mathcal C_v}\\\|\\beta\_{cv}\\\|\_{2} \\\\
& \\beta^{(t+1)} = \\underset {\\beta} {arg min}  \\big(pn \\tilde{l}\_{DTM} (\\beta)\\big) \\\\
& \\eta^{(t+1)} = \\beta^{(t+1)} + \\frac {1-\\alpha_t} {\\alpha_t} \\alpha\_{t+1} (\\beta^{(t+1)} - \\beta^{(t)}) \\\\
& \\alpha_t = 2/(t+2)
\\end{split}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bsplit%7D%0A%26%20-%20%5Ctilde%20l_%7BDTM%7D%28%5Cbeta%3B%20%5Ceta%5E%7B%28t%29%7D%29%20%3D%20-l_%7BDTM%7D%28%5Ceta%5E%7B%28t%29%7D%29%20-%20%5Clangle%5Cbeta%20-%20%5Ceta%5E%7B%28t%29%7D%2C%20%5Cnabla%20l_%7BDTM%7D%20%28%5Ceta%5E%7B%28t%29%7D%29%5Crangle%20%2B%20%5Cfrac%20%7BC%7D%20%7B2%7D%20%5C%7C%5Cbeta%20-%20%5Ceta%5E%7B%28t%29%7D%5C%7C%5E2_2%20%5C%5C%0A%26%20pn%20%5Ctilde%7Bl%7D_%7BDMT%7D%28%5Cbeta%29%20%3D%20-%20%5Ctilde%20l_%7BDTM%7D%20%28%5Cbeta%3B%20%5Ceta%5E%7B%28t%29%7D%29%20%2B%20%5Clambda_1%0A%5Csum_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%20%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%5C%7C%5Cbeta_%7Bcv%7D%5C%7C_%7B1%7D%0A%20%2B%20%5Clambda_2%5Csum_%7Bv%20%5Cin%20%5Cmathcal%20V%7D%20%5Csum_%7Bc%20%5Cin%20%5Cmathcal%20C_v%7D%5C%7C%5Cbeta_%7Bcv%7D%5C%7C_%7B2%7D%20%5C%5C%0A%26%20%5Cbeta%5E%7B%28t%2B1%29%7D%20%3D%20%5Cunderset%20%7B%5Cbeta%7D%20%7Barg%20min%7D%20%20%5Cbig%28pn%20%5Ctilde%7Bl%7D_%7BDTM%7D%20%28%5Cbeta%29%5Cbig%29%20%5C%5C%0A%26%20%5Ceta%5E%7B%28t%2B1%29%7D%20%3D%20%5Cbeta%5E%7B%28t%2B1%29%7D%20%2B%20%5Cfrac%20%7B1-%5Calpha_t%7D%20%7B%5Calpha_t%7D%20%5Calpha_%7Bt%2B1%7D%20%28%5Cbeta%5E%7B%28t%2B1%29%7D%20-%20%5Cbeta%5E%7B%28t%29%7D%29%20%5C%5C%0A%26%20%5Calpha_t%20%3D%202%2F%28t%2B2%29%0A%5Cend%7Bsplit%7D%0A "
\begin{split}
& - \tilde l_{DTM}(\beta; \eta^{(t)}) = -l_{DTM}(\eta^{(t)}) - \langle\beta - \eta^{(t)}, \nabla l_{DTM} (\eta^{(t)})\rangle + \frac {C} {2} \|\beta - \eta^{(t)}\|^2_2 \\
& pn \tilde{l}_{DMT}(\beta) = - \tilde l_{DTM} (\beta; \eta^{(t)}) + \lambda_1
\sum_{v \in \mathcal V} \sum_{c \in \mathcal C_v}\|\beta_{cv}\|_{1}
 + \lambda_2\sum_{v \in \mathcal V} \sum_{c \in \mathcal C_v}\|\beta_{cv}\|_{2} \\
& \beta^{(t+1)} = \underset {\beta} {arg min}  \big(pn \tilde{l}_{DTM} (\beta)\big) \\
& \eta^{(t+1)} = \beta^{(t+1)} + \frac {1-\alpha_t} {\alpha_t} \alpha_{t+1} (\beta^{(t+1)} - \beta^{(t)}) \\
& \alpha_t = 2/(t+2)
\end{split}
")

![
\\begin{split}
& \\beta\_{vc}^{(t+1)}  = \\underset {\\beta\_{vc}} {argmin} \\Bigg\[ \\frac 1 2
\\bigg \\\|
\\beta\_{vc} - \\Big\\{\\eta\_{vc}^{(t)} - \\frac 1 C \\nabla l\_{vc} (\\eta_v^{(t)}) \\Big\\}
\\bigg\\\|\_2^2 + \\frac {\\lambda_1} {C}\\\|\\beta\_{vc}\\\|\_1 + \\frac {\\lambda_2} {C} \\\|\\beta\_{vc}\\\|\_2
\\Bigg\]\\\\
 & \\text {s.t.} \\ \\ \\beta\_{vc}^{(t+1)}  = \\underset {\\beta\_{vc}} {argmin} \\Bigg\[ \\frac 1 2
\\bigg \\\|
\\beta\_{vc} +\\frac 1 C \\nabla l\_{vc} (0)
\\bigg\\\|\_2^2 + \\lambda \\bigg(\\frac {1-\\gamma} {C}\\\|\\beta\_{vc}\\\|\_1 + \\frac {\\gamma} {C} \\\|\\beta\_{vc}\\\|\_2\\bigg)
\\Bigg\]\\\\
& \\tilde \\beta\_{vc}^{(t+1)}  = sgn \\bigg\\{ \\eta\_{vc}^{(t)} - \\frac 1 C \\nabla l\_{vc} (\\eta\_{vc}^{(t)})\\bigg\\} \\ \\max \\bigg\\{0, \\Big\|\\eta\_{vc}^{(t)} - \\frac 1 C \\nabla l\_{vc}(\\eta_v^{(t)}) \\Big\| - \\frac {\\lambda_1} {C}\\bigg\\}\\\\
& \\beta\_{vc}^{(t+1)}  = \\frac {\\tilde \\beta\_{vc}^{(t+1)}}{\\\|\\tilde \\beta\_{vc}^{(t+1)}\\\|\_2}
\\max \\bigg\\{0, \\Big\\\|\\tilde \\beta\_{vc}^{(t+1)} \\Big\\\|\_2 - \\frac {\\lambda_1} {C}\\bigg\\} 
\\end{split}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bsplit%7D%0A%26%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%20%20%3D%20%5Cunderset%20%7B%5Cbeta_%7Bvc%7D%7D%20%7Bargmin%7D%20%5CBigg%5B%20%5Cfrac%201%202%0A%5Cbigg%20%5C%7C%0A%5Cbeta_%7Bvc%7D%20-%20%5CBig%5C%7B%5Ceta_%7Bvc%7D%5E%7B%28t%29%7D%20-%20%5Cfrac%201%20C%20%5Cnabla%20l_%7Bvc%7D%20%28%5Ceta_v%5E%7B%28t%29%7D%29%20%5CBig%5C%7D%0A%5Cbigg%5C%7C_2%5E2%20%2B%20%5Cfrac%20%7B%5Clambda_1%7D%20%7BC%7D%5C%7C%5Cbeta_%7Bvc%7D%5C%7C_1%20%2B%20%5Cfrac%20%7B%5Clambda_2%7D%20%7BC%7D%20%5C%7C%5Cbeta_%7Bvc%7D%5C%7C_2%0A%5CBigg%5D%5C%5C%0A%20%26%20%5Ctext%20%7Bs.t.%7D%20%5C%20%5C%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%20%20%3D%20%5Cunderset%20%7B%5Cbeta_%7Bvc%7D%7D%20%7Bargmin%7D%20%5CBigg%5B%20%5Cfrac%201%202%0A%5Cbigg%20%5C%7C%0A%5Cbeta_%7Bvc%7D%20%2B%5Cfrac%201%20C%20%5Cnabla%20l_%7Bvc%7D%20%280%29%0A%5Cbigg%5C%7C_2%5E2%20%2B%20%5Clambda%20%5Cbigg%28%5Cfrac%20%7B1-%5Cgamma%7D%20%7BC%7D%5C%7C%5Cbeta_%7Bvc%7D%5C%7C_1%20%2B%20%5Cfrac%20%7B%5Cgamma%7D%20%7BC%7D%20%5C%7C%5Cbeta_%7Bvc%7D%5C%7C_2%5Cbigg%29%0A%5CBigg%5D%5C%5C%0A%26%20%5Ctilde%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%20%20%3D%20sgn%20%5Cbigg%5C%7B%20%5Ceta_%7Bvc%7D%5E%7B%28t%29%7D%20-%20%5Cfrac%201%20C%20%5Cnabla%20l_%7Bvc%7D%20%28%5Ceta_%7Bvc%7D%5E%7B%28t%29%7D%29%5Cbigg%5C%7D%20%5C%20%5Cmax%20%5Cbigg%5C%7B0%2C%20%5CBig%7C%5Ceta_%7Bvc%7D%5E%7B%28t%29%7D%20-%20%5Cfrac%201%20C%20%5Cnabla%20l_%7Bvc%7D%28%5Ceta_v%5E%7B%28t%29%7D%29%20%5CBig%7C%20-%20%5Cfrac%20%7B%5Clambda_1%7D%20%7BC%7D%5Cbigg%5C%7D%5C%5C%0A%26%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%20%20%3D%20%5Cfrac%20%7B%5Ctilde%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%7D%7B%5C%7C%5Ctilde%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%5C%7C_2%7D%0A%5Cmax%20%5Cbigg%5C%7B0%2C%20%5CBig%5C%7C%5Ctilde%20%5Cbeta_%7Bvc%7D%5E%7B%28t%2B1%29%7D%20%5CBig%5C%7C_2%20-%20%5Cfrac%20%7B%5Clambda_1%7D%20%7BC%7D%5Cbigg%5C%7D%20%0A%5Cend%7Bsplit%7D%0A "
\begin{split}
& \beta_{vc}^{(t+1)}  = \underset {\beta_{vc}} {argmin} \Bigg[ \frac 1 2
\bigg \|
\beta_{vc} - \Big\{\eta_{vc}^{(t)} - \frac 1 C \nabla l_{vc} (\eta_v^{(t)}) \Big\}
\bigg\|_2^2 + \frac {\lambda_1} {C}\|\beta_{vc}\|_1 + \frac {\lambda_2} {C} \|\beta_{vc}\|_2
\Bigg]\\
 & \text {s.t.} \ \ \beta_{vc}^{(t+1)}  = \underset {\beta_{vc}} {argmin} \Bigg[ \frac 1 2
\bigg \|
\beta_{vc} +\frac 1 C \nabla l_{vc} (0)
\bigg\|_2^2 + \lambda \bigg(\frac {1-\gamma} {C}\|\beta_{vc}\|_1 + \frac {\gamma} {C} \|\beta_{vc}\|_2\bigg)
\Bigg]\\
& \tilde \beta_{vc}^{(t+1)}  = sgn \bigg\{ \eta_{vc}^{(t)} - \frac 1 C \nabla l_{vc} (\eta_{vc}^{(t)})\bigg\} \ \max \bigg\{0, \Big|\eta_{vc}^{(t)} - \frac 1 C \nabla l_{vc}(\eta_v^{(t)}) \Big| - \frac {\lambda_1} {C}\bigg\}\\
& \beta_{vc}^{(t+1)}  = \frac {\tilde \beta_{vc}^{(t+1)}}{\|\tilde \beta_{vc}^{(t+1)}\|_2}
\max \bigg\{0, \Big\|\tilde \beta_{vc}^{(t+1)} \Big\|_2 - \frac {\lambda_1} {C}\bigg\} 
\end{split}
")
