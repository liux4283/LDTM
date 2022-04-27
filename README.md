
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

$
$

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

$
\_k + 1}) ({\_k})} \\end{split} $

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

$ p\_{l} = b\_{v_0 v_1^l} b\_{v_1^l v_2^l} … b\_{v\_{D_l - 1}^l} = *{vV}
*{c C} b\_{vc}<sup>{</sup>{vc}(l)} $

-   for
    ![p_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_l "p_l")
    is the product of branch probabilities from root
    ![v_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_0 "v_0")
    to final leaf
    ![l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;l "l")

**this is for the leaf level**

$ f_M(y; p) = f_M(y; b_v, vV) = *{v V} {*{c C_v} ({y\_{vc} + 1})} *{c
C_v} b*{vc}^{y_vc} $

![b_v = \\{b\_{vc}, c\\in \\mathcal C_v\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_v%20%3D%20%5C%7Bb_%7Bvc%7D%2C%20c%5Cin%20%5Cmathcal%20C_v%5C%7D "b_v = \{b_{vc}, c\in \mathcal C_v\}")

**this is for the internal nodes**
![b_v, v \\in \\mathcal V](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_v%2C%20v%20%5Cin%20%5Cmathcal%20V "b_v, v \in \mathcal V")
the joint density of
![(b_v, v \\in \\mathcal V)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28b_v%2C%20v%20%5Cin%20%5Cmathcal%20V%29 "(b_v, v \in \mathcal V)")
has the form of:

$ f\_{DT} (u_v; *v) = *{vV} f_D(u_v; *v) = *{v V} {*{c C_v} ({*{vc}})}
*{c C_v} u*{vc}^{\_vc - 1} $

-   ![u_v = (u\_{vc}, c\\in \\mathcal C_v) \\in \\Phi^{K_v -1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_v%20%3D%20%28u_%7Bvc%7D%2C%20c%5Cin%20%5Cmathcal%20C_v%29%20%5Cin%20%5CPhi%5E%7BK_v%20-1%7D "u_v = (u_{vc}, c\in \mathcal C_v) \in \Phi^{K_v -1}")

-   ![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
    is the number of children of
    ![v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v "v")

-   ![\\alpha_v = (\\alpha_vc \> 0, c \\in \\mathcal C_v)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_v%20%3D%20%28%5Calpha_vc%20%3E%200%2C%20c%20%5Cin%20%5Cmathcal%20C_v%29 "\alpha_v = (\alpha_vc > 0, c \in \mathcal C_v)")

$
$

that each component in the product
![\\prod\_{v \\in \\mathcal V}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cprod_%7Bv%20%5Cin%20%5Cmathcal%20V%7D "\prod_{v \in \mathcal V}")
corresponds to a interior node in the tree, which is a Dirichlet
multinomial distribution based on the accumulated counts along branches
of given node.

$ E\[p_l\] = *{vV} *{cC_v} (E\[b\_{vc}\] )^{*{vc}(l)} = *{vV} *{cC_v}(
{*{cC_v} *{vc}} )^{*{vc}(l)} $

given
![\\sum\_{k=1}^K Y_k = \\sum\_{k=1}^K y_k; \\ E\[Y_l\] = E\[p_l\] \\sum\_{k=1} ^Ky_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bk%3D1%7D%5EK%20Y_k%20%3D%20%5Csum_%7Bk%3D1%7D%5EK%20y_k%3B%20%5C%20E%5BY_l%5D%20%3D%20E%5Bp_l%5D%20%5Csum_%7Bk%3D1%7D%20%5EKy_k "\sum_{k=1}^K Y_k = \sum_{k=1}^K y_k; \ E[Y_l] = E[p_l] \sum_{k=1} ^Ky_k")

Tips: 1. when
![\\mathcal V \\neq \\{v_0\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%20V%20%5Cneq%20%5C%7Bv_0%5C%7D "\mathcal V \neq \{v_0\}"),
each
![p_l](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_l "p_l")
has an independet variance

1.  two components $p\_{v_i}, p\_{v_j} $ with in the same subtree
    indexed by
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

$
$

$ l\_{DTM}() = L\_{DTM}() = \_{v V} L_v(\_v) $

$ l_v(\_v) = log L_v(*v) = *{i=1}^n $

![\\tilde \\Gamma (.) = \\log \\big(\\Gamma (.) \\big)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctilde%20%5CGamma%20%28.%29%20%3D%20%5Clog%20%5Cbig%28%5CGamma%20%28.%29%20%5Cbig%29 "\tilde \Gamma (.) = \log \big(\Gamma (.) \big)")

![\\alpha\_{ivc} = \\exp(x_i^T\\beta\_{vc})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%7Bivc%7D%20%3D%20%5Cexp%28x_i%5ET%5Cbeta_%7Bvc%7D%29 "\alpha_{ivc} = \exp(x_i^T\beta_{vc})")

![y\_{ivc} = \\sum\_{l \\in \\mathcal L} \\delta\_{vc} (l) y\_{il}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bivc%7D%20%3D%20%5Csum_%7Bl%20%5Cin%20%5Cmathcal%20L%7D%20%5Cdelta_%7Bvc%7D%20%28l%29%20y_%7Bil%7D "y_{ivc} = \sum_{l \in \mathcal L} \delta_{vc} (l) y_{il}")

$ l_v(\_v) = L_v(*v) = *{i=1}^n $

### Regularized Likelihood Estimation

parametrization for Dirichlet Tree Multinomial Regression Model

number of
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
time the number of tree branches
![\\sum\_{v\\in \\mathcal V} K_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bv%5Cin%20%5Cmathcal%20V%7D%20K_v "\sum_{v\in \mathcal V} K_v")

$ pnl\_{DTM}(; , ) = -l\_{DTM}() + { (1- )*{v V} *{c
C_v}\|*{cv}\|*{L1} + *{v V} *{c C_v}\|*{cv}\|*{L2} } $

### Algorithm of accelerated proximal gradient method

approximate
![l\_{DTM} (\\beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;l_%7BDTM%7D%20%28%5Cbeta%29 "l_{DTM} (\beta)")
at
![\\beta^{(t)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%5E%7B%28t%29%7D "\beta^{(t)}")

$
$

$
$
