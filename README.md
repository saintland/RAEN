README
================

# The Random Ensemble Variable Selection Method

The Proportional Subdistribution Hazard (PSH) model has been popular for
estimating the effects of the covariates on the cause of interest in
Competing Risks analysis. The fast accumulation of large scale datasets
has posed a challenge to classical statistical methods. Current
penalized variable selection methods show unsatisfactory performance in
ultra-high dimensional data. We propose a novel method, the Random
Approximate Elastic Net (RAEN), with a robust and generalized solution
to the variable selection problem for the PSH model. Our method shows
improved sensitivity for variable selection compared with current
methods.

# Installation

RAEN can be installed from R-CRAN

``` r
install.packages('RAEN')
```

Users can install the developmental version from Github.

``` r
library(devtools)
install_github('saintland/RAEN')
```

# Splitting correlated variables

The simulated data `toydata` contains 200 rows, time to event, censoring
status, and 1000 predictors, 60 of which are true predictors (\(X1-X20\)
and \(X40-X80\)). The variable correlation blocks are identified as the
following example.

``` r
require(RAEN,quietly = T)
```

    ## Loaded lars 1.2

``` r
data(toydata)
x=toydata[,-c(1:2)]
y=toydata[,1:2]
fgrp<-deCorr(x)
```

    ## 
    ## No: 1 cluster contains : 19 , remaining 981 
    ## 
    ## No: 2 cluster contains : 17 , remaining 964 
    ## 
    ## No: 3 cluster contains : 19 , remaining 945 
    ## 
    ## No: 4 cluster contains : 18 , remaining 927 
    ## 
    ## No: 5 cluster contains : 2 , remaining 925

The variable selection is executed via the functions `grpselect` and
`r2select`. Users can call the main function `RAEN` to run the whole
process,

``` r
library(RAEN)
myres<-RAEN(x,y, B = 50, ncore=3)
```

where `x` is the `n`x`p` predictor matrix, `y` is the time and
censoring status data frame, and `ncore` is the number of threads to use
for parallel processing. The selected variables and the regression
coefficients are returned.

``` 
  id        coef
 x1 -0.21695698
 x2 -0.41486949
 x3 -0.21307438
 x4 -0.05911834
...
...
 x75 -0.3823193
 x76 -0.6026619
 x77 -0.2440929
 x78 -0.3976471
 x80 -0.3249638
```
