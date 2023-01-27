
# merc

<!-- badges: start -->

[![R-CMD-check](https://github.com/vanessaxiaofan/merc/workflows/R-CMD-check/badge.svg)](https://github.com/vanessaxiaofan/merc/actions)
<!-- badges: end -->

The goal of merc is to implement regression calibration methods for
validation study and reliability study.

## Installation

You can install the released version of merc from the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vanessaxiaofan/merc")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(merc)
## basic example code
# # Only one mismeasured covariate, logistic model
y <- c(1,1,2,3,2,3)
x <- c(2,2.1,3,4,2.9,4.1)
s <- c(2,5,4,3,3,5)
case <- c(0,1,0,1,0,0)
age <- c(10,10,11,12,13,13)
test <- data.frame(y,x,s,case,age)

x <- c(1.1,1.2,0.8)
x2 <- c(1.2,1.2,0.9)
x3 <- c(0.9,1.0,1.3)
relib <- data.frame(x,x2,x3)

wts <- data.frame(x=1,s=1)


### Fit main study logistic model
outModel <- glm(case ~ x + s , family = binomial(link="logit"), data = test)
outcomeParam=coef(outModel)
outcomeParamVCOV=vcov(outModel)
outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
Bstar<-outcomeParam[2:length(outcomeParam)] #p' x 1
VBstar<-outcomeParamVCOV[2:length(outcomeParam),2:length(outcomeParam)] # p' x p'

# Provide point estimate
fit1 <- mercRel(supplyEstimates=TRUE, relib=relib, pointEstimates = Bstar, vcovEstimates = VBstar, sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), rr=3, ms=test, weights=wts,link = "logit",method = "glm" )
fit1
#> 
#> Call:
#> mercRel(supplyEstimates = TRUE, relib = relib, pointEstimates = Bstar, 
#>     vcovEstimates = VBstar, sur = c("x"), woe = c("s"), weri = c("x", 
#>         "x2", "x3"), rr = 3, ms = test, weights = wts, method = "glm", 
#>     link = "logit")
#> 
#> Coefficients Uncorrected Model:
#>   Weights           B   SE(B)   OR(B)     Z Value  Pr(>|Z|) lower 95%CI
#> x       1 -0.04791056 1.10300 0.95322 -0.04343659 0.9653535     0.10972
#> s       1  0.42994470 0.83453 1.53717  0.51519382 0.6064176     0.29947
#>   upper 95%CI
#> x     8.28105
#> s     7.89022
#> 
#> Coefficients Corrected Model:
#>   Weights        B   SE(B)   OR(B)     Z Value  Pr(>|Z|) lower 95%CI(OR)
#> x       1 -0.05027 1.15735 0.95097 -0.04343543 0.9653544         0.09840
#> s       1  0.43037 0.83639 1.53783  0.51455661 0.6068629         0.29851
#>   lower 95%CI(OR)
#> x         9.19018
#> s         7.92240

# Without point estimate
fit2 <- mercRel(supplyEstimates=FALSE, relib = relib, sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), outcome = c("case"), rr=3, ms=test,method = "glm", family = binomial, link = "logit", weights=wts)
fit2
#> 
#> Call:
#> mercRel(supplyEstimates = FALSE, relib = relib, sur = c("x"), 
#>     woe = c("s"), outcome = c("case"), weri = c("x", "x2", "x3"), 
#>     rr = 3, ms = test, weights = wts, method = "glm", family = binomial, 
#>     link = "logit")
#> 
#> Coefficients Uncorrected Model:
#>   Weights           B   SE(B)   OR(B)     Z Value  Pr(>|Z|) lower 95%CI
#> x       1 -0.04791056 1.10300 0.95322 -0.04343659 0.9653535     0.10972
#> s       1  0.42994470 0.83453 1.53717  0.51519382 0.6064176     0.29947
#>   upper 95%CI
#> x     8.28105
#> s     7.89022
#> 
#> Coefficients Corrected Model:
#>   Weights        B   SE(B)   OR(B)     Z Value  Pr(>|Z|) lower 95%CI(OR)
#> x       1 -0.05027 1.15735 0.95097 -0.04343543 0.9653544         0.09840
#> s       1  0.43037 0.83639 1.53783  0.51455661 0.6068629         0.29851
#>   lower 95%CI(OR)
#> x         9.19018
#> s         7.92240
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.
