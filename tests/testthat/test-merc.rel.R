# # Only one mismeasured covariate, logistic model
y <- c(1,1,2,3,2,3)
x <- c(2,2.1,3,4,2.9,4.1)
s <- c(2,5,4,3,3,5)
case <- c(0,1,0,1,0,0)
age <- c(10,10,11,12,13,13)
test <- data.frame(y,x,s,case,age)
test

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

fit1 <- merc.rel(supplyEstimates=TRUE, relib=relib, pointEstimates = Bstar, vcovEstimates = VBstar, sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), rr=3, ms=test, weights=wts,link = "logit",method = "glm" )
fit2 <- merc.rel(supplyEstimates=FALSE, relib = relib, sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), outcome = c("case"), rr=3, ms=test,method = "glm", family = binomial, link = "logit", weights=wts)


test_that("merc.rel function works", {
  expect_equal(merc.rel(supplyEstimates=TRUE, relib=relib, pointEstimates = Bstar, vcovEstimates = VBstar, sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), rr=3, ms=test, weights=wts,link = "logit",method = "glm" ), fit1)
  expect_equal(merc.rel(supplyEstimates=FALSE, relib = relib, sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), outcome = c("case"), rr=3, ms=test,method = "glm", family = binomial, link = "logit", weights=wts), fit2)
})



