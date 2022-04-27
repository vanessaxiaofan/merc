library(mgcv)
library(coxed)
## code to prepare `DATASET` dataset goes here
set.seed(1000)
#### Logistic Regression on one covariates measured with error and one confounder， 4 replicates ####
x<-rnorm(3000,0,1)
#ICC=0.7 generate
z.main <- matrix(x[1:1500]+rnorm(1500,0,sqrt(0.4)))
W<-matrix(sapply(x[1:1500], function(t){if(t>median(x)) {return(rbinom(1,1,0.5))}
  if(t<=median(x)){return(rbinom(1,1,0.3))}}))
r<-c(rep(3,700),rep(4,800))
z.rep<- cbind(x[1501:2200]+rnorm(700,0,sqrt(0.4)),
                   x[1501:2200]+rnorm(700,0,sqrt(0.4)),
                   x[1501:2200]+rnorm(700,0,sqrt(0.4)),
                   x[1501:2200]+rnorm(700,0,sqrt(0.4)))

#prevalence about 0.105
p<-exp(-2.4+log(1.5)*x[1:1500]+log(1.5)*W)/
  (1+exp(-2.4+log(1.5)*x[1:1500]+log(1.5)*W))
Y<-sapply(p,function(x) rbinom(1,1,x))

mainLogRel <- as.data.frame(cbind(z.main,W,Y))
colnames(mainLogRel) <- c("x","s", "Y")
relibLog <- as.data.frame(z.rep)
colnames(relibLog) <- c("x1","x2","x3","x4")

#### Linear Regression on two covariates measured with error and two confounders , 2 replicates####

x<-rmvn(3000,c(0,0,0),matrix(c(1,0.3,0.2,0.3,1,0.5,0.2,0.5,1),nrow=3))
w2<-sapply(x[,1], function(t){if(t>median(x[,1])) {return(rbinom(1,1,0.5))}
  if(t<=median(x[,1])){return(rbinom(1,1,0.3))}})
#ICC=0.7 generate z
r<-c(rep(2,1500))
W<-cbind(x[1:1500,3],w2[1:1500])

z.main = x[1:1500,1:2]+rnorm(1500,0,sqrt(0.4))

z.rep<-list(rbind(cbind(x[1501:2000,1]+rnorm(500,0,sqrt(0.4)),
                        x[1501:2000,1]+rnorm(500,0,sqrt(0.4))),
                  cbind(x[2001:2400,1]+rnorm(400,0,sqrt(0.4)),
                        x[2001:2400,1]+rnorm(400,0,sqrt(0.4)))),
            rbind(cbind(x[1501:2000,2]+rnorm(500,0,sqrt(0.4)),
                        x[1501:2000,2]+rnorm(500,0,sqrt(0.4))),
                  cbind(x[2001:2400,2]+rnorm(400,0,sqrt(0.4)),
                        x[2001:2400,2]+rnorm(400,0,sqrt(0.4)))))
#prevalence about 0.105
Y <- -2.7+1.5*rowSums(x[1:1500,1:3])+1.5*w2[1:1500]

mainLinearRel <- as.data.frame(cbind(z.main,W,Y))
colnames(mainLinearRel) <- c("x","z","s1","s2", "Y")
t1 <- as.data.frame(z.rep[[1]])
t2 <- as.data.frame(z.rep[[2]])
relibLinear <- as.data.frame(cbind(t1$V1,t2$V1,t1$V2,t2$V2))
colnames(relibLinear) <- c("x1","z1","x2","z2")

#### Cox Regression on one covariate measured with error and two confounders , 4 replicates####
simdata <- coxed::sim.survdata(N=1500, T=250, xvars=3, censor=.2, num.data.frames = 1)
#summary(simdata$data)

mainCoxRel <- as.data.frame(simdata$data)
names(mainCoxRel) <- c("x","s1","s2","Y","failed")

x <- mainCoxRel$x

z.rep<-rbind(cbind(x[1:750]+rnorm(750,0,sqrt(0.4)),
                   x[1:750]+rnorm(750,0,sqrt(0.4)),
                   x[1:750]+rnorm(750,0,sqrt(0.4)),
                   x[1:750]+rnorm(750,0,sqrt(0.4))))
relibCox <- as.data.frame(z.rep)
colnames(relibCox) <- c("x1","x2","x3","x4")


usethis::use_data(mainLogRel,overwrite=TRUE)
usethis::use_data(mainLinearRel,overwrite=TRUE)
usethis::use_data(mainCoxRel,overwrite=TRUE)
usethis::use_data(relibLog,overwrite=TRUE)
usethis::use_data(relibLinear,overwrite=TRUE)
usethis::use_data(relibCox,overwrite=TRUE)
