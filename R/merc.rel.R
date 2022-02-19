#####################
# merc.rel function #
#####################

#' @title merc.rel
#' @author Xiaofan Liu and Xin Zhou
#' @description This funciton calculates regression coefficients, their standard
#'   errors, and odds ratios, when relevant, and 95% confidence intervals for a
#'   biologically meaningful difference specified by model covariates. Logistic
#'   model are implemented. A reliability study is required to empirically
#'   charaterize the measurement error model. Details are given in Rosner et
#'   al.(1989), Rosner et al.(1990), and Rosner et al.(1992) including "real
#'   data" examples
#' @references Rosner B, Spiegelman D, Willett WC. Correction of logistic
#'   regression relative risk estimates and confidence intervals for measurement
#'   error: the case of multiple covariates measured with error. Am J Epidemiol.
#'   1990 Oct;132(4):734-45. doi: 10.1093/oxfordjournals.aje.a115715. PMID:
#'   2403114.
#' @param supplyEstimates Indicates whether uncorrected estimates will be
#'   supplied by the user.
#' @param relib Name of the reliability dataset
#' @param ms name of main study data set
#' @param pointEstimates A numeric vector of point estimates from uncorrected
#'   standard regression analysis. It must be in the order of the variables
#'   indicated in sur and woe Intercept estimate must be removed. Only binary
#'   and numeric variables are accepted, all other classes of variables must be
#'   transformed (this is automatically done by most regression softwares). Must
#'   include names for each point estimates corresponding to the (transformed)
#'   names from `covCalib` followed by `covOutcome`. Must be supplied if
#'   supplyEstimates = TRUE.
#' @param vcovEstimates A p by p Variance-covariance matrix estimates from
#'   uncorrected standard regression analysis. Intercept estimates must be
#'   removed. Only binary and numeric variables are accepted, all other classes
#'   of variables must be transformed (this is automatically done by most
#'   regression softwares). Must include column names (excluding intercept) for
#'   the estimates corresponding to the (transformed) names from `covCalib`
#'   followed by `covOutcome`. Must be supplied if supplyEstimates = TRUE.
#' @param sur character vector of mismeasured exposure and covariates (i.e.
#'   surrogates) in the main study dataset
#' @param woe Character vector of names of perfectly measured covariates
#' @param outcome Outcome variable
#' @param weri  Character vector of names of variables which have reliability
#'   measures. These should have a set of variable names for each reliability
#'   measure (e.g. x1 y1 z1 x2 y2 z2 for measures). The variables must be in the
#'   same order as the names in the names and increments dataset specified in
#'   `pointEstimates`
#' @param rr Number of reliability measures taken (e.g. 2).
#' @param weights Name of the dataset containing the variable names and
#'   increments for the odds ratio or regression slopes.Including mismeasured
#'   and perfectly measured variables
#' @param r1 Number of replicates in main study
#' @param method Methods for modeling, currently only `lm` or `glm` methods are
#'   available. Required.
#' @param family Supply family parameter to pass to glm function. Not a
#'   character. Required if method="glm".
#' @param link Supply link parameter to pass to glm function. Should be
#'   character. Required if method="glm".
#' @return printable dataframe from standard regression results (when
#'   supplyEstimates==FALSE) as well as corrected results
#' @examples
#' # # Only one mismeasured covariate, linear model
#' y <- c(1,1,2,3,2,3)
#' x <- c(2,2.1,3,4,2.9,4.1)
#' s <- c(2,5,4,3,3,5)
#' case <- c(0,1,0,1,0,0)
#' age <- c(10,10,11,12,13,13)
#' test <- data.frame(y,x,s,case,age)
#'
#' x <- c(1.1,1.2,0.8)
#' x2 <- c(1.2,1.2,0.9)
#' x3 <- c(0.9,1.0,1.3)
#' relib <- data.frame(x,x2,x3)
#' wts <- data.frame(x=1,s=1)
#'
#' ### Fit main study linear model
#' outModel <- lm(y ~ x + s , data = test)
#' outcomeParam=coef(outModel)
#' outcomeParamVCOV=vcov(outModel)
#' outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
#' Bstar<-outcomeParam[2:length(outcomeParam)] #p' x 1
#' VBstar<-outcomeParamVCOV[2:length(outcomeParam),2:length(outcomeParam)]
#' merc.rel(supplyEstimates=TRUE, relib=relib, pointEstimates = Bstar,
#'          vcovEstimates = VBstar,sur = c("x"), woe = c("s"), weri = c("x","x2","x3"),
#'          rr=3, ms=test, weights=wts, method = "lm" )
#' merc.rel(supplyEstimates=FALSE, relib = relib, sur = c("x"), woe = c("s"),
#'          weri = c("x","x2","x3"), outcome = c("y"), rr=3,
#'          ms=test,method = "lm",weights=wts)
#'
#' # # Only one mismeasured covariate, logistic model
#' ### Fit main study logistic model
#' outModel <- glm(case ~ x + s , family = binomial(link="logit"), data = test)
#' outcomeParam=coef(outModel)
#' outcomeParamVCOV=vcov(outModel)
#' outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
#' Bstar<-outcomeParam[2:length(outcomeParam)] #p' x 1
#' VBstar<-outcomeParamVCOV[2:length(outcomeParam),2:length(outcomeParam)] # p' x p'
#' merc.rel(supplyEstimates=TRUE, relib=relib, pointEstimates = Bstar, vcovEstimates = VBstar,
#'          sur = c("x"), woe = c("s"), weri = c("x","x2","x3"), rr=3, ms=test, weights=wts,
#'          link = "logit",method = "glm" )
#' merc.rel(supplyEstimates=FALSE, relib = relib, sur = c("x"), woe = c("s"),
#'          weri = c("x","x2","x3"),outcome = c("case"), rr=3, ms=test,
#'          method = "glm", family = binomial, link = "logit", weights=wts)
#' @importFrom stats as.formula coef cov na.omit vcov
#' @export

merc.rel <- function(supplyEstimates=FALSE, relib, pointEstimates=NA, vcovEstimates=NA, sur, woe=NA,
                     outcome=NA, weri, rr, ms, weights, r1=1, method="lm",family=NA, link=NA){
  # require("Matrix")
  # require("gdata")
  # require("dplyr")
  ###################
  # check arguments #
  ###################
  if(missing(relib)){
    stop("Input reliability data is not supplied.")
  }else if(class(relib)!="data.frame"){
    stop("Input reliability data must be of data.frame class.")
  }

  if(supplyEstimates==FALSE){
    if(missing(ms)){
      stop("Input data ms not supplied.")
    }else if(class(ms)!="data.frame"){
      stop("Input data ms must be of data.frame class.")
    }
  }

  if(missing(sur)){
    stop("Missing variables measured with error in main study dataset.")
  }else if(class(sur)!="character"){
    stop("Variables measured with error in main study dataset are not supplied with character vector.")
  }

  if(sum(is.na(woe))==0 & class(woe)!="character"){
    stop("Variables measured without error are not supplied with character vector.")
  }

  if(missing(weights)){
    stop("Missing weights for calculations")
  }else if(class(weights)!="data.frame"){
    stop("Input weights must be of data.frame class.")
  }

  if(supplyEstimates==FALSE){
    if(missing(outcome)){
      stop("Outcome is missing.")
    }else if(class(outcome)!="character"|outcome==""|outcome==" "){
      stop("outcome is not supplied with appropriate character.")
    }
  }

  ## check if main dataset contains data indicated by sur, woe
  if(supplyEstimates==FALSE){
    if(length(woe)>0){
      MSVars_spec<-sur
    }else if(length(woe)==0){
      MSVars_spec<-(c(sur,woe))
    }
    MSVars<-colnames(ms)
    inMSVars<-(MSVars_spec%in%MSVars)
    if(sum(inMSVars)!=length(MSVars_spec)){
      stop("Main study dataset does not contain all the necessary variables specified in one of the following parameter: id, sur, woe.")
    }
  }

  ## check if reliability dataset contains data indicated by weri
  relibVars <- colnames(relib)
  relibVars_spec <- weri
  inrelibVars <- relibVars_spec %in% relibVars
  if(sum(inrelibVars)!=length(relibVars_spec)){
    stop("Reliability study dataset does not contain all the necessary variables specified in weri.")
  }

  if(supplyEstimates==TRUE){
    ## check whether there are names for the point Estimates and column names for the vcov estimates
    if(length(names(pointEstimates))==0|length(colnames(vcovEstimates))==0){
      stop("There must be names for user supplied point estiamtes and column names for variance covariance matrix!")
    }

    ## check if vcov estiamtes are a square matrix
    if(dim(vcovEstimates)[1]!=dim(vcovEstimates)[2]){
      stop("Covariance matrix supplied is not a symmetric square matrix. Please check!")
    }
  }

  ##################################
  #conduct regression when supplyEstimates=FALSE #
  ##################################
  # create formula for outcome model
  if(supplyEstimates==FALSE){
      outcomeFormula<-paste0("~",paste0(sur,collapse="+"),"+",paste0(woe,collapse="+"))
      allVars_ms<-c(outcome,sur,woe)
  }

  if(supplyEstimates==FALSE){
    ms_complete<- na.omit(dplyr::select(ms,dplyr::all_of(allVars_ms)))
    X_MS <- stats::model.matrix(object= as.formula(outcomeFormula),data=ms_complete)
    Y_MS<- ms_complete[,outcome]
    outcomeModelVarNames<-colnames(X_MS)
  }else if(supplyEstimates==TRUE){
    outcomeModelVarNames<-c(names(pointEstimates))
  }

  #step 1: outcome model
  ## RSW outcome model modeling (additive)
  if(supplyEstimates==FALSE){

    if(method=="lm"){
      outModel<-stats::lm(formula=as.formula(paste0(outcome,outcomeFormula)),data=ms_complete)
      outcomeParam=coef(outModel)
      outcomeParamVCOV=vcov(outModel)
      outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
    }else{
      outModel<-stats::glm(formula=as.formula(paste0(outcome,outcomeFormula)),
                    data=ms_complete,family=family(link=link))
      outcomeParam=coef(outModel)
      outcomeParamVCOV=vcov(outModel)
      outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
    }
  }

  if(supplyEstimates==FALSE){
    B<-t(t(outcomeParam[2:length(outcomeParam)])) #p' x 1
    V<-outcomeParamVCOV[2:length(outcomeParam),2:length(outcomeParam)] # p' x p'
  }else if(supplyEstimates==TRUE){
    B<-as.matrix(pointEstimates)
    V<-vcovEstimates
  }

  ###############################
  # Caculate adjusted point estimates #
  ###############################

  WT <-  weights
  ###### Compute the corrected point estimate ######
  q <-  length(woe)  # number of without measurement error variables
  p <-  length(sur) # number of measurement error variables

  pq <- p+q # total numbers of varibles in main study
  pqsq <- pq*pq

  #Use the main study
  Z <-  ms[, c(sur,woe)]
  ## obtain sigmaZ
  SZ <- cov(Z) #checked

  #Use the reliability study
  r <-  rr # number of replicates
  S <-  matrix(data=0,nrow=p,ncol=p)
  SR <-  matrix(data=0,nrow=p,ncol=p)
  X <-  matrix(data=0,nrow=p,ncol=p)

  for(i in 1:nrow(relib)){
    X <-  matrix(data=as.numeric(relib[i,]),nrow=r,ncol=p,byrow = TRUE)
    SR <-  cov(X)/nrow(relib)
    S <-  S + SR
  }

  remove(list=(c("X","SR")))

  SE <-  as.matrix(Matrix::bdiag(S,matrix(0,nrow=q, ncol=q))) # q is the number of variables without measurement error
  SX <-  SZ-(SE/r1)
  SXV <-  solve(SX)
  RV <-  (SX + SE/r1) %*% SXV
  BLINR <-  t(B) %*% RV

  ###############################
  # Caculate adjusted variance estimate #
  ###############################

  tem1 <-  t(c(1:pq)%*%t(rep(1,pq)))
  gdata::upperTriangle(tem1) <- gdata::lowerTriangle(tem1, byrow=TRUE)
  tem2 <-  c(1:pq)%*%t(rep(1,pq))
  gdata::upperTriangle(tem2) <- gdata::lowerTriangle(tem2, byrow=TRUE)
  DESIGN1 <-  as.vector(tem1)
  DESIGN2 <-  as.vector(tem2)

  n <-  nrow(relib) # number of persons in reliability study
  n1 <-  nrow(ms) # number of persons in main study

  dnr <-  1/(n*(r-1))
  COVSE <-  (SE[DESIGN1,DESIGN1] * SE[DESIGN2,DESIGN2]
           + SE[DESIGN1,DESIGN2] * SE[DESIGN2,DESIGN1])*dnr
  #checked
  dk1 <- 1/(n1-1)
  COVSX <- (SZ[DESIGN1,DESIGN1] * SZ[DESIGN2,DESIGN2]
         +SZ[DESIGN1,DESIGN2] * SZ[DESIGN2,DESIGN1])*dk1
  COVSX  <-  COVSX + COVSE

  COVSXV <- diag(pqsq)
  COVSXE <- matrix(data=0,nrow=pqsq,ncol=pqsq)

  for (i in 1:pq){
    for(j in i:pq){
      rowSXVij <- kronecker(SXV[i,],SXV[,j])
      x <- j+(i-1)*pq
      x2 <- i+(j-1)*pq
      for(l in i:pq){
        for(k in l:pq){
          #sumv=sum(kronecker(t(rowSXVij),(kronecker(SXV[k,],SXV[,l]))) * COVSX)
          sumv <- sum((t(rowSXVij) %x% (SXV[k,] %x% SXV[,l])) * COVSX)
          y <- l+(k-1)*pq
          y2 <- k+(l-1)*pq
          COVSXV[x,y] <- sumv
          COVSXV[y,x] <- sumv
          COVSXV[x2,y] <- sumv
          COVSXV[y,x2] <- sumv
          COVSXV[x2,y2] <- sumv
          COVSXV[y2,x2] <- sumv
          COVSXV[x,y2] <- sumv
          COVSXV[y2,x] <- sumv

        }
      }
      for(k in 1:p){
        for(l in k:p){
          y <- l+(k-1)*pq
          y2 <- k+(l-1)*pq
          sume <- -sum(rowSXVij*COVSE[,y])
          COVSXE[x,y] <- sume
          COVSXE[x,y2] <- sume
          COVSXE[x2,y] <- sume
          COVSXE[x2,y2] <- sume

        }
      }
    }
  }


  #####
  COVR <- matrix(data=0,nrow=pqsq,ncol=pqsq)
  GN <- matrix(data=0,nrow=pq,ncol=pq)
  WR <- diag(pq)

  for(j in 1:pq){
    jstripes <- seq(j,(pq-1)*pq+j,pq)
    for(l in j:pq){
      lstripes <- seq(l,(pq-1)*pq+l,pq)
      for(i in 1:p){
        iband <- c(((i-1)*pq+1) : (i*pq))
        x <- j+(i-1)*pq
        for(k in max((j==l)*i,1):pq){
          kband <- c(((k-1)*pq+1):(k*pq))
          y <- l+(k-1)*pq

          sumr <- sum((SXV[,j] %x% t(SXV[l,])) * COVSE[iband,kband])+sum( (SE[,k] %x% t(SXV[j,])) * COVSXE[lstripes,iband])+sum(( SE[,i] %x% t(SXV[l,])) * COVSXE[jstripes,kband])+sum(( SE[,i] %x% t(SE[k,])) * COVSXV[jstripes,lstripes])
          COVR[x,y] <- sumr
          COVR[y,x] <- sumr

          GN[i,k] <- sumr
          GN[k,i] <- COVR[(k-1)*pq+j,(i-1)*pq+l]

        }
      }
      WR[j,l] <- t(B) %*% GN %*% (B) ####not sure
      WR[l,j] <- WR[j,l]
    }
  }
  VARA  <-  t(RV) %*% V %*% RV+WR/(r1*r1) # no r1

  ###############################
  # Compose output table #
  ###############################
  # logistic model
  if(method=="glm"&link=="logit"){

    SED <- round(sqrt(diag(V)),5)
    ODDR <- round(exp(WT*B),5)
    UB <- round(exp(WT*(B+1.96*SED)),5)
    LB <- round(exp(WT*(B-1.96*SED)),5)
    #CI <- c(paste0(LB,"-",UB))
    zValue <- B/SED
    pValue <- 2*pnorm(-abs(as.numeric(zValue)))

    BLINR  <-  round(t(BLINR),5)
    SEDN <- round(sqrt(diag(VARA)),5)
    ODDRN <- round(exp(WT*BLINR),5)
    UBN <- round(exp(WT*(BLINR+1.96*SEDN)),5)
    LBN <- round(exp(WT*(BLINR-1.96*SEDN)),5)
    zValueN <- BLINR/SEDN
    pValueN <- 2*pnorm(-abs(as.numeric(zValueN)))
    #CIN <- c(paste0(LBN,"-",UBN))
    # compose output tables
    Uncorrected <-  data.frame(t(WT), B, SED, t(ODDR),zValue, pValue, t(LB), t(UB))
    colnames(Uncorrected)<-c("Weights","B","SE(B)","OR(B)","Z Value","Pr(>|Z|)","lower 95%CI","upper 95%CI")
    #rownames(Uncorrected)<-outcomeModelVarNames
    Corrected <-  data.frame(t(WT), BLINR, SEDN, t(ODDRN), zValueN, pValueN, t(LBN), t(UBN))
    colnames(Corrected)<-c("Weights","B","SE(B)","OR(B)","Z Value","Pr(>|Z|)", "lower 95%CI(OR)","lower 95%CI(OR)")
    #rownames(Corrected)<-outcomeModelVarNames
    outputList<-list(Uncorrected,Corrected)
    names(outputList)<-c("Uncorrected","Corrected")

  }
  else{

    SED <- round(sqrt(diag(V)),5)
    UB <- round(WT*(B+1.96*SED),5)
    LB <- round(WT*(B-1.96*SED),5)
    zValue <- B/SED
    pValue <- 2*pnorm(-abs(as.numeric(zValue)))

    BLINR  <-  round(t(BLINR),5)
    SEDN <- round(sqrt(diag(VARA)),5)
    UBN <- round((WT*(BLINR+1.96*SEDN)),5)
    LBN <- round((WT*(BLINR-1.96*SEDN)),5)
    zValueN <- BLINR/SEDN
    pValueN <- 2*pnorm(-abs(as.numeric(zValueN)))
    # compose output tables
    Uncorrected <-  data.frame(t(WT), B, SED, zValue, pValue, t(LB), t(UB))
    colnames(Uncorrected)<-c("Weights","B","SE(B)","Z Value","Pr(>|Z|)","lower 95%CI","upper 95%CI")
    #rownames(Uncorrected)<-outcomeModelVarNames
    Corrected <-  data.frame(t(WT), BLINR, SEDN, zValueN, pValueN, t(LBN), t(UBN))
    colnames(Corrected)<-c("Weights","B","SE(B)","Z Value","Pr(>|Z|)","lower 95%CI","lower 95%CI")
    #rownames(Corrected)<-names(B)
    outputList<-list(Uncorrected,Corrected)
    names(outputList)<-c("Uncorrected","Corrected")

  }
  return(outputList)
}



