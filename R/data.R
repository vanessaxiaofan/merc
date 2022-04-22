#' This is a simulated Reliability study.
#' Logistic Regression on one covariate measured with error and one confounder， 4 replicates [reliability study]
#'
#' This is the reliability dataset (relibLog)
#' @format A data frame with 1500 rows and 4 reliability measures for the mismeasured variable x:
#' \describe{
#'   \item{x1}{first reliability measure for the covariate measured with error}
#'   \item{x2}{second reliability measure for the covariate measured with error}
#'   \item{x3}{third reliability measure for the covariate measured with error}
#'   \item{x4}{fourth reliability measure for the covariate measured with error}
#'   }
#' @examples
#' data("relibLog", package = "merc")
"relibLog"

#'This is the main study dataset (mainLogRel) of logistic regression
#' @format A data frame with 1500 rows and 3 variables:
#' \describe{
#'   \item{x}{mismeasured variable, continuous}
#'   \item{s}{confounder, binary, perfectly measured }
#'   \item{Y}{outcome, continuous}
#' }
#' @examples
#' data("mainLogRel", package = "merc")
"mainLogRel"

#' This is a simulated Reliability study.
#' Linear Regression on two covariates measured with error and two confounders , 2 replicates
#' This is the reliability dataset (relibLinear)
#' @format A data frame with 1500 rows and 2 reliability measures for each 2 mismeasured variables
#' \describe{
#'   \item{x1}{first reliability measure for the mismeasured x}
#'   \item{z1}{first reliability measure for the mismeasured z}
#'   \item{x2}{second reliability measure for the mismeasured x}
#'   \item{z2}{second reliability measure for the mismeasured z}
#'   }
#' @examples
#' data("relibLinear", package = "merc")
"relibLinear"

#'This is the main study dataset (mainLinearRel)
#' @format A data frame with 1500 rows and 5 variables:
#' \describe{
#'   \item{x}{mismeasured variable x , continuous}
#'   \item{z}{mismeasured variable z , continuous}
#'   \item{s1}{confounder, continuous, perfectly measured }
#'   \item{s2}{confounder, binary, perfectly measured }
#'   \item{Y}{outcome, continuous}
#' }
#' @examples
#' data("mainLinearRel", package = "merc")
"mainLinearRel"


#' This is a simulated Reliability study.
#' Cox Regression on one covariate measured with error and two confounders , 4 replicates
#' This is the reliability dataset (relibCox)
#' @format A data frame with 1500 rows and 4 reliability measures for one mismeasured variable
#' \describe{
#'   \item{x1}{first reliability measure for the covariate measured with error}
#'   \item{x2}{second reliability measure for the covariate measured with error}
#'   \item{x3}{third reliability measure for the covariate measured with error}
#'   \item{x4}{fourth reliability measure for the covariate measured with error}
#'   }
#' @examples
#' data("relibCox", package = "merc")
"relibCox"

#'This is the main study dataset (mainCoxRel)
#' @format A data frame with 1500 rows and 5 variables:
#' \describe{
#'   \item{x}{mismeasured variable x , continuous}
#'   \item{s1}{confounder, continuous, perfectly measured }
#'   \item{s2}{confounder, continuous,perfectly measured }
#'   \item{Y}{time, continuous}
#'   \item{failed}{event, binary}
#' }
#' @examples
#' data("mainCoxRel", package = "merc")
"mainCoxRel"



