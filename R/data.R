#' This is a simulated Reliability study.
#' Logistic Regression on one covariate measured with error and one confounder， 4 replicates [reliability study]
#'
#' This is the reliability dataset (relib1)
#' @format A data frame with 1500 rows and 4 reliability measures for the mismeasured variable x:
#' \describe{
#'   \item{x1}{first reliability measure for the covariate measured with error}
#'   \item{x2}{second reliability measure for the covariate measured with error}
#'   \item{x3}{third reliability measure for the covariate measured with error}
#'   \item{x4}{fourth reliability measure for the covariate measured with error}
#'   }
#' @examples
#' data("relib1", package = "merc")
"relib1"

#'This is the main study dataset (mainRelib1)
#' @format A data frame with 5000 rows and 14 variables:
#' \describe{
#'   \item{x}{mismeasured variable, continuous}
#'   \item{s}{confounder, binary, perfectly measured }
#'   \item{Y}{outcome, continuous}
#' }
#' @examples
#' data("mainRelib1", package = "merc")
"mainRelib1"

#' This is a simulated Reliability study.
#' Linear Regression on two covariates measured with error and two confounders , 2 replicates
#' This is the reliability dataset (relib2)
#' @format A data frame with 1500 rows and 2 reliability measures for each 2 mismeasured variables
#' \describe{
#'   \item{x1}{first reliability measure for the mismeasured x}
#'   \item{z1}{first reliability measure for the mismeasured z}
#'   \item{x2}{second reliability measure for the mismeasured x}
#'   \item{z2}{fourth reliability measure for the mismeasured z}
#'   }
#' @examples
#' data("relib2", package = "merc")
"relib2"

#'This is the main study dataset (mainRelib1)
#' @format A data frame with 5000 rows and 14 variables:
#' \describe{
#'   \item{x}{mismeasured variable x , continuous}
#'   \item{z}{mismeasured variable z , continuous}
#'   \item{s1}{confounder, continuous, perfectly measured }
#'   \item{s2}{confounder, binary, perfectly measured }
#'   \item{Y}{outcome, continuous}
#' }
#' @examples
#' data("mainRelib2", package = "merc")
"mainRelib2"
