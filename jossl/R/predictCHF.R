### predictCHF.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Jul  5 2025 (16:11) 
## Version: 
## Last-Updated: Jul  7 2025 (15:40) 
##           By: Thomas Alexander Gerds
##     Update #: 11
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Predicting cumulative hazard functions
#'
#' Predicting cumulative hazard functions
#' @param object Object
#' @param newdata newdata 
#' @param times times
#' @param ... additional arguments
#' @export
predictCHF <- function(object, newdata, ...){ 
    UseMethod("predictCHF",object)
}
##' @export
#' @rdname predictCHF
#' @method predictCHF coxph
#'
predictCHF.coxph <- function(object,newdata,times,...){
    riskRegression::predictCox(object, newdata, times = times, type = "cumhazard",...)$cumhazard 
}
##' @export
#' @rdname predictCHF
#' @method predictCHF GLMnet
predictCHF.GLMnet <- function(object, newdata, times,...){ 
    riskRegression::predictCox(object,newdata,times = times,type = "cumhazard",...)$cumhazard 
}
##' @export
#' @rdname predictCHF
#' @method predictCHF rfsrc
predictCHF.rfsrc <- function(object, newdata, times, ...){
    if (any((them <- match(object$yvar.names,names(newdata),nomatch = 0))>0)){
        newdata[[object$yvar.names[[1]]]] <- NULL
        newdata[[object$yvar.names[[2]]]] <- NULL
    }
    risk_pred = riskRegression::predictRisk(object = object, newdata = newdata, times = times,cause = 1, ...)
    -log(1-risk_pred)
}


######################################################################
### predictCHF.R ends here
