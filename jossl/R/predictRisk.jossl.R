### predictRisk.jossl.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  5 2025 (16:09) 
## Version: 
## Last-Updated: Jul  8 2025 (08:30) 
##           By: Thomas Alexander Gerds
##     Update #: 45
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Predict method for joint survival super learner
##'
##' Predicts the risks of all causes and also the censoring distribution.
##' @title Predict method for joint survival super learner
##' @param object Object obtained with \code{jossl} 
##' @param newdata Newdata
##' @param times Times
##' @param cause Cause
##' @param ... Not used
##' @return matrix with one column for each value in \code{times} and
##' one row for each individual in \code{newdata}
##' @seealso \link[riskRegression]{predictRisk}
##' @examples
##' set.seed(8)
##' library(survival)
##' d <- riskRegression::sampleData(87)
##' lol <- list(
##'   cox_lasso = list("cox", x_form = ~X2+X7),
##'   aalen_johansen =  list("cox", x_form = ~ 1))
##' x = jossl(list(cause1 = lol,cause2 = lol,censor = lol),
##'           data = d,
##'           time = 36,
##'           integrate = TRUE,
##'           verbose = TRUE,
##'           split.method = "cv5",
##'           B = 1,
##'           time_name = "time",
##'           status_name = "event",
##'           cause_codes = c("1" = 1, "2" = 2, "c" = 0),
##'           vars = c("X1","X2","X3","X4","X5","X6","X7"))
##' predictRisk.jossl(x,times = 3:8,newdata = d[1:3],cause = 1)
##' @method predictRisk jossl
##' @export predictRisk.jossl
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
predictRisk.jossl <- function(object,
                              newdata,
                              times,
                              cause,
                              ...){
    if(missing(times)){
        time = object$times
    }else{
        stopifnot(all(times<max(object$times)))
    }
    ## etimes <- object$times[object$times <= max(times)]
    etimes <- sort(unique(object$times))
    etimes <- etimes[etimes <= max(times)]
    if (cause == 0){
        Gamma <- predictCHF(object$jossl$censor,newdata = newdata,times = etimes,...)
        risk <- 1-exp(-Gamma)
    }else{
        A1 <- predictCHF(object$jossl$cause1,newdata = newdata,times = etimes,...)
        A2 <- predictCHF(object$jossl$cause2,newdata = newdata,times = etimes,...)
        NT = NCOL(A1)
        if (cause == 1){
            risk <- riskRegression::rowCumSum(cbind(1,exp(-(A1+A2))[,-NT,drop = FALSE]) * t(apply(cbind(0,A1),1,diff)))
        }else{
            risk <- riskRegression::rowCumSum(cbind(1,exp(-(A1+A2))[,-NT,drop = FALSE]) * t(apply(cbind(0,A2),1,diff)))
        }
    }
    out <- cbind(0,risk)[,prodlim::sindex(eval.times = times,jump.times = etimes)+1,drop = FALSE]
    # enforce the probability range
    out <- pmax(pmin(out,1),0)
    out
}


######################################################################
### predictRisk.jossl.R ends here
