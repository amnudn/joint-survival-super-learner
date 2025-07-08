### predictCHF.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Jul  5 2025 (16:11) 
## Version: 
## Last-Updated: Jul  7 2025 (08:30) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## Cumulative hazard functions
predictCHF <- function(object, newdata, ...){
    UseMethod("predictCHF",object)
}
predictCHF.coxph <- function(object, newdata, times, type = "cumhazard",...){
    pred_obj = riskRegression::predictCox(object, newdata, times = times, type = type,...)
    chf = pred_obj$cumhazard
    return(chf)
}
predictCHF.GLMnet <- function(object, newdata, times,...){
    ## Hackidy-hack...
    risk_pred = riskRegression:::predictRisk.GLMnet(object = object, newdata = newdata, times = times, product.limit = 0, ...)
    chf = -log(1-risk_pred)
    return(chf)
}
predictCHF.rfsrc <- function(object, newdata, times, ...){
    ## Unsafe hack...
    if("time" %in% names(newdata))
        wd = copy(newdata)[, -c("time"), with = FALSE]
    else
        wd = copy(newdata)
    jump_times = object$time.interest
    cschf = predict(object, newdata = wd)$chf[, , 1]
    if(jump_times[1] != 0){
        jump_times = c(0, jump_times)
        cschf = cbind(0, cschf)
    }
    ii = prodlim::sindex(jump.times = jump_times,eval.times = times)
    out = matrix(cschf[, ii], ncol = length(times))
    ## Need to fix problem if prediction for one observation!
    return(out)
}


######################################################################
### predictCHF.R ends here
