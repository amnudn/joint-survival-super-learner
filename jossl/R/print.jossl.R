### print.jossl.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  7 2025 (14:47) 
## Version: 
## Last-Updated: Jul  8 2025 (08:45) 
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
##' Print method for jossl
##'
##' @title Print jossl object
##' @param x Object obtained with \code{Score.list}
##' @param digits Number of digits
##' @param ... passed to print
#'
#' @method print jossl
#' @export print.jossl
#' @export
print.jossl <- function(x,...){
    name_winners <- as.character(unlist(x$cv_fit[1,c("cause1","cause2","censor"),with = FALSE]))
    cat("Joint survival super learner fitted with call:\n")
    print(x$call)
    cat("\nThe discrete jossl selects:\n")
    cat("\nCumulative hazard cause 1:",name_winners[[1]])
    cat("\nCumulative hazard cause 2:",name_winners[[2]])
    cat("\nCumulative hazard censor:",name_winners[[3]])
    cat("\n\nUse predictRisk(x) to obtain predicted values.\n")
}


######################################################################
### print.jossl.R ends here
