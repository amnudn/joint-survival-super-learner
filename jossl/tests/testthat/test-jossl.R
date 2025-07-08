### test-jossl.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  7 2025 (10:35) 
## Version: 
## Last-Updated: Jul  7 2025 (12:20) 
##           By: Thomas Alexander Gerds
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(riskRegression)
testthat::test_that("comparison with predictRisk.CauseSpecificCox",{
    set.seed(8)
    d <- sampleData(1134)
    nd <- data.table(X7 = 47,X1 = factor(1,levels = 0:1))
    ttt <- 0:10
    F <- CSC(Hist(time,event)~X1+X7,data = d,fitter = "glmnet")
    f1 <- GLMnet(Surv(time,event == 1)~X1+X7,data = d)
    f2 <- GLMnet(Surv(time,event == 2)~X1+X7,data = d)
    expect_equal(c(predictRisk(F,cause = 2,times = ttt,newdata = nd)),
                 c(predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd,cause = 2,times = ttt)),
                 tolerance = 0.001)
})
testthat::test_that("jossl runs and predicts",{
    set.seed(8)
    d <- sampleData(1134)
    lol <- list(
        cox_lasso = list("cox", x_form = ~X1+X2+X3+X4+X5+X6+X7),
        aalen_johansen =  list("cox", x_form = ~ 1)
    )
    Build(jossl,Source = "~/research/SuperVision/Anders/statelearner/")
    x = jossl(list(cause1 = lol,cause2 = lol,censor = lol),
              data = d,
              time = 36,
              integrate = TRUE,
              verbose = TRUE,
              split.method = "cv5",
              B = 1,
              time_name = "time",
              status_name = "event",
              cause_codes = c("1" = 1, "2" = 2, "c" = 0),
              vars = c("X1","X2","X3","X4","X5","X6","X7"))
    predictRisk(x,times = 10,newdata = d[1:3],cause = 1)
})


######################################################################
### test-jossl.R ends here
