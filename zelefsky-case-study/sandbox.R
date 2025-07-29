### sandbox.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  5 2025 (09:54) 
## Version: 
## Last-Updated: Jul  7 2025 (21:39) 
##           By: Thomas Alexander Gerds
##     Update #: 20
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
library(survival)
library(riskRegression)
library(targets)
library(glmnet)
tar_source("~/research/SuperVision/Anders/statelearner/R-code/functions/")
d <- sampleData(34)
setorder(d,time,-event)
nd <- sampleData(3)
lol <- list(
    cox_lasso = list("GLMnet", x_form = ~X1+X2+X3+X4+X5+X6+X7),
    aalen_johansen =  list("cox", x_form = ~ 1)
)

f1 <- GLMnet(Surv(time,event == 1)~X1+X2+X3+X4+X5+X6+X7,data = d)
f1 <- coxph(Surv(time,event == 1)~X7,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X7,data = d,x = TRUE,y = TRUE)
A1 <- predictCox(f1,newdata = nd,c(0,sort(unique(d$time))),type = c("hazard","cumhazard"))
A2 <- predictCox(f2,newdata = nd,c(0,sort(unique(d$time))),type = c("hazard","cumhazard"))
data.table(time = A1$time,S = exp(-(A1$cumhazard+A2$cumhazard)),A1$hazard)

tmp <- list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2)
J1 <- predictRisk.jossl(tmp,newdata = nd[1,],cause = 1)
predictRisk.jossl(tmp,newdata = nd,cause = 2)

F <- CSC(Hist(time,event)~X7,data = d)
F1 <- predictRisk(F,times = sort(unique(d$time)),cause = 1,prodlim = TRUE,newdata = nd[1,])


p <- predictCox(f,newdata = nd,times = c(0,sort(unique(d$time))),type = c("hazard","cumhazard"))
p$hazard
p$cumhazard
all.equal(c(p$cumhazard),cumsum(p$hazard))
predictCHF(f,newdata = nd[1,],times = c(0,sort(unique(d$time))))

Build(jossl,Source = "~/research/SuperVision/Anders/statelearner/")
x = jossl(list(cause1 = lol,cause2 = lol,censor = lol),data = d,time = 36,integrate = TRUE,verbose = TRUE,split.method = "cv5",B = 1,time_name = "time",status_name = "event",cause_codes = c("1" = 1, "2" = 2, "c" = 0))

Lambda1 <- predictCHF(x$jossl$cause1,newdata = nd,times = 1:5)
Lambda2 <- predictCHF(x$jossl$cause2,newdata = nd,times = 1:5)

exp(-rowCumSum(Lambda1+Lambda2))*t(apply(Lambda1,1,function(x)diff(c(0,x))))

library(data.table)
library(jossl)
library(survival)
library(riskRegression)
library(targets)
library(randomForestSRC)
library(glmnet)
setwd("~/research/SuperVision/Anders/statelearner/zelefsky-case-study/")
## Build(jossl,Source = "~/research/SuperVision/Anders/statelearner/")
tar_load(zelefsky)
tar_load(train_zelefsky)
tar_load(test_zelefsky)
tar_load(jossl_train)
tar_load(jossl_predictions)
train_zelefsky[,cause1 := event == 1]
r <- rfsrc(Surv(time,cause1)~logPSA+stage+ggtot+sDose+hormones,data = train_zelefsky)
predictRisk(r,newdata = test_zelefsky,times = 36)

x <- riskRegression::Score(list(RF = rfsrc(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,data = zelefsky),
                                CSC = CSC(Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,data = zelefsky)),
                           times = 1:36,
                           cause = 1,
                           split.method = "cv5",
                           cens.model = "rfsrc",
                           formula = Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,
                           data = zelefsky)
plotBrier(x)

x <- riskRegression::Score(list(Jossl = jossl_predictions,
                                CSC(Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,data = train_zelefsky)),
                           times = 1:36,
                           cause = 1,
                           cens.model = "rfsrc",
                           formula = Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,
                           data = test_zelefsky)

x <- riskRegression::Score(list(Jossl = jossl_train,
                                CSC(Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,data = train_zelefsky)),
                           times = 1:36,
                           cause = 1,
                           split.method = "cv5",
                           cens.model = "rfsrc",
                           formula = Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,
                           data = zelefsky)




zelefsky[,event := (vital == "Dead")*2]
zelefsky[status == 1,event := 1]
zel_learner <- list(
    cox_lasso = list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones),
    cox_elnet = list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5),
    cox_unpenalized = list("cox", x_form = ~logPSA+stage+ggtot+sDose+hormones),
    aj =  list("cox", x_form = ~ 1),
    rf = list("rfsrc", ntree = 500)
)
set.seed(17)
inbag <- sample(1:NROW(zelefsky),size = 0.632*NROW(zelefsky))
train_zelefsky <- zelefsky[inbag]
test_zelefsky <- zelefsky[setdiff(1:NROW(zelefsky),inbag)]
fit <- jossl(learners = list(cause1 = zel_learner,cause2 = zel_learner,censor = zel_learner),
             data = train_zelefsky,
             time = 36,
             integrate = TRUE,
             split.method = "cv5",
             B = 1,
             verbose = FALSE,
             time_name = "time",
             status_name = "event",
             cause_codes = c("1" = 1, "2" = 2, "c" = 0),
             vars = c("logPSA","stage","ggtot","sDose","hormones"))
x <- riskRegression::Score(list(Jossl = fit),
                           times = 36,
                           formula = Hist(time,event)~1,
                           data = test_zelefsky)



######################################################################
### sandbox.R ends here


