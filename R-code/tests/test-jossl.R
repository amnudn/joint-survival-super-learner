### test-jossl.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  7 2025 (08:28) 
## Version: 
## Last-Updated: Jul  7 2025 (10:14) 
##           By: Thomas Alexander Gerds
##     Update #: 12
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
library(randomForestSRC)
tar_source("~/research/SuperVision/Anders/statelearner/R-code/functions/")

d <- sampleData(1134)
nd <- data.table(X7 = 50)
ttt <- 10
F <- CSC(Hist(time,event)~X7,data = d)
f1 <- coxph(Surv(time,event == 1)~X7,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X7,data = d,x = TRUE,y = TRUE)
c(predictRisk(F,cause = 1,times = ttt,newdata = nd),
  predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd,cause = 1,times = ttt))

d <- sampleData(1134)
nd <- data.table(X1 = factor(1,levels = 0:1))
ttt <- 10
F <- CSC(Hist(time,event)~X1,data = d)
f1 <- coxph(Surv(time,event == 1)~X1,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X1,data = d,x = TRUE,y = TRUE)
c(predictRisk(F,cause = 1,times = ttt,newdata = nd),
  predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd,cause = 1,times = ttt))


d <- sampleData(1134)
nd <- data.table(X7 = 47,X1 = factor(1,levels = 0:1))
ttt <- 0:10
F <- CSC(Hist(time,event)~X1+X7,data = d)
f1 <- coxph(Surv(time,event == 1)~X1+X7,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X1+X7,data = d,x = TRUE,y = TRUE)
cbind(c(predictRisk(F,cause = 2,times = ttt,newdata = nd)),
      c(predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd,cause = 2,times = ttt)))

d <- sampleData(1134)
nd <- data.table(X7 = 47,X1 = factor(1,levels = 0:1))
ttt <- 0:10
F <- CSC(Hist(time,event)~X1+X7,data = d,fitter = "glmnet")
f1 <- GLMnet(Surv(time,event == 1)~X1+X7,data = d)
f2 <- GLMnet(Surv(time,event == 2)~X1+X7,data = d)
cbind(c(predictRisk(F,cause = 2,times = ttt,newdata = nd)),
      c(predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd,cause = 2,times = ttt)))




set.seed(9)
d <- sampleData(34)
setorder(d,time,-event)
nd <- data.table(X1 = 1,X7 = 50)
ttt <- 10
F <- CSC(Hist(time,event)~X7,data = d)
f1 <- coxph(Surv(time,event == 1)~X7,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X7,data = d,x = TRUE,y = TRUE)
F1 <- predictRisk(F,cause = 1,times = ttt,newdata = nd)
uu <- readRDS("~/tmp/uu.rds")

j1 <- predictRisk.jossl(list(times = ttt,cause1 = f1,cause2 = f2),newdata = nd,cause = 1,centered = FALSE)
J1 <- predictRisk.jossl(list(times = ttt,cause1 = f1,cause2 = f2),newdata = nd,cause = 1)
plot(J1,F1,xlim = c(0,1),ylim = c(0,1));abline(0,1,col = 2)



set.seed(9)
d <- sampleData(113)
nd <- data.table(X1 = 1,X7 = 50)
F <- CSC(Hist(time,event)~X7,data = d)
predictRisk(F,newdata = nd,times = 5)


f1 <- coxph(Surv(time,event == 1)~X1,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X1,data = d,x = TRUE,y = TRUE)
J1 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 1)
J2 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 2)
F <- CSC(Hist(time,event)~X1,data = d)
F1 <- predictRisk(F,times = sort(unique(d$time)),cause = 1,prodlim = FALSE,newdata = nd[1,])
F2 <- predictRisk(F,times = sort(unique(d$time)),cause = 2,prodlim = FALSE,newdata = nd[1,])
cbind(c(J1),c(F1))
plot(J1,F1,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)
cbind(c(J2),c(F2))
plot(J2,F2,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)

set.seed(0)
# calculating risk in uncensored data
d <- sampleData(1134)[,.(time,event,X7)]



nd <- data.frame(X1 = factor(1,levels = 0:1),X2 = factor(1,levels = 0:1),X7 = 50)
aj <- prodlim(Hist(time,event)~X7,data = d)
R1 <- predictRisk(aj,cause = 1,newdata = nd[1,],times = sort(unique(d$time)))
f1 <- coxph(Surv(time,event == 1)~X7,data = d,x = TRUE,y = TRUE)
f2 <- coxph(Surv(time,event == 2)~X7,data = d,x = TRUE,y = TRUE)
J1 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 1)
J2 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 2)
F <- CSC(Hist(time,event)~X1,data = d)
F1 <- predictRisk(F,times = sort(unique(d$time)),cause = 1,prodlim = FALSE,newdata = nd[1,])
F2 <- predictRisk(F,times = sort(unique(d$time)),cause = 2,prodlim = FALSE,newdata = nd[1,])
cbind(R1,c(J1),c(F1))
plot(J1,F1,xlim = c(0,1),ylim = c(0,1))
plot(R1,J1,xlim = c(0,1),ylim = c(0,1))
plot(R1,F1,xlim = c(0,1),ylim = c(0,1))

plot(aj1,F1,xlim = c(0,1),ylim = c(0,1))

abline(0,1,col = 2)
cbind(c(J2),c(F2))
plot(J2,F2,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)

f1 <- GLMnet(Surv(time,event == 1)~X1+X2+X7,data = d)
f2 <- GLMnet(Surv(time,event == 2)~X1+X2+X7,data = d)
nd <- data.frame(X1 = 1,X2 = 1,X7 = 50)
J1 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 1)
J2 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 2)
F <- CSC(Hist(time,event)~X1,data = d)
F1 <- predictRisk(F,times = sort(unique(d$time)),cause = 1,prodlim = FALSE,newdata = nd[1,])
F2 <- predictRisk(F,times = sort(unique(d$time)),cause = 2,prodlim = FALSE,newdata = nd[1,])
cbind(c(J1),c(F1))
plot(J1,F1,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)
cbind(c(J2),c(F2))
plot(J2,F2,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)

f1 <- GLMnet(Surv(time,event == 1)~X1+X2+X7,data = d,x = TRUE,y = TRUE)
f2 <- GLMnet(Surv(time,event == 2)~X1+X2+X7,data = d,x = TRUE,y = TRUE)
nd <- data.frame(X1 = 1,X2 = 1,X7 = 50)
J1 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 1)
J2 <- predictRisk.jossl(list(times = sort(unique(d$time)),cause1 = f1,cause2 = f2),newdata = nd[1,],cause = 2)
F <- CSC(Hist(time,event)~X1,data = d)
F1 <- predictRisk(F,times = sort(unique(d$time)),cause = 1,prodlim = FALSE,newdata = nd[1,])
F2 <- predictRisk(F,times = sort(unique(d$time)),cause = 2,prodlim = FALSE,newdata = nd[1,])
cbind(c(J1),c(F1))
plot(J1,F1,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)
cbind(c(J2),c(F2))
plot(J2,F2,xlim = c(0,1),ylim = c(0,1))
abline(0,1,col = 2)


######################################################################
### test-jossl.R ends here
