library(targets)
library(tarchetypes)
if (!file.exists("~/research/SuperVision/Anders/statelearner/zelefsky-case-study")){
    library(here)
    library(parallel)
    try(setwd(here("zelefsky-case-study")),silent = TRUE)
    tar_source(here("R-code/functions"))
}else{
    try(setwd("~/research/SuperVision/Anders/statelearner/zelefsky-case-study/"),silent = TRUE)
    ## tar_source("../R-code/functions")
}

tar_option_set(packages = c("glmnet","lava","foreach","data.table","prodlim","survival","riskRegression","ranger","jossl","randomForestSRC"))

list(
    tar_target(zelefsky,{
        try(zelefsky <- setDT(get(load(here("data/zelefsky.rda")))),silent = TRUE)
        try(zelefsky <- setDT(get(load("data/zelefsky.rda"))),silent = TRUE)
        zelefsky$status=as.numeric(as.character(factor(zelefsky$recur,levels=c("No","Yes"),labels=c("0","1"))))
        zelefsky$logPSA <- log(zelefsky$psa)
        zelefsky$sDose <- as.vector(scale(zelefsky$dose))
        zelefsky$hormonesYes <- zelefsky$hormones=="Yes"
        for (v in levels(zelefsky$stage)){
            zelefsky[,paste("stage",v,sep="")] <- as.numeric(zelefsky$stage==v)
        }
        zelefsky[,event := (vital == "Dead")*2]
        zelefsky[status == 1,event := 1]
        z = zelefsky[, .(time,status,event,logPSA,stage,ggtot,sDose,hormones,vital)]
        z[]
    }),
    tar_target(random_split_inbag,{
        set.seed(17)
        sample(1:NROW(zelefsky),size = 0.632*NROW(zelefsky))
    }),
    tar_target(name = train_zelefsky,
               command = {
                   zelefsky[random_split_inbag]
               }),
    tar_target(name = test_zelefsky,
               command = {
                   zelefsky[setdiff(1:NROW(zelefsky),random_split_inbag)]
               }),
    tar_target(random_split_inbag_2,{
        set.seed(2)
        sample(1:NROW(zelefsky),size = 0.632*NROW(zelefsky))
    }),
    tar_target(name = train_zelefsky_2,
               command = {
                   zelefsky[random_split_inbag_2]
               }),
    tar_target(name = test_zelefsky_2,
               command = {
                   zelefsky[setdiff(1:NROW(zelefsky),random_split_inbag_2)]
               }),
    tar_target(name = jossl_zelefsky,
               command = {
                   zel_learner <- list(
                       cox_lasso = list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones),
                       cox_elnet = list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5),
                       cox_unpenalized = list("cox", x_form = ~logPSA+stage+ggtot+sDose+hormones),
                       aj =  list("cox", x_form = ~ 1),
                       rf = list("rfsrc", ntree = 500,x_form = ~logPSA+stage+ggtot+sDose+hormones)
                   )
                   fit <- jossl(learners = list(cause1 = zel_learner,cause2 = zel_learner,censor = zel_learner),
                                data = zelefsky,
                                time = 36,
                                integrate = TRUE,
                                split.method = "cv5",
                                B = 5,
                                verbose = FALSE,
                                time_name = "time",
                                status_name = "event",
                                cause_codes = c("1" = 1, "2" = 2, "c" = 0),
                                vars = c("logPSA","stage","ggtot","sDose","hormones"))
               }),
    tar_target(name = jossl_train,
               command = {
                   zel_learner <- list(
                       cox_lasso = list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones),
                       cox_elnet = list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5),
                       cox_unpenalized = list("cox", x_form = ~logPSA+stage+ggtot+sDose+hormones),
                       aj =  list("cox", x_form = ~ 1),
                       rf = list("rfsrc", ntree = 500,x_form = ~logPSA+stage+ggtot+sDose+hormones)
                   )
                   fit <- jossl(learners = list(cause1 = zel_learner,cause2 = zel_learner,censor = zel_learner),
                                data = train_zelefsky,
                                time = 36,
                                integrate = TRUE,
                                split.method = "cv5",
                                B = 5,
                                verbose = FALSE,
                                time_name = "time",
                                status_name = "event",
                                cause_codes = c("1" = 1, "2" = 2, "c" = 0),
                                vars = c("logPSA","stage","ggtot","sDose","hormones"))
               }),
    tar_target(name = jossl_predictions,{
        command = 
            predictRisk.jossl(jossl_train,newdata = test_zelefsky,times = 1:36,cause = 1)
    }),
    tar_target(name = score_jossl_1,
               command ={
                   x <- riskRegression::Score(list(Jossl = jossl_train,
                                                   CSC = CSC(Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,data = train_zelefsky)),
                                              times = 1:36,
                                              cause = 1,
                                              cens.model = "rfsrc",
                                              summary = "risks",
                                              plots = c("roc","cal"),
                                              formula = Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,
                                              data = test_zelefsky)
               }),
    tar_target(name = score_jossl_2,
               command ={
                   x <- riskRegression::Score(list(Jossl = jossl_train,
                                                   CSC = CSC(Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,data = train_zelefsky)),
                                              times = 1:36,
                                              cause = 2,
                                              cens.model = "rfsrc",
                                              summary = "risks",
                                              plots = c("roc","cal"),
                                              formula = Hist(time,event)~logPSA+stage+ggtot+sDose+hormones,
                                              data = test_zelefsky)
               }),
    tar_target(name = score_jossl_0,
               command ={
                   Gjossl = predictRisk(jossl_train,newdata = test_zelefsky,cause = 0,times = 1:36)
                   Gcox = coxph(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,data = train_zelefsky,x = TRUE,y = TRUE)
                   # Score does not work for cause=0 hence we change the roles of censored and event
                   x <- riskRegression::Score(list(Jossl = Gjossl,
                                                   Cox = Gcox),
                                              times = 1:36,
                                              cause = 1,
                                              cens.model = "rfsrc",
                                              plots = c("roc","cal"),
                                              summary = "risks",
                                              formula = Hist(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,
                                              data = test_zelefsky)
               })
)
