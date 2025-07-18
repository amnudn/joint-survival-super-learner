### jossl.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Jul  5 2025 (16:09) 
## Version: 
## Last-Updated: Jul  8 2025 (08:39) 
##           By: Thomas Alexander Gerds
##     Update #: 36
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Joint survival super learner 
##'
##' Argument names may change in the future
##' @title Joint Survival Super Learner
##' @param learners library of learners
##' @param data data
##' @param time time 
##' @param integrate logical. If \code{TRUE} then integrate
##' @param split.method Passed on to riskRegression::getSplitMethod
##' @param collapse logical. If \code{TRUE} then average
##' @param B Passed on to riskRegression::getSplitMethod
##' @param verbose logical. If \code{TRUE} then talk a lot
##' @param time_grid_length Time grid to approximate the integrated Brier score
##' @param time_name Name of time variable
##' @param status_name Name of status variable
##' @param cause_codes Labels 
##' @param vars Additional variables
##' @return Fitted object with class jossl
##' ##' @examples
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
##' @export 
##' @author Anders Munch 
jossl <- function(learners,
                  data,
                  time,
                  integrate = TRUE,
                  split.method = "cv5",
                  collapse = TRUE,
                  B=1,
                  verbose = FALSE,
                  time_grid_length = 100,
                  time_name = "time",
                  status_name = "status",
                  cause_codes = c("1" = 1, "2" = 2, "c" = 0),
                  vars = NULL){
    ## Fit model on canonical data form
    fit_cause_model <- function(model, data, cause, x_form = NULL, ...){
        cause_names = c("cause1", "cause2", "censor")
        cause_names = cause_names[cause_names %in% names(data)]
        stopifnot(cause %in% cause_names, !("internal_event_var" %in% names(data)))
        wd = copy(data)[, -cause_names[cause_names != cause], with = FALSE]
        setnames(wd, old = cause, new = "internal_event_var")
        ## not so elegant -- should to this in a better way and export to separate functions
        out = NULL
        if(is.null(x_form))
            form = Surv(time, internal_event_var) ~ .
        else
            form = update(x_form, Surv(time, internal_event_var) ~ .)
        if(model == "cox"){
            wd[, internal_event_var := 1*(internal_event_var == 1)] ## To make this cause of interest
            out = survival::coxph(form, data = wd,x = TRUE,y = TRUE, ...)
        }
        if(model == "GLMnet"){
            wd[, internal_event_var := 1*(internal_event_var == 1)] ## To make this cause of interest
            out = riskRegression::GLMnet(form, data = wd, ...)
        }
        if(model == "rfsrc"){
            out = rfsrc(form, data = wd, ...)
        }
        if(model == "rfsrc.fast"){
            out = rfsrc.fast(form, data = wd, forest = TRUE, ...)
        }
        if(is.null(out))
            stop(paste("The model", model, "is not implemented in statelearner."))
        return(out)
    }
    ## util for rewrapping list of lists 
    rewrap <- function(ll, comb_fun = rbind){
        elements_in_list = length(ll[[1]])
        names_of_elements = names(ll[[1]])
        if(class(comb_fun) == "function")
            rep(c(comb_fun), times = elements_in_list)
        stopifnot(elements_in_list == length(comb_fun))
        out = lapply(1:elements_in_list, function(ii){
            do.call(comb_fun[[ii]], lapply(ll, function(xx){
                xx[[ii]]
            }))
        })
        names(out) = names_of_elements
        return(out)
    }
    ## Calculate F from CSCHFs
    abs_risk_from_cschf <- function(...){
        chfs = list(...)
        S = exp(-Reduce("+", chfs))
        Sminus = cbind(1,S[,-ncol(S)])
        abs_risk = lapply(chfs, function(xx) {
            xx[xx == Inf | is.na(xx)] <- max(xx[xx != Inf], na.rm = TRUE) ## Hack to avoid NaN's. Make better fix
            cs_haz = t(apply(cbind(0,xx), 1, diff))
            abs_risk_diff = cs_haz*Sminus
            return(riskRegression::rowCumSum(abs_risk_diff))
        })
        return(abs_risk)
    }
    comp_event_present = !is.null(learners$cause2)
    ## NB: Expecting that time variable is called time
    partition = riskRegression::getSplitMethod(split.method=split.method,
                                               B=B,
                                               N=NROW(data))
    ## Set data on canonical form
    wd = data.table::as.data.table(data.table::copy(data)[, c(time_name, status_name, vars), with = FALSE])
    data.table::setnames(wd, old = c(time_name, status_name), new = c("time", "status"))
    ## Make dummy variable for cause of interest = 1 and combine all other causes = 2
    wd[, cause1 := 1+1*(status != cause_codes[["1"]] )]
    wd[, censor := 1+1*(status != cause_codes[["c"]] )]
    if(wd[, length(unique(status))>2])
        wd[, cause2 := 1+1*(status != cause_codes[["2"]] )]
    stopifnot(length(time) == 1)
    time_grid = seq(0, time, length.out = time_grid_length)
    ## fix names
    if(is.null(names(learners$cause1)))
        names(learners$cause1) = paste0("cause1_", 1:length(learners$cause1))
    if(is.null(names(learners$censor)))
        names(learners$censor) = paste0("censor_", 1:length(learners$censor))
    if(comp_event_present){
        if(is.null(names(learners$cause2)))
            names(learners$cause2) = paste0("cause2_", 1:length(learners$cause2))
    }
    list_cv_fits = lapply(1:partition$B, function(bb){
        if(verbose) message("--- Running on split ", bb, " out of ", B, " ---")
        raw_out0 = lapply(1:partition$k, function(fold_k){
            if(verbose) message("------ Running on split ", fold_k, " out of ", partition$k)
            train = wd[partition$index(bb)!=fold_k]
            test = wd[partition$index(bb)==fold_k]
            ## Construct time grid to evaluate hazard function
            eval_times = sort(unique(c(0, train[, sort(unique(time))])))
            eval_times = eval_times[eval_times <= max(time_grid)]
            ## Construct list of cumulative hazard predictions
            ## This might give memory problems...
            time_fit = Sys.time()
            super_list_ch = lapply(names(learners), function(l_name){
                lapply(learners[[l_name]], function(mm){
                    fit = do.call(fit_cause_model,c(mm, list(data = train, cause = l_name)))
                    ch = predictCHF(fit, newdata = test, times = eval_times)
                    return(ch)
                })
            })
            ## ## Test
            ## fit_cause_model("GLMnet", train, "cause1")
            ## ## Test end            
            names(super_list_ch) = names(learners)
            list_ch_cause1 = super_list_ch[["cause1"]]
            list_ch_censor = super_list_ch[["censor"]]
            if(comp_event_present)
                list_ch_cause2 = super_list_ch[["cause2"]]
            time_fit = Sys.time() - time_fit
            ## Calculating Brier scores in holds out samples
            ## Construct matrices of events
            time_brier_calc = Sys.time()
            grid_mat = matrix(time_grid,ncol = length(time_grid),nrow = test[, .N],byrow = TRUE)
            event_mat = apply(grid_mat, 2, function(x) 1*(test[, time] <= x))
            cause1_event = apply(event_mat, 2, function(x) x*(test[, cause1 == 1]))
            censor_event = apply(event_mat, 2, function(x) x*(test[, censor == 1]))
            if(comp_event_present)
                cause2_event = apply(event_mat, 2, function(x) x*(test[, cause2 == 1]))
            if(comp_event_present)
                model_grid = expand.grid(cause1 = names((list_ch_cause1)),
                                         cause2 = names((list_ch_cause2)),
                                         censor = names((list_ch_censor)))
            else
                model_grid = expand.grid(cause1 = names((list_ch_cause1)),
                                         censor = names((list_ch_censor)))
            data.table::setDT(model_grid)
            brier_scores = do.call(rbind, lapply(1:nrow(model_grid), function(ii){
                base_dt = model_grid[ii]
                if(comp_event_present){
                    abs_risks = abs_risk_from_cschf(list_ch_cause1[[base_dt[, cause1]]],
                                                    list_ch_cause2[[base_dt[, cause2]]],
                                                    list_ch_censor[[base_dt[, censor]]])
                    event_counts = list(cause1_event, cause2_event, censor_event)
                }else{
                    abs_risks = abs_risk_from_cschf(list_ch_cause1[[base_dt[, cause1]]],
                                                    list_ch_censor[[base_dt[, censor]]])
                    event_counts = list(cause1_event, censor_event)
                }
                loss = sum(as.numeric(sapply(1:length(abs_risks), function(ll){
                    risk_pred = abs_risks[[ll]][, prodlim::sindex(jump.times = eval_times,eval.times = time_grid)]
                    obs_event = event_counts[[ll]]
                    ## Check this ok!
                    brier = (risk_pred-obs_event)^2
                    if (integrate)
                        loss0 = mean(rowSums(brier*diff(time_grid)[1]))
                    else
                        loss0 = mean(brier[, dim(brier)[2]])
                    return(loss0)
                })))
                out = copy(base_dt)
                out[, loss := loss]
                return(out)                
            }))
            time_brier_calc = Sys.time() - time_brier_calc
            return(list(brier_scores = brier_scores,
                        time_fit =  time_fit,
                        time_brier_calc = time_brier_calc))
        })
        out0 = rewrap(raw_out0, c(rbind, sum, sum))
        if(comp_event_present)
            out0$brier_scores = out0$brier_scores[, .(loss=mean(loss)), .(cause1, cause2, censor)]
        else
            out0$brier_scores = out0$brier_scores[, .(loss=mean(loss)), .(cause1, censor)]
        out0$brier_scores[, b:=bb]
        return(out0)
    })
    rw_list_cv_fits = rewrap(list_cv_fits, c(rbind, sum, sum))
    cv_fits = rw_list_cv_fits$brier_scores
    if(verbose){
        all_time_fit = rw_list_cv_fits$time_fit
        all_time_brier_calc = rw_list_cv_fits$time_brier_calc
        message("Spend ", round(as.numeric(all_time_fit)/60, digits = 2), " minutes fitting models\n",
                "and   ", round(as.numeric(all_time_brier_calc)/60, digits = 2), " minutes calculating Brier scores.")
    }        
    if (B>1){
        if(comp_event_present)
            ave_cv_fit = cv_fits[, .(loss=mean(loss), sd = sd(loss)), .(cause1, cause2, censor)]
        else
            ave_cv_fit = cv_fits[, .(loss=mean(loss), sd = sd(loss)), .(cause1, censor)]
    }else{
        ave_cv_fit = cv_fits
    }
    if(collapse == TRUE)
        cv_return <- ave_cv_fit
    else
        cv_return <- cv_fits
    setkey(cv_return,loss)
    ## Fit winner models
    cause1_fit = do.call(fit_cause_model, c(learners$cause1[[ave_cv_fit[1, cause1]]], list(data = wd, cause = "cause1")))
    censor_fit = do.call(fit_cause_model, c(learners$censor[[ave_cv_fit[1, censor]]], list(data = wd, cause = "censor")))
    if(comp_event_present){
        cause2_fit = do.call(fit_cause_model, c(learners$cause2[[ave_cv_fit[1, cause2]]], list(data = wd, cause = "cause2")))
        winners = list(cause1 = cause1_fit,cause2 = cause2_fit,censor = censor_fit)
    }else{
        winners = list(cause1 = cause1_fit,censor = censor_fit)
    }
    out = list(call = match.call(),
               cv_fit = cv_return,
               jossl = winners,
               times = sort(unique(wd$time)))
    class(out) = "jossl"
    return(out)
}


######################################################################
### jossl.R ends here
