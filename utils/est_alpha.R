############################################
## estimate miscoverage rate 
## using estimated quantile of T
############################################
alpha_qt <- function(newdata,
                     # data_fit, 
                     data_calib,
                     covariates_fit, follow.up.time_fit,
                     pred_event_censor_obj,
                     xnames, alpha, len_x,
                     # mdl0, 
                     cens_rt){
  
  # v_list = v_pts_qt(pred_event_censor_obj,
  #                     # mdl, data_fit, data_calib,
  #                     xnames, alpha, cens_rt)
  # v_list = sort(unique(as.numeric(v_list)))
  v_list <- seq(0.01, 0.2, by = 0.002)
  
  ## Obtain the final confidence interval
  lower_bnd_l <- rep(0,len_x)
  lower_bnd_g <- rep(0,len_x)
  
  alpha_v_list <- sapply(v_list, est_alpha_qt, 
                    # mdl = mdl, 
                    data_calib = data_calib, 
                    pred_event_censor_obj=pred_event_censor_obj,
                    xnames = xnames, alpha = alpha, 
                    cens_rt = cens_rt, 
                    # mdl0 = mdl0, 
                    newdata = newdata)

  # monotonize alpha
  # alpha_v <- monot(alpha_v_list)
  # return 0 if alpha is above target level
  if(sum(alpha_v_list<=alpha)==0){
    v_hat_l = NULL
    v_hat_g = NULL
  }else{
    v_hat_l <- max(v_list[alpha_v_list <= alpha])
    # v_hat_l <- min(v_list[alpha_v_list <= alpha])
  }

  

  time.var <- "censored_T"
  event.var <- "event"

  event.surv.data <- cbind(covariates_fit, follow.up.time_fit)
  event.surv.test.data <- cbind(newdata, follow.up.time_fit)

  fmla_event <- as.formula(paste("Surv(",time.var,",",event.var,")",
                         paste(as.character("~ ."),collapse=""),
                         collapse=""))
  event.control=fit_surv_option(
      option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names),
                  cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names)))
  
  fit_surv_arg<-c(
      list(method="survSuperLearner",formula=fmla_event,data=event.surv.data,newdata=event.surv.test.data,time.var=time.var,event.var=event.var),
      event.control
  )
  pred_event_obj<-do.call(fit_surv, fit_surv_arg)
  
  # for(i in 1:len_x){
  #   nxi <- data.frame(newdata[i,])
  #   colnames(nxi) <- xnames
  #   lower_bnd_l[i] <- lv_qt(pred_event_censor_obj$event$surv, pred_event_censor_obj$event$time, v_hat_l)  
  # }

  lower_bnd_l <- apply(pred_event_obj$surv, 1, lv_qt, times = pred_event_obj$time, quantile = v_hat_l)
  # browser()
  
  return(list(lower_bnd_l = lower_bnd_l, 
              lower_bnd_g = lower_bnd_g))
  
}

############################################
## compute estimated miscoverage rate
############################################
est_alpha_qt <- function(
  pred_event_censor_obj,
  data_calib, 
  xnames, v, alpha, cens_rt, 
  # mdl0, 
  newdata){
  ## Check the dimensionality of the input
  if(is.null(dim(newdata)[1])){
    len_x <- length(newdata)
    p <- 1
  }else{
    len_x <- dim(newdata)[1]
    p <- dim(newdata)[2]
  }
  calib_x <- data.frame(X = data_calib[,names(data_calib) %in% xnames])
  names(calib_x) = xnames
  
  lv_calib <- apply(pred_event_censor_obj$event$surv, 1, lv_qt, times = pred_event_censor_obj$event$time, quantile = v)

  # lv_calib <- lv_qt(pred_event_censor_obj$event$surv, pred_event_censor_obj$event$time, v)

  # todo4: get P(C ≥ qa(X) | X = x)
  pr_calib <- apply(pred_event_censor_obj$censor$surv, 1, cens_prob, times = pred_event_censor_obj$censor$time, lv_calib)
  # cens_prob(pred_event_censor_obj$censor$surv, pred_event_censor_obj$censor$time, lv_calib)
  # browser()
  # cens = cens_prob(mdl0,data_calib,NULL,
  #                  method = "gpr",
  #                  xnames = xnames,
  #                  c = lv_calib)
  
  # pr_calib = cens$pr_calib
  weight_calib <- 1/pr_calib
  
  w_new = 0
  
  ind1 = (data_calib$censored_T < lv_calib) & (data_calib$C >= lv_calib)
  ind2 = (data_calib$C >= lv_calib)
  
  sum_num <- sum(weight_calib[ind1]) + w_new
  sum_den <- sum(weight_calib[ind2]) + w_new
  if((length(ind2)==0)||is.na(sum_num/sum_den)){
    alphav = 1
  }else{
    alphav <- sum_num / sum_den
  }
  
  return(alphav)
} 

lv_qt <- function(surv_probs, times, quantile) {
    index <- which.min(surv_probs >= 1-quantile)
    return(times[index])
}

cens_prob <- function(surv_probs, times, lv_calib) {
    index <- which.max(times >= lv_calib)
    # browser()
    return(surv_probs[index])
}

lv_qct <- function(pred_event_censor_obj, v, cens_rt){
  if(length(v)==0){
    return(lv_calib = 0)
  }
  lv1_calib <- apply(pred_event_censor_obj$censor$surv, 1, lv_qt, times = pred_event_censor_obj$censor$time, quantile = cens_rt)
  # lv_qt(pred_event_censor_obj$censor$surv, pred_event_censor_obj$censor$time, cens_rt)
  lv2_calib <- apply(pred_event_censor_obj$event$surv, 1, lv_qt, times = pred_event_censor_obj$event$time, quantile = v)
  # lv_qt(pred_event_censor_obj$event$surv, pred_event_censor_obj$event$time, v)

  # todo5: predict v quantile for both T and C
  # lv1_calib <- predict(qc_mdl, calib_x, cens_rt)
  # lv2_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1-v)
  # lv1_calib <- (lv1_calib$predictions)[,1]
  lv_calib = pmin(lv1_calib, lv2_calib)
  return(lv_calib)
}



################################################################
## The function returns  predictive intervals resulting
## from the L_v defined based on integrtaed quantiles
################################################################
alpha_qct <- function(newdata,
                     # data_fit, 
                     data_calib,
                     covariates_fit, follow.up.time_fit,
                     pred_event_censor_obj,
                     xnames, alpha, len_x,
                     # mdl0, 
                     cens_rt){

  # v_list = v_pts_qct(mdl, qc_mdl,
  #                     data_fit, data_calib,
  #                     xnames, alpha, cens_rt)
  # v_list = sort(unique(as.numeric(v_list)))
  v_list <- seq(0.01, 0.2, by = 0.002)

  ## Obtain the final confidence interval
  lower_bnd_l <- rep(0,len_x)
  lower_bnd_g <- rep(0,len_x)
  
  # alpha_v_list <- sapply(v_list, est_alpha_qct, mdl = mdl, qc_mdl = qc_mdl,
  #                   data_calib = data_calib, xnames = xnames, alpha = alpha, 
  #                   cens_rt = cens_rt, mdl0 = mdl0, newdata = newdata)

  alpha_v_list <- sapply(v_list, est_alpha_qct, 
                    # mdl = mdl, 
                    data_calib = data_calib, 
                    pred_event_censor_obj=pred_event_censor_obj,
                    xnames = xnames, alpha = alpha, 
                    cens_rt = cens_rt, 
                    # mdl0 = mdl0, 
                    newdata = newdata)
  # monotonize alpha
  # alpha_v <- monot(alpha_v_list)
  # return 0 if alpha is above target level
  if(sum(alpha_v_list<=alpha)==0){
    v_hat_l = NULL
    v_hat_g = NULL
  }else{
    v_hat_l <- min(v_list[alpha_v_list <= alpha])
  }
  time.var <- "censored_T"
  event.var <- "event"

  event.surv.data <- cbind(covariates_fit, follow.up.time_fit)
  event.surv.test.data <- cbind(newdata, follow.up.time_fit)

  fmla_event <- as.formula(paste("Surv(",time.var,",",event.var,")",
                         paste(as.character("~ ."),collapse=""),
                         collapse=""))
  event.control=fit_surv_option(
      option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names),
                  cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names)))
  
  fit_surv_arg<-c(
      list(method="survSuperLearner",formula=fmla_event,data=event.surv.data,newdata=event.surv.test.data,time.var=time.var,event.var=event.var),
      event.control
  )
  pred_event_obj<-do.call(fit_surv, fit_surv_arg)
  
  # for(i in 1:len_x){
  #   nxi <- data.frame(newdata[i,])
  #   colnames(nxi) <- xnames
  #   lower_bnd_l[i] <- lv_qt(pred_event_censor_obj$event$surv, pred_event_censor_obj$event$time, v_hat_l)  
  # }
  # lower_bnd_l <- apply(pred_event_obj$surv, 1, lv_qct, times = pred_event_obj$time, quantile = v_hat_l)
  lower_bnd_l <- lv_qct(pred_event_censor_obj, v_hat_l, cens_rt)
  
  # for(i in 1:len_x){
  #   nxi <- data.frame(newdata[i,])
  #   colnames(nxi) <- xnames

  #   lower_bnd_l[i] <- lv_qct(mdl, qc_mdl, nxi, v_hat_l, alpha, cens_rt)
  # }

  return(list(lower_bnd_l = lower_bnd_l, 
              lower_bnd_g = lower_bnd_g))

}

################################################################
## Compute estimated alpha using integrated quantiles
################################################################
est_alpha_qct <- function(
  pred_event_censor_obj,
  data_calib, 
  xnames, v, alpha, cens_rt, 
  # mdl0, 
  newdata){
  # mdl, qc_mdl, data_calib, xnames, v, alpha, cens_rt, mdl0, newdata){
  ## Check the dimensionality of the input
  if(is.null(dim(newdata)[1])){
    len_x <- length(newdata)
    p <- 1
  }else{
    len_x <- dim(newdata)[1]
    p <- dim(newdata)[2]
  }
  calib_x <- data.frame(X = data_calib[,names(data_calib) %in% xnames])
  names(calib_x) = xnames

  # lv_calib = lv_qct(mdl, qc_mdl, calib_x, v, alpha, cens_rt)
  
  # # todo4: get P(C ≥ qa(X) | X = x)
  # cens = cens_prob(mdl0, data_calib,NULL,
  #                  method = "gpr",
  #                  xnames = xnames,
  #                  c = lv_calib)

  lv_calib <- lv_qct(pred_event_censor_obj, v, cens_rt)

  # lv_calib <- lv_qt(pred_event_censor_obj$event$surv, pred_event_censor_obj$event$time, v)

  # todo4: get P(C ≥ qa(X) | X = x)
  pr_calib <- apply(pred_event_censor_obj$censor$surv, 1, cens_prob, times = pred_event_censor_obj$censor$time, lv_calib)
  # cens_prob(pred_event_c
  
  # pr_calib = cens$pr_calib
  weight_calib <- 1/pr_calib
  
  w_new = 0
  
  ind1 = (data_calib$censored_T < lv_calib) & (data_calib$C >= lv_calib)
  ind2 = (data_calib$C >= lv_calib)
  
  sum_num <- sum(weight_calib[ind1]) + w_new
  sum_den <- sum(weight_calib[ind2]) + w_new
  if((length(ind2)==0)||is.na(sum_num/sum_den)){
    alphav = 1
  }else{
    alphav <- sum_num / sum_den
  }

  return(alphav)
} 