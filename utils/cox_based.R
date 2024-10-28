############################################
## Lower prediction bound based on Cox model
## with adaptive cutoffs
############################################

cox_based <- function(x,alpha,
                      data_fit,
                      data_calib,
                      dist, 
                      # mdl0,
                      event.control,
                     censor.control){
  
  ## Check the dimensionality of the input
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  time.var <- "censored_T"
  event.var <- "event"

  n <- nrow(data_calib)
  ## Fit the survival model
  xnames <- paste0("X",1:p)
  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  # browser()
  covariates <- data_fit[, c(xnames), drop=FALSE]
  follow.up.time <- data_fit[, c("censored_T", "event", "C")]

  cal.covariates <- data_calib[, c(xnames), drop=FALSE]
  cal.follow.up.time <- data_calib[, c("censored_T", "event", "C")]

  event.surv.data <- cbind(covariates, follow.up.time)
  # cal.event.follow.up.time<-follow.up.time
  event.surv.cal.data <- cbind(cal.covariates, cal.follow.up.time)
  
  fmla_event <- as.formula(paste("Surv(",time.var,",",event.var,")",
                         paste(as.character(~ .),collapse=""),
                         collapse=""))

  fit_surv_arg<-c(
      list(method="survSuperLearner",formula=fmla_event,data=event.surv.data,newdata=event.surv.cal.data,time.var=time.var,event.var=event.var),
      event.control
  )

  pred_event_obj<-do.call(fit_surv, fit_surv_arg)

  # todo1: seems like no need to left.shift
  censor.follow.up.time<-follow.up.time%>%
                mutate("{event.var}":=1-.data[[event.var]],
                       "{time.var}":=left.shift.censoring(.data[["C"]],.data[[event.var]]))
  censor.surv.data <- cbind(covariates, censor.follow.up.time)
  cal.censor.follow.up.time<-cal.follow.up.time%>%
              mutate("{event.var}":=1-.data[[event.var]],
                     "{time.var}":=left.shift.censoring(.data[["C"]],.data[[event.var]]))
  censor.surv.cal.data <- cbind(cal.covariates, cal.censor.follow.up.time)

  fmla_censor <- as.formula(paste("Surv(",time.var,",",event.var,")",
                         paste(as.character(~ .),collapse=""),
                         collapse=""))
  
  fit_surv_arg<-c(
      list(method="survSuperLearner",formula=fmla_censor,data=censor.surv.data,newdata=censor.surv.cal.data,time.var=time.var,event.var=event.var),
      censor.control
  )
  pred_censor_obj<-do.call(fit_surv, fit_surv_arg)    

  pred_event_censor_obj<-list(event=pred_event_obj,censor=pred_censor_obj)

  # # todo2: fit T|X
  # fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  # mdl <- survreg(fmla, data = data_fit, dist= "weibull")

  ## The fitted quantile for the calibration data
  # xdf <- data.frame(data_calib[,names(data_calib) %in% xnames])
  # colnames(xdf) = colnames(newdata)
  
  # cutoff for quantile of C
  # to bound weights
  cens_rt = 1-1/log(n)
  
  start_time = proc.time()[3]
  qt_res <- alpha_qt(newdata, data_calib,
                    covariates, follow.up.time,
                     # data_fit, 
                     pred_event_censor_obj,
                     xnames, alpha, len_x,
                     cens_rt = cens_rt)
                     # mdl0 = mdl0)
  lower_bnd_qtg <- qt_res$lower_bnd_g
  lower_bnd_qtl <- qt_res$lower_bnd_l
  end_time <- proc.time()[3]
  time_qt <- end_time - start_time
  
  start_time = proc.time()[3]
  ## Fit the model for C with quantile_forest (now only supports 1d)
  # fit_X <- data.frame(X = data_fit[,names(data_fit) %in% xnames])


  # qc_mdl <- quantile_forest(fit_X, as.vector(data_fit$C))

  # browser()
  
  # qct_res <- alpha_qct(newdata, data_calib,
  #                   covariates, follow.up.time,
  #                    # data_fit, data_calib,
  #                    pred_event_censor_obj,
  #                    xnames, alpha, len_x,
  #                    cens_rt = cens_rt)

  # lower_bnd_qctg <- qct_res$lower_bnd_g
  # lower_bnd_qctl <- qct_res$lower_bnd_l
  # end_time <- proc.time()[3]
  # time_qct <- end_time - start_time
  

  return(list(lower_bnd_qtg =  lower_bnd_qtg,
              lower_bnd_qtl = lower_bnd_qtl,
              # lower_bnd_qctg = lower_bnd_qctg,
              # lower_bnd_qctl = lower_bnd_qctl,
              time_qt = time_qt
              # time_qct = time_qct
              ))
}



