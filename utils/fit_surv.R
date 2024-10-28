#' @title Options of machine learning methods' wrappers for fitting conditional survival curves
#' @name fit_surv_option
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to the wrapped machine learning function. Will be used in a command like `do.call(machine.learning, option)` where `machine.learning` is the machine learning function being called. `formula` and `data` should not be specified. For \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}, if `tune=TRUE`, then `mtry` and `nodesize` should not be specified either.
#' @param oob whether to use out-of-bag (OOB) fitted values from random forests (\code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}, \code{\link[party:cforest]{party::cforest}}) and \code{\link[grf:survival_forest]{grf::survival_forest}}) when sample splitting is not used (`nfold=1`). Ignored otherwise.
#' @param tune whether to tune `mtry` and `nodesize` for \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. Ignored for other methods.
#' @param tune.option a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} is used and `tune=TRUE`; ignored otherwise. `doBest` should not be specified.
#' @param lambda bandwidth parameter for uniform smoothing kernel in nearest neighbours estimation for method `"akritas"`. The default value of 0.5 is arbitrary and should be chosen by the user
#' @export
fit_surv_option<-function(nfold=1,option=list(),oob=TRUE,tune=TRUE,tune.option=list(),lambda=0.5){
    assert_that(is.count(nfold))
    assert_that(is.flag(oob))
    assert_that(is.flag(tune))
    assert_that(is.number(lambda),lambda>0)
    out<-list(nfold=nfold,option=option,oob=oob,tune=tune,tune.option=tune.option,lambda=lambda)
    class(out)<-"fit_surv_option"
    out
}

fit_surv<-function(method=c("survSuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","survival_forest","no_event"),...){
    method<-match.arg(method)
    if(method=="survSuperLearner"){
        fit_survSuperLearner(...)
    }
}

#' @title Wrapper of `survSuperLearner::survSuperLearner`
#' @name fit_survSuperLearner
#' @param formula formula containing all covariates to be used
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survSuperLearner:survSuperLearner]{survSuperLearner::survSuperLearner}}. We encourage using a named list. Will be passed to \code{\link[survSuperLearner:survSuperLearner]{survSuperLearner::survSuperLearner}} by running a command like `do.call(survSuperLearner, option)`. The user should not specify `time`, `event`, `X`, or `newX`. We encourage the user to specify `event.SL.library` and `cens.SL.library`.
#' @param ... ignored
#' @return a \code{\link{pred_event_censor}} class containing fitted survival curves for individuals in `data`
#' @export
fit_survSuperLearner<-function(formula,data,newdata,time.var,event.var,nfold=1,option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc")),...){
    # .requireNamespace("survSuperLearner")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("time","event","X","newX") %in% names(option))){
        stop("option specifies time, event, X, or newX")
    }
    
    # formula<-as.formula(paste(as.character(formula)[-2],collapse=" "))
    if(nfold==1){
        time<-data%>%pull(time.var)
        event<-data%>%pull(event.var)
        # browser()
        
        # X<-model.frame(formula,data=data%>%select(!c(.data[[time.var]],.data[[event.var]])))
        X<-data%>%select(-all_of(c(time.var, event.var, "C")))
        # newX<-model.frame(formula,data=newdata%>%select(!c(.data[[time.var]],.data[[event.var]])))
        newX<-newdata%>%select(-all_of(c(time.var, event.var, "C")))

        new.time<-newdata%>%pull(time.var)
        # new.times<-sort(unique(new.time))
        new.times<-seq(0,max(new.time),length.out=1000) #t grid
        
        arg<-c(list(time=time,event=event,X=X,newX=newX,new.times=new.times),option)
        model<-do.call(survSuperLearner::survSuperLearner,arg)
        
        event.pred<-model$event.SL.predict
        # browser()
        # row.names(event.pred)<-data%>%pull(id.var)
        # pred_surv(time=new.times,surv=event.pred)
        # return the event prob at new.times
        return(list(time=new.times,surv=event.pred))
    }else{
        all.times<-data%>%pull(.data[[time.var]])%>%unique%>%sort
        
        all.times<-seq(min(all.times),max(all.times),length.out=1000) #t grid
        
        folds<-create.folds(pull(data,.data[[id.var]]),pull(data,.data[[event.var]]),nfold)
        
        # pred_event_censor.list
        pred_event.list<-lapply(folds,function(fold){
            d<-data%>%filter(!(.data[[id.var]] %in% .env$fold))
            test.d<-data%>%filter(.data[[id.var]] %in% .env$fold)
            
            time<-d%>%pull(time.var)
            event<-d%>%pull(event.var)
            
            if(all(event==0)){
                event.pred<-matrix(1,nrow=length(fold),ncol=length(all.times))
            }else{
                X<-model.frame(formula,data=d%>%select(!c(.data[[id.var]],.data[[time.var]],.data[[event.var]])))
                newX<-model.frame(formula,data=test.d%>%select(!c(.data[[id.var]],.data[[time.var]],.data[[event.var]])))
                new.times<-all.times
                
                arg<-c(list(time=time,event=event,X=X,newX=newX,new.times=new.times),option)
                model<-do.call(survSuperLearner::survSuperLearner,arg)
                
                event.pred<-model$event.SL.predict
            }
            
            row.names(event.pred)<-fold
            pred_surv(all.times,event.pred)
        })
        
        event.pred<-lapply(pred_event.list,function(x){
            x$surv
        })%>%do.call(what=rbind)
        event.pred<-event.pred[order(rownames(event.pred)),,drop=FALSE]
        pred_surv(all.times,event.pred)
    }
}
