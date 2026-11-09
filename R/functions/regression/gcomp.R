gcomp=function(data=NULL,indices,treatment=NULL,time=NULL,estimand="ATE",estimate="surv_hr",formula=NULL,landmark=365){
  # for warnings and errors
  converged=TRUE
  result=NA_real_
  # create index for bootstrap function
  data=data[indices,]
  # detect model type from formula LHS
  lhs=deparse(formula[[2]])
  # define target population
  target=if (estimand == "ATE") {
    data
  } else if (estimand=="ATT") {
    data[data[[treatment]]==1,,drop=FALSE]
  } else if (estimand=="ATC"){
    data[data[[treatment]]==0,,drop=FALSE]
  }
  # assign target to dataset which will be used for predictions
  newdata0=target
  newdata1=target
  # assign treatment value
  newdata0[[treatment]]=0
  newdata1[[treatment]]=1
  # time to event branch
  if(grepl('Surv',lhs)){
    tryCatch({
      # fit outcome model
      model=coxph(formula,data=data,x=T,ties='breslow',model=T)
      # predict
      predtrt=survfit(model,newdata=newdata1); meantrt=rowMeans(as.matrix(predtrt$surv))
      predctrl=survfit(model,newdata=newdata0); meanctrl=rowMeans(as.matrix(predctrl$surv))
      # find closest time point
      idx=which.min(abs(predtrt$time-landmark))
      # get risk ratio
      rr=risk_ratio(data=data,
                    model = model,
                    estimand=estimand,
                    treatment = treatment,
                    landmark = landmark)
      # get log HR
      loghr=log(-log(meantrt[idx]))-log(-log(meanctrl[idx]))
      # define result
      if(estimate=="surv_rr"){
        result=rr
      } else if(estimate=="surv_hr"){
        result=as.numeric(coef(model)[treatment])
      }
    },error=function(e){
      converged <<- FALSE
      result <<- NA_real_
    })
    # longitudinal branch
  } else {
    # fit outcome model
    fit=lmer(formula,data=data)
    # assign landmark evaluation point
    marginal=emmeans(fit, ~ trt_overall | timen,
                     at = list(timen=2),
                     cov.reduce=list(yoo=mean))
    # see contrast
    result=summary(contrast(marginal, method="trt.vs.ctrl", ref=1))$estimate[1]
  }
  c(estimate=result,converged=as.numeric(converged))
}