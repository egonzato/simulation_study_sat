analyse_call=function(data,landmark=730,methods){
  # branch based on methods
  ## iptw
  if (grepl("iptw",methods)){
    weights_quo=enquo(weights)
    # run cox model
    model=tryCatch(
      withCallingHandlers(
        coxph(Surv(time_trial,event)~trt+calendar,x=T,data=data,weights=w),
        warning = function(w) {
          if (grepl("coefficient may be infinite|did not converge",
                    conditionMessage(w))) {
            stop("convergence_warning")
          }
        }
      ),
      error = function(e) {
        if (conditionMessage(e) == "convergence_warning") {
          coxph(Surv(time_trial,event)~ridge(trt,theta=1)+calendar,x=T,data=data,weights=w)
        } else {
          stop(e)
        }
      }
    )    # calculate risk ratio
    rr=risk_ratio(data,model=model,treatment="trt",estimand='ATT',landmark=landmark)
  } else if (methods=="raw"){
    model=tryCatch(
      withCallingHandlers(
        coxph(Surv(time_trial,event)~trt+calendar,x=T,data=data),
        warning = function(w) {
          if (grepl("coefficient may be infinite|did not converge",
                    conditionMessage(w))) {
            stop("convergence_warning")
          }
        }
      ),
      error = function(e) {
        if (conditionMessage(e) == "convergence_warning") {
          coxph(Surv(time_trial,event)~ridge(trt,theta=1)+calendar,x=T,data=data)
        } else {
          stop(e)
        }
      }
    )
    # calculate risk ratio
    rr=risk_ratio(data,model=model,treatment="trt",estimand='ATT',landmark=landmark)
  } else if (methods=="gcomp"){
    ## rr
    rr=tryCatch(
      withCallingHandlers(
        gcomp(data=data,
              estimand=c("ATT"),
              treatment="trt",
              estimate="surv_rr",
              formula=Surv(time_trial,event)~trt+calendar+ageonset+sex,
              landmark=landmark),
        warning = function(w) {
          if (grepl("coefficient may be infinite|did not converge",
                    conditionMessage(w))) {
            stop("convergence_warning")
          }
        }
      ),
      error = function(e) {
        if (conditionMessage(e) == "convergence_warning") {
          gcomp(data=data,
                estimand=c("ATT"),
                treatment="trt",
                estimate="surv_rr",
                formula=Surv(time_trial,event)~ridge(trt,theta=1)+calendar+ageonset+sex,
                landmark=landmark)
        } else {
          stop(e)
        }
      }
    )
  } else {
    # run model
    model=tryCatch(
      withCallingHandlers(
        coxph(Surv(time_trial,event)~trt+calendar,x=T,data=data,cluster=cluster,robust=T),
        warning = function(w) {
          if (grepl("coefficient may be infinite|did not converge",
                    conditionMessage(w))) {
            stop("convergence_warning")
          }
        }
      ),
      error = function(e) {
        if (conditionMessage(e) == "convergence_warning") {
          coxph(Surv(time_trial,event)~ridge(trt,theta=1)+calendar,data=data,cluster=cluster,robust=T)
        } else {
          stop(e)
        }
      }
    )
    # calculate risk ratio
    rr=risk_ratio(data,model=model,treatment="trt",estimand='ATT',landmark=landmark)
  }
  # return
  rr
}