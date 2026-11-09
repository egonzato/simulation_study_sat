risk_ratio = function(model, data, treatment, estimand="ATT", landmark=730){
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
  if (class(model)[1]=="coxph"){
    # calculate quantities to predict
    bh=basehaz(model,centered=TRUE)
    H0=approx(bh$time,bh$hazard,xout=landmark,method="constant",rule=2)$y
    S0=exp(-H0)
  } else {
    # calculate quantities to predict
    y=model$y
    time=y[, "time"]
    status=y[, "status"]
    lp=model$linear.predictors
    # breslow baseline cumulative hazard
    event_times=sort(unique(time[status == 1]))
    # calculate dh per time
    dH=vapply(event_times, function(t) {
      sum(status == 1 & time == t) /
        sum(exp(lp[time >= t]))
    }, numeric(1))
    # calculate baseline hazard and surival given landmark
    H0=sum(dH[event_times <= landmark])  
    S0=exp(-H0)
  }
  # predict
  lp0=predict(model, newdata=newdata0, type="lp")
  lp1=predict(model, newdata=newdata1, type="lp")
  # calcualte both risks
  risk0=1-mean(S0^exp(lp0))
  risk1=1-mean(S0^exp(lp1))
  # calculate risk ratio
  rr=risk1 / risk0
  # return risk ratio
  return(list(rr=rr,model=class(model)[1]))
}