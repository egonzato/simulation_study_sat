analyse_call=function(data,landmark=2,methods){
  # branch based on methods
  ## iptw
  if (grepl("iptw",methods)){
    # linear mixed model
    fit=suppressMessages(lmer(hui_linear_trt~(1|id)+hui_baseline+trt*timen+yoo,data=data,weights=w))
  } else if (methods=="raw"){
    # linear mixed model
    fit=suppressMessages(lmer(hui_linear_trt~(1|id)+hui_baseline+trt*timen+yoo,data=data))
  } else if (methods=="gcomp"){
    # linear mixed model
    fit=suppressMessages(lmer(hui_linear_trt~(1|id)+hui_baseline+trt*timen+yoo+aoo+sex,data=data))
  } else {
    # enter (duration at alignment) is needed next to yoo: among the treated yoo tracks the duration almost exactly
    fit=suppressMessages(lmer(hui_linear_trt~(1|cluster)+(1|cluster:id)+hui_baseline+trt*time_aligned+yoo+enter,data=data))
    }
  # marginalize
  marginal=emmeans(fit, ~ trt | time_aligned,
                   at = list(time_aligned=landmark),
                   cov.reduce=list(yoo=mean))
  # see contrast
  # treated minus untreated (pairs() would give untreated minus treated)
  ate=summary(contrast(marginal, method="trt.vs.ctrl", ref=1))$estimate[1]
  # return
  list(ate=ate,
       singular=isSingular(fit,tol=1e-4))
  }