analyse_call=function(data,landmark=2,methods){
  # branch based on methods
  ## iptw
  if (grepl("iptw",methods)){
    # linear mixed model
    fit=suppressMessages(lmer(hui_trt~(1|id)+hui_baseline+trt*timen+yoo,data=data,weights=w))
  } else if (methods=="raw"){
    # linear mixed model
    fit=suppressMessages(lmer(hui_trt~(1|id)+hui_baseline+trt*timen+yoo,data=data))
  } else if (methods=="gcomp"){
    # linear mixed model
    fit=suppressMessages(lmer(hui_trt~(1|id)+hui_baseline+trt*timen+yoo+aoo+sex,data=data))
  } else {
    fit=suppressMessages(lmer(hui_trt~(1|cluster)+(1|cluster:id)+hui_baseline+trt*timen+yoo,data=data))
    }
  # marginalize
  marginal=emmeans(fit, ~ trt | timen,
                   at = list(timen=landmark),
                   cov.reduce=list(yoo=mean))
  # see contrast
  ate=summary(contrast(marginal, method="trt.vs.ctrl", ref=1))$estimate
  # return
  list(ate=ate,
       singular=isSingular(fit,tol=1e-4))
  }