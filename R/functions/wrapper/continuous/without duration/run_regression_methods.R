run_regression_methods=function(data=NULL,methods=NULL,ratio=NULL,id=NULL,landmark=2){
  # define quantities 
  message=NA_character_
  # extract id
  id_enquo=enquo(id)
  # call functions 
  out=tryCatch({
    withCallingHandlers({
      ## match
      dataset=regression_call(data,id,time="timen",methods)
      ## one row per patient
      baseline=dataset %>% 
        group_by(!!id_enquo) %>% 
        filter(row_number()==1) %>% 
        ungroup()
      ## run model
      res=analyse_call(data=dataset,landmark=landmark,methods=methods)
      # define quantities for ess
      if (grepl("iptw",methods)) {
        # weighted means
        smdagea=smd(baseline$aoo,baseline$trt,w=baseline$w)$estimate
        smdsexa=smd(baseline$sex,baseline$trt,w=baseline$w)$estimate
        # effectice sample size
        esst=sum(baseline$w[baseline$trt==1])^2/sum(baseline$w[baseline$trt==1]^2)
        essc=sum(baseline$w[baseline$trt==0])^2/sum(baseline$w[baseline$trt==0]^2)
        ess=esst+essc
      } else {
        # weighted means
        smdagea=NA_real_
        smdsexa=NA_real_
        # effective sample size
        esst=NA_real_
        essc=NA_real_
        ess=NA_real_
      }
      ## results
      data.frame(method=methods,
                 ratio=NA_real_,
                 result=res$ate,
                 singular=res$singular,                 
                 smdagea=smdagea,
                 smdsexa=smdsexa,
                 esst=esst,
                 essc=essc,
                 ess=ess,
                 event_post=NA_real_,
                 ntrt_post=NA_real_,
                 nctrl_post=NA_real_,
                 converged=TRUE,
                 message=message)
    },warning=function(w){
      # store message
      message<<-conditionMessage(w)
      # invoke restart
      invokeRestart("muffleWarning")
    })
  },error=function(e){
    ## results
    data.frame(method=methods,
               ratio=NA_real_,
               result=NA_real_,
               singular=NA_real_,
               smdagea=NA_real_,
               smdsexa=NA_real_,
               esst=NA_real_,
               essc=NA_real_,
               ess=NA_real_,
               event_post=NA_real_,
               ntrt_post=NA_real_,
               nctrl_post=NA_real_,
               converged=FALSE,
               message=ifelse(is.na(message),
                              conditionMessage(e),
                              paste(message,conditionMessage(e),sep=" | ")))
  })
  out
}