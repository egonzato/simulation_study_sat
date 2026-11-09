run_regression_methods=function(data=NULL,methods=NULL,landmark=730){
  # define quantities 
  message=NA_character_
  smdagea=NA_real_
  smdsexa=NA_real_
  esst=NA_real_
  essc=NA_real_
  ess=NA_real_
  # call functions 
  out=tryCatch({
    withCallingHandlers({
      ## match
      dataset=regression_call(data,methods)
      ## run model
      rr=analyse_call(data=dataset,landmark=landmark,methods=methods)
      # define quantities for ess
      if (grepl("iptw",methods)) {
        # weighted means
        smdagea=smd(dataset$ageonset,dataset$trt,w=dataset$w)$estimate
        smdsexa=smd(dataset$sex,dataset$trt,w=dataset$w)$estimate
        # effectice sample size
        esst=sum(dataset$w[dataset$trt==1])^2/sum(dataset$w[dataset$trt==1]^2)
        essc=sum(dataset$w[dataset$trt==0])^2/sum(dataset$w[dataset$trt==0]^2)
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
                 result=rr[[1]],
                 smdagea=smdagea,
                 smdsexa=smdsexa,
                 esst=esst,
                 essc=essc,
                 ess=ess,
                 event_post=unname(table(dataset$event)[2]),
                 ntrt_event_post=sum(dataset$event[dataset$trt==1]==1,na.rm=TRUE),
                 nctrl_event_post=sum(dataset$event[dataset$trt==0]==1,na.rm=TRUE),
                 ntrt_post=unname(table(dataset$trt)[2]),
                 nctrl_post=unname(table(dataset$trt)[1]),
                 converged=TRUE,
                 class=rr[[2]],
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
               smdagea=smdagea,
               smdsexa=smdsexa,
               esst=esst,
               essc=essc,
               ess=ess,
               event_post=NA_real_,
               ntrt_event_post=NA_real_,
               nctrl_event_post=NA_real_,
               ntrt_post=NA_real_,
               nctrl_post=NA_real_,
               converged=FALSE,
               class=NA_character_,
               message=ifelse(is.na(message),
                              conditionMessage(e),
                              paste(message,conditionMessage(e),sep=" | ")))
  })
  out
}