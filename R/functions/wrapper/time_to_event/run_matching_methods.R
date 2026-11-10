run_matching_methods=function(data=NULL,methods=NULL,ratio=NULL,landmark=730){
  # define quantities 
   message=NA_character_
   # define quantities for ess
   # call functions 
  out=tryCatch({
    withCallingHandlers({
      ## match
      matched=matching_call(data=data,ratio=ratio,methods=methods,landmark=landmark)
      ## run model
      rr=analyse_call(matched,landmark=landmark,methods=methods)
      ## results
      data.frame(method=methods,
                 ratio=ratio,
                 result=rr[[1]],
                 smdagea = smd(matched$ageonset, matched$trt)$estimate,
                 smdsexa = smd(matched$sex, matched$trt)$estimate,
                 esst=NA_real_,
                 essc=NA_real_,
                 ess=NA_real_,
                 event_post=sum(matched$event==1,na.rm=TRUE),
                 ntrt_event_post=sum(matched$event[matched$trt==1]==1,na.rm=TRUE),
                 nctrl_event_post=sum(matched$event[matched$trt==0]==1,na.rm=TRUE),
                 ntrt_post=unname(table(matched$trt)[2]),
                 nctrl_post=unname(table(matched$trt)[1]),
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
               ratio=ratio,
               result=NA_real_,
               smdagea=NA_real_,
               smdsexa=NA_real_,
               esst=NA_real_,
               essc=NA_real_,
               ess=NA_real_,
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