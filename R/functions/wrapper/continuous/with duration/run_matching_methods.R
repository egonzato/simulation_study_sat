run_matching_methods=function(data=NULL,methods=NULL,ratio=NULL,id=NULL,landmark=2,interpolation=T){
  # define quantities 
   message=NA_character_
   # extract id
   id_enquo=enquo(id)
   # call functions 
  out=tryCatch({
    withCallingHandlers({
      ## match
      matched=matching_call(data=data,ratio=ratio,methods=methods,id=id,interpolation=interpolation)
      ## one row per patient
      baseline=matched %>% 
        group_by(!!id_enquo) %>% 
        filter(row_number()==1) %>% 
        ungroup()
      ## run model
      res=analyse_call(data=matched,landmark=landmark,methods=methods)
      ## results
      data.frame(method=methods,
                 ratio=ratio,
                 result=res$ate,
                 singular=res$singular,
                 smdagea=smd(baseline$aoo,baseline$trt)$estimate,
                 smdsexa=smd(baseline$sex,baseline$trt)$estimate,
                 esst=NA_real_,
                 essc=NA_real_,
                 ess=NA_real_,
                 ntrt_post=unname(table(baseline$trt)[2]),
                 nctrl_post=unname(table(baseline$trt)[1]),
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
               ratio=ratio,
               result=NA_real_,
               singular=NA_real_,
               smdagea=NA_real_,
               smdsexa=NA_real_,
               esst=NA_real_,
               essc=NA_real_,
               ess=NA_real_,
               ntrt_post=NA_real_,
               nctrl_post=NA_real_,
               converged=FALSE,
               message=ifelse(is.na(message),
                              conditionMessage(e),
                              paste(message,conditionMessage(e),sep=" | ")))
  })
  out
}