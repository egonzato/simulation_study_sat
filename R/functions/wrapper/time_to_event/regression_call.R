regression_call=function(data,methods,landmark=730){
  # branch based on methods
  if (methods=="raw"){
    dataset=data
  } else if (methods=="iptw, unstabilised"){
    dataset=iptw(data=data,type='unstabilised',estimand='ATT',
                 formula=trt~ageonset+sex)[[1]]
  }else if (methods=="iptw, stabilised"){
    dataset=iptw(data=data,type='stabilised',estimand='ATT',
                 formula=trt~ageonset+sex)[[1]]
  }else if (methods=="iptw, stabilised, truncated 0.025"){
    dataset=iptw(data=data,type='stabilised',estimand='ATT',truncate=0.025,
                 formula=trt~ageonset+sex)[[1]]
  }else if (methods=="iptw, stabilised, truncated 0.05"){
    dataset=iptw(data=data,type='stabilised',estimand='ATT',truncate=0.05,
                 formula=trt~ageonset+sex)[[1]]
  } else if (methods=="gcomp"){
    dataset=data
  }
  dataset=dataset %>% 
    mutate(time_trial=timetotal,
           event=ifelse(time_trial>landmark,0,event),
           time_trial=ifelse(time_trial>landmark,landmark,time_trial))
}