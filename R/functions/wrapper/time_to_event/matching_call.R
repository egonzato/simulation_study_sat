matching_call=function(data,ratio,methods,landmark=730){
  # branche based on run
  duration_var=if ("duration" %in% names(data)) {rlang::sym("duration")} else {quote(NULL)}
  time_var=rlang::sym(if ("duration" %in% names(data)) "time" else "timetotal")  
  # branch based on methods
    if (methods=="exact_greedy"){
      matched=rlang::inject(matching(data=data,treatment=trt,ageofonset=ageonset,sex=sex,
                                     duration=!!duration_var,time=!!time_var,method="greedy",
                                     caliper=0.5,step=0.1,maxcaliper=0.5,ratio=ratio,
                                     replacement=F,history=F,warning=F))
    } else if (methods=="exact_optimal"){
      matched=rlang::inject(matching(data=data,treatment=trt,ageofonset=ageonset,sex=sex,
                                     duration=!!duration_var,time=!!time_var,method="optimal",
                                     caliper=NULL,step=NULL,maxcaliper=NULL,ratio=ratio,
                                     replacement=F,history=F,warning=F))
    } else if (methods=="psm_p_greedy_noc"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="glm",
                           method="greedy",
                           caliper=NULL,
                           duration=!!duration_var,
                           time=!!time_var))
    } else if (methods=="psm_p_greedy_caliper_01"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="glm",
                           method="greedy",
                           caliper=0.1,
                           duration=!!duration_var,
                           time=!!time_var))
    } else if (methods=="psm_p_greedy_caliper_005"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="glm",
                           method="greedy",
                           caliper=0.05,
                           duration=!!duration_var,
                           time=!!time_var))
    }else if (methods=="psm_m_greedy_noc"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="mahalanobis",
                           method="greedy",
                           caliper=NULL,
                           duration=!!duration_var,
                           time=!!time_var))
    }else if (methods=="psm_m_greedy_caliper_01"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="mahalanobis",
                           method="greedy",
                           caliper=0.1,
                           duration=!!duration_var,
                           time=!!time_var))
    }else if (methods=="psm_m_greedy_caliper_005"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="mahalanobis",
                           method="greedy",
                           caliper=0.05,
                           duration=!!duration_var,
                           time=!!time_var))
    }else if (methods=="psm_p_optimal"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="glm",
                           method="optimal",
                           caliper=NULL,
                           duration=!!duration_var,
                           time=!!time_var))
    } else if (methods=="psm_m_optimal"){
      matched=rlang::inject(psm_matching(data=data,
                           formula=trt~ageonset+sex,
                           ratio=ratio,
                           replacement=F,
                           distance="mahalanobis",
                           method="optimal",
                           caliper=NULL,
                           duration=!!duration_var,
                           time=!!time_var))
    } else if (methods=="duration_matching"){
      matched=rlang::inject(duration_matching(data=data,
                                              treatment=trt, 
                                              duration=!!duration_var,
                                              time=!!time_var,
                                              ratio=ratio))
    }
  # extract dataset
  matcheddf=matched$dataset %>% 
    mutate(event=ifelse(time_trial>landmark,0,event),
           time_trial=ifelse(time_trial>landmark,landmark,time_trial))
  # return
  matcheddf
}





