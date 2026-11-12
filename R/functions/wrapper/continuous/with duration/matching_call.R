matching_call=function(data=NULL,id=NULL,ratio,methods,interpolation=T){
  # branche based on run
  duration_var=if ("duration" %in% names(data)) {rlang::sym("duration")} else {quote(NULL)}
  time_var=rlang::sym(if ("duration" %in% names(data)) "timen" else "timen")  
  id_enquo=enquo(id)
  # branch based on methods
    if (methods=="duration_matching"){
      matched=rlang::inject(duration_matching_interpolation(data=data,treatment=trt_overall,identifier=id,y=hui_trt,
                                                            duration=!!duration_var,time=!!time_var,ratio=ratio,
                                                            replacement=F,history=F,warning=F))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    }  else if (methods=="exact_greedy"){
      matched=rlang::inject(matching_interpolation(data=data,treatment=trt_overall,ageofonset=aoo,sex=sex,
                                                   identifier=id,y=hui_trt,
                                                   duration=!!duration_var,time=!!time_var,method='greedy',
                                                   caliper=0.5,step=0.1,maxcaliper=0.5,ratio=ratio,
                                                   replacement=F,history=F,warning=F))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    }  else if (methods=="exact_optimal"){
      matched=rlang::inject(matching_interpolation(data=data,treatment=trt_overall,ageofonset=aoo,sex=sex,
                                                   identifier=id,y=hui_trt,
                                                   duration=!!duration_var,time=!!time_var,method='optimal',
                                                   caliper=NULL,step=NULL,maxcaliper=NULL,ratio=ratio,
                                                   replacement=F,history=F,warning=F))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    } else if (methods=="psm_p_greedy_noc"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='glm',
                                         method='greedy',
                                         caliper=NULL,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    } else if (methods=="psm_p_greedy_caliper_01"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='glm',
                                         method='greedy',
                                         caliper=0.1,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    } else if (methods=="psm_p_greedy_caliper_005"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='glm',
                                         method='greedy',
                                         caliper=0.05,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    }else if (methods=="psm_m_greedy_noc"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='greedy',
                                         caliper=NULL,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    }else if (methods=="psm_m_greedy_caliper_01"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='greedy',
                                         caliper=0.1,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    }else if (methods=="psm_m_greedy_caliper_005"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='greedy',
                                         caliper=0.05,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    } else if (methods=="psm_p_optimal"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='glm',
                                         method='optimal',
                                         caliper=NULL,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    } else if (methods=="psm_m_optimal"){
      matched=rlang::inject(psm_matching_interpolation(data=data,
                                         formula=trt_overall~aoo+sex,
                                         ratio=ratio,
                                         identifier=id,
                                         outcome=hui_trt,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='optimal',
                                         caliper=NULL,
                                         warning=F,
                                         duration=!!duration_var,
                                         time=!!time_var))
      if (interpolation==T){
        # extract dataset
        trialdf_matched=matched$interpolated %>% 
          filter(time_aligned<=2) %>% 
          mutate(hui_trt=y) %>% 
          select(-y)
      } else {
        # extract dataset
        trialdf_matched=matched$aligned %>% 
          filter(time_aligned>=0)
      }
    }
  # return
  trialdf_matched %>% 
    group_by(!!id_enquo) %>% 
    arrange(!!enquo(time_var)) %>% 
    mutate(hui_baseline = first(hui_trt)) %>%
    filter(row_number()!=1) %>% 
    ungroup()
}





