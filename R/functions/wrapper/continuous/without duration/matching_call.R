matching_call=function(data=NULL,id=NULL,ratio,methods){
  # branche based on run
  duration_var=if ("duration" %in% names(data)) {rlang::sym("duration")} else {quote(NULL)}
  time_var=rlang::sym(if ("duration" %in% names(data)) "timen" else "timen")  
  id_enquo=enquo(id)
  id_name=rlang::as_name(id_enquo)
  # create baseline
  baseline=data %>% 
    group_by(!!id_enquo) %>% 
    filter(row_number()==1) %>% 
    ungroup()
  # extract matches
  extract_matches = function(matched) {
    matched$dataset %>%
      ungroup() %>%
      transmute(!!id_name := baseline[[id_name]][id],
                cluster    = baseline[[id_name]][cluster])
  }
  # branch based on methods
    if (methods=="exact_greedy"){
      matched=rlang::inject(matching(data=baseline,treatment=trt_overall,ageofonset=aoo,sex=sex,
                                     duration=!!duration_var,time=!!time_var,method='greedy',
                                     caliper=0.5,step=0.1,maxcaliper=0.5,ratio=ratio,
                                     replacement=F,history=F,warning=F))
        # merge with full data
      matches = extract_matches(matched)
    } else if (methods=="exact_optimal"){
      matched=rlang::inject(matching(data=baseline,treatment=trt_overall,ageofonset=aoo,sex=sex,
                                     duration=!!duration_var,time=!!time_var,method='optimal',
                                     caliper=NULL,step=NULL,maxcaliper=NULL,ratio=ratio,
                                     replacement=F,history=F,warning=F))
      # merge with full data
      matches = extract_matches(matched)
    } else if (methods=="psm_p_greedy_noc"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='glm',
                                         method='greedy',
                                         caliper=NULL,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    } else if (methods=="psm_p_greedy_caliper_01"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='glm',
                                         method='greedy',
                                         caliper=0.1,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    } else if (methods=="psm_p_greedy_caliper_005"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='glm',
                                         method='greedy',
                                         caliper=0.05,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    }else if (methods=="psm_m_greedy_noc"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='greedy',
                                         caliper=NULL,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    }else if (methods=="psm_m_greedy_caliper_01"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='greedy',
                                         caliper=0.1,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    }else if (methods=="psm_m_greedy_caliper_005"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='greedy',
                                         caliper=0.05,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    }else if (methods=="psm_p_optimal"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='glm',
                                         method='optimal',
                                         caliper=NULL,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    }else if (methods=="psm_m_optimal"){
      matched=rlang::inject(psm_matching(data=baseline,
                                         formula=trt~aoo+sex,
                                         ratio=ratio,
                                         replacement=F,
                                         distance='mahalanobis',
                                         method='optimal',
                                         caliper=NULL,
                                         duration=!!duration_var,
                                         time=!!time_var))
      # merge with full data
      matches = extract_matches(matched)
    }
  # merge
  trialdf_matched = merge(data, matches, by = id_name)
  # add baseline value
  trialdf_matched %>%
    group_by(!!id_enquo) %>%
    arrange(!!enquo(time_var)) %>% 
    mutate(hui_baseline = first(hui_trt)) %>%
    filter(row_number()!=1) %>% 
    ungroup()
}





