# matching function
duration_matching_interpolation=function(data=NULL,identifier=NULL,treatment=NULL,y=NULL,duration=NULL,time=NULL,ratio=NULL,replacement=FALSE,warning=TRUE,history=FALSE){
  #  warnings to avoid wrong input
  if(!(replacement %in% c(TRUE,FALSE))){stop("Error: replacement must be either TRUE or FALSE")}
  if(replacement==FALSE & is.null(ratio)){stop("Please specify ratio when replacement = FALSE")}
  if(!is.null(ratio)&&!(is.numeric(ratio)&&length(ratio)==1&&ratio==as.integer(ratio)&&ratio>0)){stop("Ratio must be an integer or set to NULL")}
  # create id column to identify patients from the initial dataset
  id_quo=enquo(identifier)
  data=data %>% 
    group_by(!!id_quo) %>% 
    mutate(idm=cur_group_id())
  # baseline covs
  baseline=data %>% 
    group_by(idm) %>% 
    filter(row_number()==1) %>% 
    ungroup()
  # column names to look into the dataset
  y_quo=enquo(y); duration_quo=enquo(duration); treatment_quo=enquo(treatment); time_quo=enquo(time)
  # initialize matrix for the matched dataset and vector for unmatched
  matched=data.frame(); unmatched=c()
  # divide into treated and untreated
  treated=baseline %>% filter(!!treatment_quo==1) %>% mutate(idx=row_number())
  untreated=baseline %>% filter(!!treatment_quo==0) %>% mutate(idx=row_number())
  # pull before entering the loop to avoid inefficiency
  treated_duration=treated %>% pull(!!duration_quo)
  # pull for untreated, needed for distance matrix
  untreated_duration=untreated %>% pull(!!duration_quo)
  # start with methods
  filter_expr=expr((!!duration_quo)>current_duration)
  # look for matches
  if(replacement==FALSE){
    # first loop looks for controls sequentially, since all treated patients need at least a control and then the procedure starts again for the j-th
    for(j in 1:ratio){
      for(i in 1:nrow(treated)){
        # at the beginning of each search the caliper is reset to its initial value
        filtered=NULL; controls=NULL; treatedp=NULL
        # extract values for current individual, for the matching criteria
        current_duration=treated_duration[i]; 
        # while loop runs until: i) one control is found or ii) current caliper is higher than the maximum allowed
        # filter for the matching criteria
        filtered=untreated %>% 
          filter(!!filter_expr)
        # if set to TRUE, prints the number of available controls for each treated participant (given the criteria), at each sequential search
        if(history==TRUE){print(paste('Patient n.:',i,'Control n.:',j,'Available:',nrow(filtered),sep=' '))}
        if(nrow(filtered)==0){
          # if the number of rows is equal to zero (i.e. no control was found), then the id is added to the unmatched vector
          # this can happen also if the first control was found but not the second one, so all individuals with a number of controls
          # lower than the ratio specified will appear in this vector
          unmatched=c(unmatched,treated$idm[i])
          next
        }
        # if controls are found, then the loop stops, otherwise, the caliper is updated and the search starts with the new caliper
        # if caliper becomes higher than the maximum allowed, then the loop will break for individual i
        # at the end of the loop, whether stopped because a number of controls was found
        # or because the caliper was greater than the maximum allowed, the number of rows is evaluated
        # since the filtered dataset is arranged in descending order for year, then the first control is the most recent
        controls=filtered[sample.int(nrow(filtered), 1),]
        # create cluster variable, which in this case is the id value of the treated participant
        controls$cluster=treated$idm[i]; treatedp=treated[i,]; treatedp$cluster=treated$idm[i]
        # create disease duration for both controls and treated
        controls$enter=current_duration; treatedp$enter=current_duration
        # bind together, both treated and control if this is the first iteration, otherwise only the control
        if(j==1){matched=rbind(matched,treatedp,controls)} else{matched=rbind(matched,controls)}
        # get control id to delete from the pool of available controls
        idxs=controls$idx; untreated=untreated %>% filter(!(idx %in% idxs))
      }
    }
  } else {
    for(i in 1:nrow(treated)){
      # at the beginning of each search the caliper is reset to its initial value
      current_caliper=caliper; filtered=NULL; controls=NULL; treatedp=NULL
      # while loop runs until: i) one control is found or ii) current caliper is higher than the maximum allowed
      # extract values for current individual
      current_duration=treated_duration[i]
      # filter for the matching criteria
      filtered=untreated %>% 
        filter(!!filter_expr)
      # if set to TRUE, prints the number of available controls for each treated participant (given the criteria), at each sequential search
      if(history==TRUE){print(paste('Patient n.:',i,'Available:',nrow(filtered),sep=' '))}
      # given that the argument ratio can be an integer or set to NULL, in matching with replacement it is needed to specify
      # the update of caliper in following way:
      # i) 1 to ratio, update if the number of controls is lower than ratio, up to the maximum allowed
      # ii) 1 to n, update until at least a control is found
      ratio.rep=ifelse(is.null(ratio),1,ratio)
      # if controls are found, then the loop stops, otherwise, the caliper is updated and the search starts with the new caliper
      if(nrow(filtered)>=ratio.rep){break} else{current_caliper=current_caliper+step}
      # if the number of rows is equal to zero (i.e. no control was found), then the id is added to the unmatched vector
      # if the number of controls is lower than the ratio specified (when it is not null), then add to the unmatched vector
      if(nrow(filtered)==0|(!is.null(ratio)&nrow(filtered)<ratio.rep)){unmatched=c(unmatched,treated$idm[i])}
      # however, since every treated patients is priceless in terms of information, only those who have zero controls will be discarded
      if(nrow(filtered)==0){next}
      # if ratio == NULL, then keep all controls, otherwise pick just the number specified in the function
      n_controls=ifelse(is.null(ratio),nrow(filtered),min(ratio.rep,nrow(filtered)))
      controls=filtered[1:n_controls,]
      # create cluster variable, which in this case is the id value of the treated participant
      controls$cluster=treated$idm[i]; treatedp=treated[i,]; treatedp$cluster=treated$idm[i]
      # create disease duration for both controls and treated
      controls$enter=current_duration; treatedp$enter=current_duration
      # bind together, both treated and controls if this is the first iteration, otherwise only the controls
      matched=rbind(matched,treatedp,controls)
    }
  }
  # from matched, keep only id, cluster and enter
  matched=matched %>% 
    group_by(cluster) %>% 
    mutate(count=max(row_number()),
           nmcluster=count-1,
           wmcluster=1/nmcluster) %>% 
    ungroup() %>% 
    select(-count,-nmcluster) %>% 
    select(idm,cluster,enter,wmcluster) 
  # merge initial dataset to merged
  matched_merged=left_join(matched,data,by='idm')
  # create timeline and delete previous duration for trial participants
  dftrial=matched_merged %>% 
    mutate(eweight=ifelse((!!treatment_quo)==0,wmcluster,1),
           time_aligned=(!!time_quo)-enter)
  ## interpolated
  ### before time zero
  before_zero=dftrial %>% 
      filter(!((!!treatment_quo)==1 & time_aligned<0)) %>%
      filter(time_aligned<0) %>%
      group_by(idm) %>%
      filter(time_aligned==max(time_aligned)) %>%
      rename(time_before=time_aligned,y_before=(!!y_quo))  %>% 
      ungroup()
  ### after zero
  after_zero=dftrial %>% 
    filter(!((!!treatment_quo)==1 & time_aligned<0)) %>%
    filter(time_aligned>=0) %>%
    group_by(idm) %>%
    filter(time_aligned==min(time_aligned)) %>%
    rename(time_after=time_aligned, 
           y_after=(!!y_quo)) %>%  
    select(idm,time_after,y_after) %>% 
    ungroup()
  ### merge
  interpolation_points=before_zero %>%
      inner_join(after_zero, by = "idm")
  ### interpolate
  interpolation_points2=interpolation_points %>%
      mutate(y_zero=y_before+(0-time_before)/(time_after-time_before)*(y_after - y_before),
             time_zero=0) %>%
      rename(time_aligned=time_zero,y=y_zero)
  ### get after zero matched dataset
  matched_after=dftrial %>% 
      filter(time_aligned>=0) %>% 
      mutate(y=(!!y_quo)) %>% 
      select(-(!!y_quo))
  ### row bind
  dftrial_aligned_interpolated=bind_rows(matched_after,interpolation_points2) %>% 
      arrange(idm,time_aligned)
  ## aligned
    dftrial_filtered=dftrial %>% 
      filter(time_aligned>=0) %>% 
      arrange(idm,time_aligned)
  # delete weights if replacement == FALSE
  if (replacement==FALSE){dftrial=dftrial %>% select(-eweight)}
  # keep unique values of unmatched patients 
  unmatched=unique(unmatched)
  # set unmatched to null if empty
  discarded=if(length(unmatched)>0) unmatched else NULL
  # warning unmatched
  if (warning && length(unmatched)>0) {warning("Some treated patients were unmatched")}
  # list the robjects to return
  return(list(aligned=dftrial,filtered=dftrial_filtered,interpolated=dftrial_aligned_interpolated,discarded=discarded,ratio=ratio,replacement=replacement))
}

