# matching function
psm_matching=function(data=NULL,formula=NULL,duration=NULL,time=NULL,distance='glm',method='greedy',tune=FALSE,caliper=NULL,ratio=1,replacement=FALSE,seed=181299,warning=FALSE,history=FALSE){
  #  warnings to avoid wrong input
  if(!(replacement %in% c(TRUE,FALSE))){stop("Error: replacement must be either TRUE or FALSE")}
  if(replacement==FALSE & is.null(ratio)){stop("Please specify ratio when replacement = FALSE")}
  if(replacement==TRUE & method!='greedy'){stop("Optimal matching is only available without replacement")}
  if(replacement==TRUE & is.null(caliper) & is.null(ratio)){stop("Caliper has to be specified when replacement == TRUE and ratio==NULL")}
  if(!is.null(ratio)&&!(is.numeric(ratio)&&length(ratio)==1&&ratio==as.integer(ratio)&&ratio>0)){stop("Ratio must be an integer or set to NULL")}
  if(!(method %in% c('greedy','optimal','optimal.ilp'))){stop("Method specified not available,please choose one between 'greedy', 'optimal' and 'optimal.ilp'")}
  if(!(distance %in% c('glm','mahalanobis','random forest'))){stop("Method specified not available,please choose one between 'glm','mahalanobis'and 'random forest'")}
  if(tune==TRUE & distance!='random forest'){warning("tune=TRUE is only used when distance='random forest', ignoring")}
  # create id column to identify patients from the initial dataset
  data=data %>% mutate(id=1:n())
  #rownames(data) = as.character(data$id)
  # get variables
  if(!is.null(enexpr(duration))){duration_quo=enquo(duration)}; 
  time_quo=enquo(time)
  treatment_var=all.vars(formula)[1]
  treatment_quo=sym(treatment_var)
  # calculate probability for caliper difference
  data$ps=predict(glm(formula, data=data, family='binomial'), type="response")
  data$ps=pmax(pmin(data$ps, 1-1e-6), 1e-6)
  data$logit_ps=log(data$ps / (1 - data$ps))
  # define propensity score model with the distance given
  if (distance=='glm'){
    if (is.null(caliper)){
      data$distance=data$ps
    } else {
      data$distance=data$logit_ps
    }
    } else if (distance=='random forest'){
    if (tune==TRUE){
      data[[treatment_var]]=as.factor(data[[treatment_var]])
      set.seed(seed)
      tune_result=tuneRF(x=data[,all.vars(formula)[-1]], 
                           y=data[[treatment_var]],
                           ntreeTry = 500,
                           stepFactor = 1.5,
                           improve = 0.01,
                           trace = FALSE,
                           plot = FALSE)
      best_mtry=tune_result[which.min(tune_result[, 2]), 1]
      rf_model=randomForest(formula,data=data,ntree=500,mtry=best_mtry)
      data$distance=rf_model$votes[,2]
      data[[treatment_var]]=as.integer(as.character(data[[treatment_var]]))
    } else {
      data[[treatment_var]]=as.factor(data[[treatment_var]])
      set.seed(seed)
      rf_model=randomForest(formula,data=data,ntree=500)
      data$distance=rf_model$votes[,2]
      data[[treatment_var]]=as.integer(as.character(data[[treatment_var]]))
    }
  }
  # initialize matrix for the matched dataset and vector for unmatched
  matched=data.frame(); unmatched=c()
  # get treatment out of the formula 
  # divide into treated and untreated
  treated=data %>% filter(!!treatment_quo==1) %>% mutate(idx=row_number())
  untreated=data %>% filter(!!treatment_quo==0) %>% mutate(idx=row_number())
  # pull duration and distance before entering
  if(!is.null(enexpr(duration))){treated_duration=treated %>% pull(!!duration_quo)}
  if(!is.null(enexpr(duration))){untreated_duration=untreated %>% pull(!!duration_quo)}
  # depending on the distance, store values
  if(distance!='mahalanobis'){
    treated_distance=treated$distance
    untreated_distance=untreated$distance
  } else {
    cov_data = model.matrix(formula, data=data)[, -1, drop=FALSE]  # drop intercept column
    cov_inv  = solve(cov(cov_data))
    treated_cov   = model.matrix(formula, data=treated)[, -1, drop=FALSE]
    untreated_cov = model.matrix(formula, data=untreated)[, -1, drop=FALSE]
    }  
  # clipping distance
  if(distance != 'mahalanobis'){
    data$distance = pmax(pmin(data$distance,1-1e-6),1e-6)
  }
  # calculate caliper
  sdcaliper = if(!is.null(caliper)) caliper * sd(data$logit_ps) else NULL
  # start with methods
  if(method=='greedy'){
  # build filter expression dynamically: duration and caliper always included, sex conditional on exact
  if (!is.null(enexpr(duration))){
    if(!is.null(caliper)){
      filter_expr = expr((!!duration_quo) > current_duration &
                           logit_ps <= current_logit_ps + sdcaliper &
                           logit_ps >= current_logit_ps - sdcaliper)
    } else {
      filter_expr = expr((!!duration_quo) > current_duration)
    }
  } else{
    if(!is.null(caliper)){
      filter_expr = expr(logit_ps <= current_logit_ps + sdcaliper &
                           logit_ps >= current_logit_ps - sdcaliper)
    } else {
      filter_expr=TRUE
    }
  }
    # the first condition set identifies which direction the procedure takes, whether with or without replacement
    if(replacement==FALSE){
      # first loop looks for controls sequentially, since all treated patients need at least a control and then the procedure starts again for the j-th
      for(j in 1:ratio){
        for(i in 1:nrow(treated)){
          # at the beginning of each search the caliper is reset to its initial value
          filtered=NULL; controls=NULL; treatedp=NULL
          # extract values for current probability measure
          if(!is.null(enexpr(duration))){current_duration=treated_duration[i]}; current_logit_ps=treated$logit_ps[i]
          if(distance != 'mahalanobis') current_distance = treated$distance[i]
          # while loop runs until: i) one control is found or ii) current caliper is higher than the maximum allowed
          # filter for the matching criteria
          filtered = untreated %>%
            filter(!!filter_expr) %>%
            { if(distance != 'mahalanobis')
              mutate(., diff = abs(distance - current_distance))
              else {
                d = mutate(., diff = mahalanobis(untreated_cov[.$idx, ],
                                                 treated_cov[i, ],
                                                 cov_inv, inverted=TRUE))
                d
              }
            } %>%
            arrange(diff) %>%
            select(-diff)
          # if set to TRUE, prints the number of available controls for each treated participant (given the criteria), at each sequential search
          if(history==TRUE){print(paste('Patient n.:',i,'Control n.:',j,'Available:',nrow(filtered),sep=' '))}
          if(nrow(filtered)==0){
            # if the number of rows is equal to zero (i.e. no control was found), then the id is added to the unmatched vector
            # this can happen also if the first control was found but not the second one, so all individuals with a number of controls
            # lower than the ratio specified will appear in this vector
            unmatched=c(unmatched,treated$id[i])
            next
          }
          # if controls are found, then the loop stops, otherwise, the caliper is updated and the search starts with the new caliper
          # if caliper becomes higher than the maximum allowed, then the loop will break for individual i
          # at the end of the loop, whether stopped because a number of controls was found
          # or because the caliper was greater than the maximum allowed, the number of rows is evaluated
          # since the filtered dataset is arranged in descending order for year, then the first control is the most recent
          controls=filtered[1,]
          # create cluster variable, which in this case is the id value of the treated participant
          controls$cluster=treated$id[i]; treatedp=treated[i,]; treatedp$cluster=treated$id[i]
          # create disease duration for both controls and treated
          if(!is.null(enexpr(duration))){controls$enter=current_duration; treatedp$enter=current_duration}
          # bind together, both treated and control if this is the first iteration, otherwise only the control
          if(j==1){matched=rbind(matched,treatedp,controls)} else{matched=rbind(matched,controls)}
          # get control id to delete from the pool of available controls
          idxs=controls$idx; untreated=untreated %>% filter(!(idx %in% idxs))
        }
      }
    } else {
      for(i in 1:nrow(treated)){
        # at the beginning of each search the caliper is reset to its initial value
        filtered=NULL; controls=NULL; treatedp=NULL
        # extract values for current probability measure
        if(!is.null(enexpr(duration))){current_duration=treated_duration[i]}; current_logit_ps=treated$logit_ps[i]
        if(distance != 'mahalanobis') current_distance = treated$distance[i]
        # while loop runs until: i) one control is found or ii) current caliper is higher than the maximum allowed
        # filter for the matching criteria
        filtered = untreated %>%
          filter(!!filter_expr) %>%
          { if(distance != 'mahalanobis')
            mutate(., diff = abs(distance - current_distance))
            else {
              d = mutate(., diff = mahalanobis(untreated_cov[.$idx, ],
                                               treated_cov[i, ],
                                               cov_inv, inverted=TRUE))
              d
            }
          } %>%
          arrange(diff) %>%
          select(-diff)
        # if set to TRUE, prints the number of available controls for each treated participant (given the criteria), at each sequential search
        if(history==TRUE){print(paste('Patient n.:',i,'Available:',nrow(filtered),sep=' '))}
        # given that the argument ratio can be an integer or set to NULL, in matching with replacement it is needed to specify
        # the update of caliper in following way:
        # i) 1 to ratio, update if the number of controls is lower than ratio, up to the maximum allowed
        # ii) 1 to n, update until at least a control is found
        ratio.rep=ifelse(is.null(ratio),1,ratio)
        # if the number of rows is equal to zero (i.e. no control was found), then the id is added to the unmatched vector
        # if the number of controls is lower than the ratio specified (when it is not null), then add to the unmatched vector
        if(nrow(filtered)==0|(!is.null(ratio)&nrow(filtered)<ratio.rep)){unmatched=c(unmatched,treated$id[i])}
        # however, since every treated patients is priceless in terms of information, only those who have zero controls will be discarded
        if(nrow(filtered)==0){next}
        # if ratio == NULL, then keep all controls, otherwise pick just the number specified in the function
        n_controls=ifelse(is.null(ratio),nrow(filtered),min(ratio.rep,nrow(filtered)))
        controls=filtered[1:n_controls,]
        # create cluster variable, which in this case is the id value of the treated participant
        controls$cluster=treated$id[i]; treatedp=treated[i,]; treatedp$cluster=treated$id[i]
        # create disease duration for both controls and treated
        if(!is.null(enexpr(duration))){controls$enter=current_duration; treatedp$enter=current_duration}
        # bind together, both treated and controls if this is the first iteration, otherwise only the controls
        matched=rbind(matched,treatedp,controls)
      }
    }
  } 
  else if (method=='optimal'){
    # checks for optimal method
    if(is.null(ratio)){stop("Ratio must be specified for optimal matching")}
    if(!is.null(caliper) & method=='optimal'){warning("caliper not used when method=='optimal'")}
    # build distance matrix depending on distance method
    if(distance != 'mahalanobis'){
      distance_matrix = outer(treated_distance, untreated_distance, FUN=function(i,j) abs(i-j))
    } else {
      distance_matrix = matrix(0, nrow=nrow(treated), ncol=nrow(untreated))
      for(i in 1:nrow(treated)){
        for(j in 1:nrow(untreated)){
          distance_matrix[i,j] = mahalanobis(untreated_cov[j,], treated_cov[i,], cov_inv, inverted=TRUE)
        }
      }
    }
    # given the presence of categorical variables (duration filtering and sex) censor controls based on individual characteristics
    for(i in 1:nrow(treated)){
      # censor if duration is present
      if(!is.null(enexpr(duration))){duration_mask=untreated_duration<=treated_duration[i]}
      # set mask 
      if(!is.null(enexpr(duration))){
        distance_matrix[i,duration_mask]=Inf
      } else {
        distance_matrix=distance_matrix
      }
    }
    # define ids
    rownames(distance_matrix) = as.character(treated$id)
    colnames(distance_matrix) = as.character(untreated$id)
    # bind treated and untreated for data argument
    full_data=as.data.frame(bind_rows(treated,untreated))
    rownames(full_data)=c(as.character(treated$id),as.character(untreated$id))
    # controls = ratio specifies how many controls per treated patient
    mean_controls = min(ratio, nrow(untreated) / nrow(treated))
    match_result = fullmatch(distance_matrix,
                             min.controls = 1,
                             max.controls = ratio,
                             mean.controls = mean_controls,
                             remove.unmatchables = TRUE,
                             data = full_data)
    # # convert match result to dataframe
    match_df=data.frame(id=as.integer(names(match_result)),subclass=as.character(match_result),stringsAsFactors=FALSE) %>% filter(!is.na(subclass))
    # join back to original data to recover full information
    matched_treated=treated %>% inner_join(match_df,by="id")
    matched_untreated=untreated %>% inner_join(match_df,by="id")
    # assign cluster as the id of the treated patient in each subclass
    cluster_map=matched_treated %>% select(subclass,cluster=id)
    # bind treated and controls together with cluster variable
    matched=bind_rows(matched_treated,matched_untreated) %>%
      left_join(cluster_map,by="subclass") %>%
      select(-subclass) %>%
      group_by(cluster)  %>%
      { if(!is.null(enexpr(duration)))
        mutate(.,enter=(!!duration_quo)[(!!treatment_quo)==1][1]) 
        else .
      }
    # recover unmatched treated patients
    unmatched=treated %>% filter(!(id %in% match_df$id)) %>% pull(id)
  }
  else if (method=='optimal.ilp'){
    if(replacement==TRUE){stop("Optimal matching is only available without replacement")}
    if(is.null(ratio)){stop("Ratio must be specified for optimal matching")}
    #if(is.null(caliper)){stop("Caliper must be specified for optimal matching")}
    # define distance matrix
    # fill in the distance matrix with the difference of age of onset, which will be used for pairing up
    # treated on the rows and untreated on the columns
    if(distance != 'mahalanobis'){
      distance_matrix = outer(treated_distance, untreated_distance, FUN=function(i,j) abs(i-j))
    } else {
      distance_matrix = matrix(0, nrow=nrow(treated), ncol=nrow(untreated))
      for(i in 1:nrow(treated)){
        for(j in 1:nrow(untreated)){
          distance_matrix[i,j] = mahalanobis(untreated_cov[j,], treated_cov[i,], cov_inv, inverted=TRUE)
        }
      }
    }
    # apply duration and caliper masks
    for(i in 1:nrow(treated)){
      if(!is.null(enexpr(duration))){duration_mask=untreated_duration<=treated_duration[i]}
      #caliper_mask = abs(untreated$logit_ps - treated$logit_ps[i]) > sdcaliper
      #distance_matrix[i, duration_mask | caliper_mask] = Inf
      if(!is.null(enexpr(duration))){
        distance_matrix[i,duration_mask]=Inf
      } else {
        distance_matrix=distance_matrix
      }   
    }
    # run two-step ILP optimiser
    match_result=ilp_matching(distance_matrix,ratio,treated_ids=treated$id)
    pairs=match_result$pairs
    # recover matched dataset from pairs
    if(nrow(pairs)>0){
      for(k in 1:nrow(pairs)){
        i=pairs$treated_idx[k]; j=pairs$control_idx[k]
        treatedp=treated[i,]; treatedp$cluster=pairs$cluster[k]; 
        controlp=untreated[j,]; controlp$cluster=pairs$cluster[k]; 
        if(!is.null(enexpr(duration))){current_duration=treated_duration[i]}
        if(!is.null(enexpr(duration))){controlp$enter=current_duration; treatedp$enter=current_duration}
        matched=rbind(matched,treatedp,controlp)
      }
      # remove duplicate treated rows keeping only first occurrence per cluster
      matched=matched %>% group_by(cluster) %>% filter(!(!!treatment_quo==1 & row_number()>1)) %>% ungroup()
    }
    # recover unmatched treated patients
    unmatched=treated$id[match_result$unmatched]
  } 
  # then calculate the weights as explained in the readme file
  dftrial=matched %>%
    { if(!is.null(enexpr(duration)))
      mutate(.,time_trial=ifelse((!!treatment_quo)==0,(!!duration_quo)-enter,(!!time_quo)))
      else
        mutate(.,time_trial=(!!time_quo))
    } %>% 
    group_by(cluster) %>%
    mutate(count=max(row_number()),nmcluster=count-1,wmcluster=1/nmcluster) %>%
    ungroup() %>%
    group_by(idx,(!!treatment_quo)) %>%
    mutate(weight=sum(wmcluster)) %>%
    ungroup() %>%
    mutate(weight=ifelse((!!treatment_quo)==0,weight,1)) %>%
    select(-count,-nmcluster,-wmcluster) %>%
    filter(!is.na(time_trial)) %>%
    arrange(cluster,desc((!!treatment_quo)))
  # delete weights if replacement == FALSE
  if(replacement==FALSE){dftrial=dftrial %>% select(-weight)}
  # keep unique values of unmatched patients
  unmatched=unique(unmatched)
  # set unmatched to null if empty
  discarded=if(length(unmatched)>0) unmatched else NULL
  # warning unmatched
  if(warning&&length(unmatched)>0){warning("Some treated patients were unmatched")}
  # define results of optimal search
  optimalres=if(method %in% c('optimal','optimal.ilp')) match_result else NULL
  # list the robjects to return
  return(list(original=data,dataset=dftrial,optimalres=optimalres,discarded=discarded,method=method,caliper=caliper,ratio=ratio,replacement=replacement))
}