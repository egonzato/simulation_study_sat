# inverse probability of weighting function
iptw=function(data=NULL,formula=NULL,type='stabilised',estimand='ATE',truncate=NULL,distance='glm',tune=NULL){
  # warnings
  if(!(distance %in% c('glm','random forest'))){stop("Method specified not available,please choose one between 'glm' and 'random forest'")}
  if(distance=='random forest' & is.null(tune)){stop("Please express preference for tuning when method == 'random forest'")}
  if(!(type %in% c('stabilised','unstabilised','overlap'))){stop("Type should be either stabilised or unstabilised")}
  # extract variables
  treatment_var=all.vars(formula)[1]
  treatment_quo=sym(treatment_var)
  # calculate probability  based on distance pointed out
  if (distance=='glm'){
    data$distance=predict(glm(formula,data=data,family='binomial'),type="response")
  } else if (distance=='random forest'){
    if (tune==TRUE){
      data[[treatment_var]]=as.factor(data[[treatment_var]])
      tune_result=tuneRF(x=data[, all.vars(formula)[-1]], 
                           y=data[[treatment_var]],
                           ntreeTry=500,
                           stepFactor=1.5,
                           improve=0.01,
                           trace=FALSE,
                           plot=FALSE)
      best_mtry=tune_result[which.min(tune_result[, 2]), 1]
      rf_model=randomForest(formula, data=data, ntree=500, mtry=best_mtry)
      data$distance = rf_model$votes[, 2]
      data[[treatment_var]] = as.integer(as.character(data[[treatment_var]]))
    } else {
      data[[treatment_var]]=as.factor(data[[treatment_var]])
      rf_model=randomForest(formula, data=data, ntree=500)
      data$distance = rf_model$votes[, 2]
      data[[treatment_var]]=as.integer(as.character(data[[treatment_var]]))
    }
  }
  # clip to avoid errors
  data$distance = pmax(pmin(data$distance, 1-1e-6), 1e-6)
  # calculate proportions of treated
  pt=sum(data[[treatment_var]]==1)/nrow(data)
  # calculate weights
  if(type=='stabilised'){
    if(estimand=='ATE'){
      data = data %>% mutate(w = ifelse(!!treatment_quo==1, pt/distance, (1-pt)/(1-distance)))
    } else if(estimand=='ATT'){
      # controls: odds weights rescaled by (1-pt)/pt so they sum to ~n_ctrl
      data = data %>% mutate(w = ifelse(!!treatment_quo==1, 1, (1-pt)*distance/(pt*(1-distance))))
    }
  } else if(type=='unstabilised'){
    if(estimand=='ATE'){
      data = data %>% mutate(w = ifelse(!!treatment_quo==1, 1/distance, 1/(1-distance)))
    } else if(estimand=='ATT'){
      data = data %>% mutate(w = ifelse(!!treatment_quo==1, 1, distance/(1-distance)))
    }
  } else if(type=='overlap'){
    data$w = data$distance * (1 - data$distance)
  }
  # truncate if asked
  ## calculate percentile
  if (!is.null(truncate)){
    pl=quantile(data$w,truncate)
    pu=quantile(data$w,1-truncate)
    ## apply truncation
    data=data %>% 
      mutate(w=case_when(w<pl~pl,
                         w>pu~pu,
                         T~w))
  }
  # return
  return(list(data=data,weights=data$w,formula=formula,type=type,truncate=truncate,distance=distance,tune=tune))
}