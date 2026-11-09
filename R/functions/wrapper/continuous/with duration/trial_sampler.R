trial_sampler=function(full_data=NULL,id=NULL,treatment=NULL,
                      ntrt=NULL,nctrl=NULL){
  # extract variables
  id_quo=enquo(id)
  treatment_quo=enquo(treatment)
  # define baseline dataset
  baseline=full_data %>% 
    group_by(!!id_quo) %>% 
    filter(row_number()==1) %>% 
    ungroup()
  # sample treated
  trt_sample=baseline %>%
    filter(!!treatment_quo==1) %>% 
    slice_sample(n = ntrt) %>% 
    mutate(id=!!id_quo) %>% 
    select(id)
  # sample untreated
  ctrl_sample=baseline %>%
    filter(!!treatment_quo==0) %>% 
    slice_sample(n=nctrl) %>% 
    mutate(id=!!id_quo) %>% 
    select(id)
  # bind together
  sampled=bind_rows(trt_sample, ctrl_sample)
  # keep selected
  sampled_long=full_data %>% 
    filter(!!id_quo %in% sampled$id)
  # return
  return(sampled_long)
}