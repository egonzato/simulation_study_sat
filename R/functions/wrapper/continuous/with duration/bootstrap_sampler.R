bootstrap_sampler=function(full_data=NULL,id=NULL,treatment=NULL,ntrt=NULL,nctrl=NULL){
  # extract variables
  id_quo=enquo(id); id_name=rlang::as_name(id_quo)
  treatment_quo=enquo(treatment)
  # define baseline dataset
  baseline=full_data %>% 
    group_by(!!id_quo) %>% 
    filter(row_number()==1) %>% 
    ungroup()
  # sample treated
  trt_sample=baseline %>%
    filter(!!treatment_quo==1) %>% 
    slice_sample(n = ntrt,replace=TRUE) %>% 
    mutate(id=!!id_quo) %>% 
    select(id)
  # sample untreated
  ctrl_sample=baseline %>%
    filter(!!treatment_quo==0) %>% 
    slice_sample(n=nctrl,replace=TRUE) %>% 
    mutate(id=!!id_quo) %>% 
    select(id)
  # bind together ids
  sampled_ids=bind_rows(trt_sample,ctrl_sample) %>% 
    mutate(boot_id=row_number()) %>% 
    select(id,boot_id)
  # merge back
  bootstrap_sample=full_data %>% 
    inner_join(sampled_ids,
               by=setNames("id",id_name),
               relationship="many-to-many") %>% 
    mutate(!!id_quo := boot_id) 
  #return
  return(bootstrap_sample)
}