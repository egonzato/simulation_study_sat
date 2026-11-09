bootstrap_sampler=function(data=NULL,treatment=NULL,ntrt=NULL,nctrl=NULL){
  # extract treatment
  trt_quo=enquo(treatment)
  # sample from treated
  boot_treated=data %>% 
    filter(!!trt_quo==1) %>% 
    slice_sample(n=ntrt,replace=TRUE)
  # sample from controls
  boot_controls=data %>% 
    filter(!!trt_quo==0) %>% 
    slice_sample(n=nctrl,replace=TRUE)
  # bind together
  bootstrap_sample=bind_rows(boot_treated,boot_controls)
  # give new ids
  bootstrap_sample$id=1:nrow(bootstrap_sample)
  #return
  return(bootstrap_sample)
}