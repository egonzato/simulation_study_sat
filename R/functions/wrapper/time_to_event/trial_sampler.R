trial_sampler=function(trtsubpop=NULL,ctrlsubpop=NULL,ntrt=NULL,nctrl=NULL){
  # sample treated
  trt_sample=trtsubpop %>%
    slice_sample(n = ntrt)
  # sample untreated
  ctrl_sample=ctrlsubpop %>%
    slice_sample(n = nctrl)
  # bind together
  sampled=bind_rows(trt_sample, ctrl_sample)
  # return
  return(sampled)
}