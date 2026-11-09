regression_call=function(data,id=NULL,time="timen",methods){
  # extract id
  id_enquo=enquo(id)
  time_var=rlang::sym(time)
  # baseline
  baseline=data %>% 
    group_by(!!id_enquo) %>% 
    filter(row_number()==1) %>% 
    ungroup()
  # branch based on methods
  if (methods=="raw"){
    dataset=data
  } else if (methods=="iptw, unstabilised"){
    # get weights
    weights=iptw(data=baseline,type='unstabilised',estimand='ATT',
                 formula=trt~aoo+sex)
    # merge with full dataset
    dataw=weights$data %>% 
      select(id,w)
    # merge dataset
    dataset=merge(data,dataw,by='id')
  }else if (methods=="iptw, stabilised"){
    weights=iptw(data=baseline,type='stabilised',estimand='ATT',
                 formula=trt~aoo+sex)
    # merge with full dataset
    dataw=weights$data %>% 
      select(id,w)
    # merge dataset
    dataset=merge(data,dataw,by='id')
  }else if (methods=="iptw, stabilised, truncated 0.01"){
    weights=iptw(data=baseline,type='stabilised',estimand='ATT',truncate=0.01,
                 formula=trt~aoo+sex)
    # merge with full dataset
    dataw=weights$data %>% 
      select(id,w)
    # merge dataset
    dataset=merge(data,dataw,by='id')
  }else if (methods=="iptw, stabilised, truncated 0.1"){
    weights=iptw(data=baseline,type='stabilised',estimand='ATT',truncate=0.1,
                 formula=trt~aoo+sex)
    # merge with full dataset
    dataw=weights$data %>% 
      select(id,w)
    # merge dataset
    dataset=merge(data,dataw,by='id')
  } else if (methods=="gcomp"){
    dataset=data
  }
  dataset %>% 
    group_by(!!id_enquo) %>% 
    arrange(!!time_var) %>% 
    mutate(hui_baseline=first(hui_linear_trt)) %>% 
    filter(row_number()!=1) %>% 
    ungroup()
}