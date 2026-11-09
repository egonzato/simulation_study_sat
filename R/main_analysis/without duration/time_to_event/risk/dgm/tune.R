## -----------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr, lubridate)


## -----------------------------------------------------------------------------------------------------------------------------------
function_dirs = c("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\wrapper\\time_to_event")
for (i in seq_along(function_dirs)) {
  files = list.files(function_dirs[i], pattern = "\\.R$", full.names = TRUE)
  sapply(files, source)
}


## -----------------------------------------------------------------------------------------------------------------------------------
options(scipen = 999)
set.seed(3111963)
# define N
n=1e5
# define start and end of registry
start_date=as.Date("2005-01-01")
end_date=as.Date("2022-12-31")
max_days=as.numeric(end_date - start_date)
# simulate dates
backward_days=rtruncexp(n=n,
                        rate=0.0002,
                        max=max_days)
# calculate onset date
onsetdatev=end_date-floor(backward_days)


## -----------------------------------------------------------------------------------------------------------------------------------
set.seed(181299)
# baseline covs
covs=data.frame(id=1:n,
                dateonset=onsetdatev,
                aoo=runif(n,0,8),
                event=1,
                u=runif(n,0,1),
                e=rnorm(n,0.5,0.02),
                sex=sample(c('male','female'),n,replace=T)) %>% 
  mutate(male=as.integer(sex=="male"),
         trialdate=as.Date("2021-05-31"),
         post=dateonset>=trialdate)


## -----------------------------------------------------------------------------------------------------------------------------------
# visits as in the DGM: date = dateonset + time + terr, terr ~ U(0,jitter) for time>0
tt=seq(0,1095,by=182.625)          # same grid as the DGM; change it here if you change it there
jitter=90
set.seed(22101957)                 # different from the seeds above on purpose
nvis_post=integer(n); first_post=rep(NA_real_,n); lastdur=rep(NA_real_,n)
for (k in seq_along(tt)){
  dk=as.numeric(covs$dateonset)+tt[k]+(if (k==1) 0 else runif(n,0,jitter))
  pre_k=dk<as.numeric(covs$trialdate)
  nvis_post=nvis_post+!pre_k
  first_post=ifelse(is.na(first_post) & !pre_k,dk,first_post)
  lastdur=ifelse(pre_k,dk-as.numeric(covs$dateonset),lastdur)   # last visit before the trial date
}
covs$nvis_post=nvis_post                                          # visits on/after the trial date
covs$wait=(first_post-as.numeric(covs$dateonset))/365.25          # years onset -> entry (0 if onset after trial)
covs$lastdur=lastdur/365.25                                       # controls: years onset -> last pre-trial visit
rm(dk,pre_k,first_post,nvis_post,lastdur)


## -----------------------------------------------------------------------------------------------------------------------------------
# filter pre
prep_long = function(base, t_int, t_age, t_sex, wait=FALSE){
  d = base %>% 
    mutate(ptrt=plogis(t_int - t_age*aoo + t_sex*male),
           conf=u < ptrt)
  if (!wait){
    d = d %>% 
      mutate(trt =as.integer(conf & post)) %>% 
      filter(!(post & !conf))   
  } else {
    d = d %>% 
      mutate(trt=as.integer(conf & nvis_post>=1),          # treated: at least one visit after the trial date
             duration=ifelse(trt==1,wait,lastdur)) %>% 
      filter(trt==1 | !post)                               # controls need onset before the trial date
  }
  d
}
get_smds_long = function(base, t_int, t_age, t_sex, wait=FALSE){
  d = prep_long(base, t_int, t_age, t_sex, wait)
  abs(c(smd::smd(d$aoo, d$trt)$estimate, smd::smd(d$male, d$trt)$estimate))
}


## -----------------------------------------------------------------------------------------------------------------------------------
solve_treatment_long = function(base, target_smd_age, target_smd_sex, t_int=2, start=c(0.65, 2.2), wait=FALSE){
  obj=function(par) sum((get_smds_long(base, t_int, par[1], par[2], wait) - c(target_smd_age, target_smd_sex))^2)
  o=optim(start, obj, method="Nelder-Mead", control=list(reltol=1e-10))
  c(t_age=o$par[1], t_sex=o$par[2])
}


## -----------------------------------------------------------------------------------------------------------------------------------
# expected covariate differences after duration matching: for each treated wait w,
# the control is drawn uniformly among those with duration > w
matched_diffs = function(d){
  tr=d[d$trt==1,]; ct=d[d$trt==0,]
  o=order(ct$duration); dc=ct$duration[o]
  suf=function(x) c(rev(cumsum(rev(x))),0)
  Sa=suf(ct$aoo[o]); Sm=suf(ct$male[o]); N=nrow(ct)
  idx=findInterval(tr$duration,dc)          # controls with duration <= w
  n_el=N-idx; ok=n_el>0; i1=idx[ok]+1
  c(diff_aoo=mean(tr$aoo[ok])-mean(Sa[i1]/n_el[ok]),
    diff_male=mean(tr$male[ok])-mean(Sm[i1]/n_el[ok]),
    dropped=sum(!ok))
}

solve_outcome_long = function(d, target_bias, share_age=0.5, effect=0.10, wait=FALSE){
  if (!wait){
    diff_aoo =mean(d$aoo[d$trt==1])  - mean(d$aoo[d$trt==0])
    diff_male=mean(d$male[d$trt==1]) - mean(d$male[d$trt==0])
  } else {
    md=matched_diffs(d)
    diff_aoo=unname(md["diff_aoo"]); diff_male=unname(md["diff_male"])
  }
  c(b_age=share_age*target_bias*effect/diff_aoo,
    b_sex=(1-share_age)*target_bias*effect/diff_male,
    diff_aoo=diff_aoo, diff_male=diff_male)
}


## -----------------------------------------------------------------------------------------------------------------------------------
calibrate_long = function(base, target_smd_age, target_smd_sex, target_bias, t_int=-1, share_age=0.5, wait=FALSE){
  trt_par=solve_treatment_long(base, target_smd_age, target_smd_sex, t_int, wait=wait)
  d=prep_long(base, t_int, trt_par["t_age"], trt_par["t_sex"], wait)
  out_par=solve_outcome_long(d, target_bias, share_age, wait=wait)
  list(parameters=c(t_int=t_int, trt_par, out_par[c("b_age","b_sex")]),
       achieved=c(get_smds_long(base, t_int, trt_par["t_age"], trt_par["t_sex"], wait),
                  bias=unname(out_par["b_age"]*out_par["diff_aoo"] + out_par["b_sex"]*out_par["diff_male"])/0.10),
       n_treated=table(d$trt))
}

# case with duration, new selection rule
calibrate_long(covs, target_smd_age=0.5, target_smd_sex=0.5, target_bias=-0.20, t_int=0, wait=TRUE)

