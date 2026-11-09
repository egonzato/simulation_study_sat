## --------------------------------------------------------------------------------
pacman::p_load(dplyr,lubridate,survival,survminer,ggsurvfit,rlang,ipw,ggplot2)


## --------------------------------------------------------------------------------
set.seed(03111963)
n=1e7
df=data.frame(id=1:n,
              onsetdate=as.Date(round(runif(n,12000,20000)),origin='1970-01-01'),
              ageonset=runif(n,0.5,6),
              event=1,
              u=runif(n,0,1),
              e=rnorm(n,0.5,0.02),
              sex=sample(c('male','female'),n,replace=T))
# disease duration
df=df %>% 
  mutate(calendar=year(onsetdate),
         linearpred=-0.004*ageonset+0.002*I(sex=='male')-0.001*calendar+e,
         hazard=exp(linearpred),
         timedays=rexp(n,hazard)*30,
         stopdate=timedays+onsetdate)
# treatment model
df=df %>% 
  mutate(trialdate=as.Date('2015-05-23'),
         trt_yoo=ifelse(stopdate>=trialdate,1,0),
         ptrt_conf=1/(1+exp(-(+ 1.7*ageonset + 1.2*I(sex=='male')))),
         trt_conf=rbinom(n,1,ptrt_conf),
         trt=ifelse(trt_conf == 1 & trt_yoo == 1,1,0))
# number of treated
table(df$trt)


## --------------------------------------------------------------------------------
target_population=df %>% mutate(trt==1)
# functions to fine tune
## compute risk ratio
compute_rr = function(k, t=730, pop=target_population){
  s0 = mean(exp(-pop$hazard * t/30))      
  s1 = mean(exp(-k * pop$hazard * t/30))
  risk0 = 1 - s0
  risk1 = 1 - s1
  risk1 / risk0
}
## find deflation factor
find_k = function(target, t=730, pop=target_population){
  if (target == 1) return(1)   
  uniroot(function(k) compute_rr(k, t, pop) - target,
          interval=c(1e-6, 1))$root
}


## --------------------------------------------------------------------------------
targets = c(1, 0.75, 0.5, 0.25)
ks = sapply(targets, find_k)
names(ks) = targets


## --------------------------------------------------------------------------------
# set path
folder="C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\with duration\\time_to_event\\population_rr_"
# save multiple dataset
for (i in 1:length(targets)){
  # apply treatment
  df=df %>% 
  mutate(hazardtrt=hazard*ks[i],
         timedaystrt=ifelse(trt==1,
                            rexp(n(),hazardtrt),
                            NA_real_)*30,
         duration=ifelse(trt==1,
                         trialdate-onsetdate,
                         round(timedays)),
         timetotal=ifelse(trt==1,timedaystrt+duration,timedays),
         start=as.Date(ifelse(trt==1,trialdate,onsetdate)),
         stop=as.Date(ifelse(trt==1,start+timedaystrt,start+timedays)),
         time=as.numeric(stop-start)) %>% 
  filter(duration>0)
  # define dataset
  path=paste(folder,targets[i],".xlsx",sep="")
  # save dataset
  saveRDS(df,
          path)
}


## ----eval=FALSE------------------------------------------------------------------
# # writexl::write_xlsx(df,"C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_dur.xlsx")
# # save dataset
# saveRDS(df,"C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_dur.xlsx")


## --------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\with duration\\time_to_event\\risk\\dgm\\effect_size_loop.qmd")

