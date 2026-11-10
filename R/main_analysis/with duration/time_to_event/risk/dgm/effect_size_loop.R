## -------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,lubridate,survival,survminer,ggsurvfit,rlang,ipw,ggplot2)


## -------------------------------------------------------------------------------------------------------------
function_dirs=c("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\wrapper\\time_to_event",
                "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\matching",
                "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\regression")
for (i in 1:length(function_dirs)){
  files=list.files(function_dirs[i], pattern = "\\.R$", full.names = TRUE)
sapply(files, source)
}


## -------------------------------------------------------------------------------------------------------------
options(scipen = 999)
set.seed(181299)
# define N
n=1e7
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


## -------------------------------------------------------------------------------------------------------------
set.seed(03111963)
n=1e7
df=data.frame(id=1:n,
              onsetdate=onsetdatev,
              ageonset=runif(n,0.5,6),
              event=1,
              u=runif(n,0,1),
              e=rnorm(n,0.5,0.02),
              sex=sample(c('male','female'),n,replace=T))
# disease duration
df=df %>% 
  mutate(calendar=year(onsetdate),
         calendar_c=calendar-mean(year(onsetdate)),
         lpcov=-0.243*ageonset+0.158*(sex=="male")-0.01*calendar_c+e,
         trialdate=as.Date('2021-05-31'))
# set parameters to calibrate mortality
pzero=0.5; taud=2*365.25
# calibrate given lpcov, pzero and taud
eta=df$lpcov
lambda_rexp=uniroot(function(lambda) {
  mean(1 - exp(-lambda * exp(eta) * taud / 30)) - pzero},
  interval = c(0, 1e15),
  tol = 1e-12)$root
# define trt effect with new 
df=df %>% 
  mutate(hazard=lambda_rexp * exp(lpcov),
         timedays=rexp(n,rate=hazard)*30,
         stopdate=timedays+onsetdate,
         ptrt_conf=1/(1+exp(-(+1 -0.93*ageonset + 2.04*I(sex=='male')))),
         trt_yoo=ifelse(stopdate>=trialdate,1,0),
         trt_conf=rbinom(n,1,ptrt_conf),
         trt=ifelse(trt_conf == 1 & trt_yoo == 1,1,0),
         duration=ifelse(trt==1,
                         trialdate-onsetdate,
                         round(timedays))) %>% 
  # treated need a waiting time (onset before the trial date); controls must be historical (onset before the trial date)
  filter(duration>0,
         !(trt==0 & onsetdate>=trialdate))
# histogram of probabilities
hist(df$ptrt_conf)
# number of treated
table(df$trt)


## -------------------------------------------------------------------------------------------------------------
# smd on ageonset and sex between treated and untreated
smd::smd(df$ageonset,df$trt)$estimate
smd::smd(df$sex,df$trt)$estimate


## -------------------------------------------------------------------------------------------------------------
target_population=df %>% filter(trt==1)
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
          interval = c(1e-15, 1),
          tol = 1e-12)$root
}


## -------------------------------------------------------------------------------------------------------------
targets = c(1, 0.75, 0.5, 0.25)
ks = sapply(targets, find_k)
names(ks) = targets


## -------------------------------------------------------------------------------------------------------------
# set path
folder="C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\main_analysis\\with duration\\time_to_event\\population_rr_"
# save multiple dataset
for (i in 1:length(targets)){
  # apply treatment
  dfi=df %>% 
  mutate(hazardtrt=hazard*ks[i],
         timedaystrt=ifelse(trt==1,
                            rexp(n(),hazardtrt),
                            NA_real_)*30,
         timetotal=ifelse(trt==1,timedaystrt+duration,timedays),
         start=as.Date(ifelse(trt==1,trialdate,onsetdate)),
         stop=as.Date(ifelse(trt==1,start+timedaystrt,start+timedays)),
         time=as.numeric(stop-start)) 
  # define dataset
  path=paste(folder,targets[i],".rds",sep="")
  # save dataset
  saveRDS(dfi,
          path)
}


## -------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\with duration\\time_to_event\\risk\\dgm\\effect_size_loop.qmd")

