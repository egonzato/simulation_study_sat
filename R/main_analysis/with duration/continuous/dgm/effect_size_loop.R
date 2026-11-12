## -----------------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,ggplot2,lubridate,ipw,emmeans,marginaleffects,
               lpSolve,randomForest,
               optmatch,lme4,lmerTest,ggdag,RColorBrewer,MatchIt,lubridate)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
base="C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\"
function_dirs=c(paste0(base,"Code\\paper\\R\\functions\\wrapper\\continuous\\with duration"),
                paste0(base,"Code\\paper\\R\\functions\\matching"),
                paste0(base,"Code\\paper\\R\\functions\\regression"))
for (i in 1:length(function_dirs)){
  files=list.files(function_dirs[i], pattern = "\\.R$", full.names = TRUE)
  sapply(files, source)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(3111963)
n=1e4
start_date=as.Date("2005-01-01")
end_date=as.Date("2022-12-31")
max_days=as.numeric(end_date - start_date)
backward_days=rtruncexp(n=n, rate=0.0002, max=max_days)
onsetdatev=end_date-floor(backward_days)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
covs=data.frame(id=1:n,
                dateonset=onsetdatev,
                aoo=runif(n,0,8),
                event=1,
                u=runif(n,0,1),
                e=rnorm(n,0.5,0.02),
                sex=sample(c('male','female'),n,replace=T))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(181299)
t=seq(0,365.25*10,by=182.625)
tid=expand.grid(id=1:n,time=t)
df=merge(covs,tid,by='id') %>%
  arrange(id,time) %>%
  mutate(terr=ifelse(time==0,0,runif(n(),0,+90)),
         date=dateonset+time+terr,
         yoo=year(dateonset),
         yoo_c=yoo-mean(yoo))
df=df %>%
  group_by(id) %>%
  mutate(row=row_number(),
         u0=rnorm(1,0,0.07),
         epsilon=rnorm(n(),0,0.08)) %>% ungroup()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
# population constants, computed once and outside any grouping
b0=0.75-(0.0093*mean(df$aoo)-0.0430*mean(df$sex=="male")+0.005*mean(df$yoo_c))

df_sim=df %>% 
  group_by(id) %>% 
  mutate(trialdate=as.Date("2021-05-31"),
         trt_yoo=ifelse(date>=trialdate,1,0),
         ptrt_conf=1/(1+exp(-(0 - 0.296*aoo + 1.30*I(sex=='male')))),
         trt_conf=rbinom(1,1,ptrt_conf),
         trt=ifelse(trt_conf == 1 & trt_yoo == 1,1,0),
         trt_overall=ifelse(sum(trt)>=1,1,0)) %>% 
  ungroup()
# selection as in tune.qmd:
# - treated: onset before the trial date and at least 9 visits after it (full protocol schedule up to 24 months)
# - controls: never eligible for treatment, onset before the trial date and at least 5 visits before it
df_sim=df_sim %>% 
  group_by(id) %>% 
  mutate(npost=sum(trt_yoo==1),
         npre=sum(trt_yoo==0),
         post=first(dateonset)>=first(trialdate),
         keep=(trt_conf==1 & npost>=9 & !post) | (trt_overall==0 & !post & npre>=5)) %>% 
  ungroup() %>% 
  filter(keep) %>% 
  # historical controls: drop their visits after the trial date
  filter(!(trt_overall==0 & trt_yoo==1)) %>% 
  select(-npost,-npre,-post,-keep)
# start of the trial and duration
df_sim=df_sim %>% 
  arrange(desc(trt_overall),id,date) %>% 
  group_by(id) %>%
  mutate(lastvisit=max(date),
         starttrial=as.Date(ifelse(any(trt == 1),
                                   min(date[trt == 1], na.rm = TRUE),
                                   as.Date(NA))),
         duration=ifelse(trt_overall==1,
                         as.numeric(starttrial-dateonset)/365.25,
                         as.numeric(lastvisit-dateonset)/365.25)) %>%
  ungroup()
# loop to redefine dates
df_sim=df_sim %>% 
  group_by(id) %>% 
  mutate(date=if_else(trt==1, starttrial+90*(cumsum(trt)-1), date)) %>% 
  ungroup()
# time (disease clock), untreated outcome, and selection of treated visits, AFTER the dates are final
df_sim=df_sim %>% 
  mutate(timen=as.numeric(date-dateonset)/365.25,
         timen_trt=ifelse(trt==0,0,as.numeric(date-starttrial)/365.25),
         hui_untrt=(b0+u0)+0.0093*aoo-0.0430*(sex=="male")+0.005*yoo_c-0.1*timen+epsilon) %>% 
  group_by(id,trt) %>% 
  filter(!(trt==1 & !(row_number() %in% c(1,2,3,5,7,9)))) %>% 
  ungroup()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(181299)
targets=c(0,0.025,0.075,0.15)
folder=paste0(base,"Data\\paper\\main_analysis\\with duration\\continuous\\population_diff_")
for (i in 1:length(targets)){
  # treatment effect accrues with time on treatment (timen_trt = 0 at trial entry)
  out=df_sim %>% mutate(hui_trt=hui_untrt+0*trt+targets[i]*trt*timen_trt)
  saveRDS(out, paste(folder,targets[i]*2,".rds",sep=""))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\with duration\\continuous\\dgm\\effect_size_loop.qmd")

