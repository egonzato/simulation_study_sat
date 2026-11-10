## --------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,lubridate,survival,survminer,ggsurvfit,rlang,ipw,ggplot2)


## --------------------------------------------------------------------------------------------------------------------------------
set.seed(03111963)
n=1e6
df=data.frame(id=1:n,
              onsetdate=as.Date(round(runif(n,1500,20000)),origin='1970-01-01'),
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
         trialdate=as.Date('2015-05-23'),
         trt_yoo=ifelse(onsetdate>=trialdate,1,0),
         ptrt_conf=1/(1+exp(-(1.3*ageonset + 0.95*I(sex=='male')))),
         trt_conf=rbinom(n,1,ptrt_conf),
         trt=ifelse(trt_conf == 1 & trt_yoo == 1,1,0))
# number of treated
table(df$trt)


## --------------------------------------------------------------------------------------------------------------------------------
df=df %>% 
  mutate(hazardtrt=hazard*0.6,
         timedaystrt=ifelse(trt==1,
                            rexp(n(),hazardtrt),
                            NA_real_)*30,
         timetotal=ifelse(trt==1,timedaystrt,timedays))


## --------------------------------------------------------------------------------------------------------------------------------
writexl::write_xlsx(df,"C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_nodur.xlsx")


## --------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\without duration\\time_to_event\\HR\\dgm\\drift.qmd")

