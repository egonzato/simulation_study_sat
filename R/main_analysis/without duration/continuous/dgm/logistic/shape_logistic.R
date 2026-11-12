## -------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,tidyr,ggplot2,lubridate)


## -------------------------------------------------------------------------------------------------------
expit=function(x) 1/(1+exp(-x))
logit=function(p) log(p/(1-p))
# bounds of the HUI multiscore
L=-0.36; U=1
hui=function(eta) L+(U-L)*expit(eta)
# latent value for a given HUI (to choose a0 and k from interpretable values)
eta_of=function(h) logit((h-L)/(U-L))
# shape: HUI of an average patient at onset and after `t_ref` years
hui_onset=0.85; hui_ref=0; t_ref=5
a0=eta_of(hui_onset)
k=(eta_of(hui_onset)-eta_of(hui_ref))/t_ref
# between-patient SD (random intercept) and measurement error, both on the latent scale
sd_u0=0.5; sd_eps=0.3
# illustrative covariate effects and treatment effect (latent scale); the final values come from tune_logistic.qmd
b_aoo=0.05; b_male=-0.2; b_yoo=0.02; gamma=0.3
round(c(a0=a0,k=k,delta=gamma/k),3)


## -------------------------------------------------------------------------------------------------------
tibble(t=seq(0,10,by=0.05)) %>%
  mutate(`logistic (sensitivity)`=hui(a0-k*t),
         `linear (main)`=0.95-0.1*t) %>%
  pivot_longer(-t,names_to="model",values_to="HUI") %>%
  ggplot(aes(x=t,y=HUI,colour=model))+
  geom_hline(yintercept=c(L,U),linetype=2,colour="grey50")+
  geom_line(linewidth=1)+
  labs(x="Years since onset",y="HUI",colour=NULL)+
  theme_bw()+theme(legend.position="bottom")


## -------------------------------------------------------------------------------------------------------
source_dir="C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\wrapper\\continuous\\without duration\\"
source(paste0(source_dir,"rtruncexp.R"))
set.seed(3111963)
n=300
end_date=as.Date("2022-12-31")
covs=tibble(id=1:n,
            dateonset=end_date-floor(rtruncexp(n=n,rate=0.0002,max=as.numeric(end_date-as.Date("2005-01-01")))),
            aoo=runif(n,0,8),
            sex=sample(c("male","female"),n,replace=TRUE),
            u0=rnorm(n,0,sd_u0)) %>%
  mutate(yoo=year(dateonset))
visits=covs %>%
  crossing(time=seq(0,365.25*10,by=182.625)) %>%
  mutate(terr=ifelse(time==0,0,runif(n(),0,90)),
         timen=(time+terr)/365.25,
         yoo_c=yoo-mean(covs$yoo),
         eta0=a0+u0+b_aoo*(aoo-mean(covs$aoo))+b_male*((sex=="male")-mean(covs$sex=="male"))+b_yoo*yoo_c-k*timen,
         eps=rnorm(n(),0,sd_eps),
         hui_untrt=hui(eta0+eps))


## -------------------------------------------------------------------------------------------------------
set.seed(181299)
ids=sample(covs$id,15)
visits %>% filter(id %in% ids) %>%
  ggplot(aes(x=timen,y=hui_untrt,group=id,colour=factor(id)))+
  geom_hline(yintercept=c(L,U),linetype=2,colour="grey50")+
  geom_point(size=0.8)+geom_line()+
  guides(colour="none")+
  labs(x="Years since onset",y="HUI (untreated)")+
  theme_bw()


## -------------------------------------------------------------------------------------------------------
visits %>%
  mutate(years=cut(timen,c(0,1,2,4,6,8,Inf),right=FALSE,labels=c("0-1","1-2","2-4","4-6","6-8","8+"))) %>%
  group_by(years) %>%
  summarise(mean_HUI=round(mean(hui_untrt),2),
            near_ceiling=round(mean(hui_untrt>0.95),3),
            near_floor=round(mean(hui_untrt< -0.3),3),.groups="drop")


## -------------------------------------------------------------------------------------------------------
visits %>% filter(id %in% ids) %>%
  mutate(hui_trt=hui(eta0+gamma*timen+eps)) %>%
  pivot_longer(c(hui_untrt,hui_trt),names_to="arm",values_to="HUI") %>%
  mutate(arm=ifelse(arm=="hui_trt","treated","untreated")) %>%
  ggplot(aes(x=timen,y=HUI,group=interaction(id,arm),colour=factor(id),linetype=arm))+
  geom_hline(yintercept=c(L,U),linetype=2,colour="grey50")+
  geom_line()+
  scale_linetype_manual(values=c(treated="solid",untreated="dashed"))+
  guides(colour="none")+
  labs(x="Years since start of treatment (= onset)",y="HUI",linetype=NULL)+
  theme_bw()+theme(legend.position="bottom")


## -------------------------------------------------------------------------------------------------------
visits %>% filter(timen>=1.9,timen<=2.1) %>%
  mutate(hui_trt=hui(eta0+gamma*timen+eps),effect=hui_trt-hui_untrt) %>%
  ggplot(aes(x=hui_untrt,y=effect))+
  geom_point(alpha=0.5)+
  labs(x="Untreated HUI at 2 years",y="Treated - untreated HUI at 2 years",
       title=paste0("Effect at 2 years for gamma = ",gamma))+
  theme_bw()


## -------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\without duration\\continuous\\dgm\\logistic\\shape_logistic.qmd")

