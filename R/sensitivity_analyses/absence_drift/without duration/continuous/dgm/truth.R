## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,tidyr,tidyverse,ggplot2,lpSolve,MatchIt,
               optmatch,ipw,lme4,lmerTest,emmeans,boot,randomForest)


## ----warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
df=readxl::read_xlsx("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_hui_like_without_dur.xlsx")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trt=df %>% filter(trt_overall==1)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# treated
countertrt=trt %>%
  mutate(idc=paste("trt",id,sep=""),
         trt=1,
         time=timen,
         hui=hui_linear_trt) %>% 
  select(idc,id,trt,time,hui,yoo)
# untreated
counteruntrt=trt %>%
  mutate(idc=paste("untrt",id,sep=""),
         trt=0,
         time=timen,
         hui=hui_linear_untrt) %>% 
  select(idc,id,trt,time,hui,yoo)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
counter=rbind(countertrt,counteruntrt)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# sample four random ids
ids=sample(unique(counter$id),4)
# plot those ids
counter %>% 
  filter(id %in% ids) %>% 
  ggplot(.,aes(x=time,y=hui,group=idc,col=factor(trt)))+
  geom_point()+
  geom_line()+
  guides(col="none")+
  facet_wrap(~id)+
  scale_color_manual(values=c('#FF0000','#0066FF'))+
  theme_bw()+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# model
trial=lmer(hui~(1|idc)+trt*time,counter)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# marginalize
estimand=emmeans(trial, ~ trt | time, at = list(time=2))
# see contrast
pairs(estimand)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\without duration\\continuous\\dgm\\truth.qmd")

