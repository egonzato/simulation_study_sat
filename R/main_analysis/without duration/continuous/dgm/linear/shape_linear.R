## -------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,ggplot2,lubridate,ipw,emmeans,marginaleffects,
               lpSolve,randomForest,smd,
               optmatch,lme4,lmerTest,ggdag,RColorBrewer,MatchIt,lubridate)


## -------------------------------------------------------------------------------------------------------
function_dirs=c("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\wrapper\\time_to_event",
                "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\matching",
                "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\regression")
for (i in 1:length(function_dirs)){
  files=list.files(function_dirs[i], pattern = "\\.R$", full.names = TRUE)
  sapply(files, source)
}


## -------------------------------------------------------------------------------------------------------
# define coordinates in the space
coord_dag=list(x = c(hui=3,aoo=0.8,yoo=3,t=3,sex=0.8,trt=1.5),
               y = c(hui=0,aoo=1,yoo=1.5,t=-1,sex=0,trt=1.5))
# define colors
col=brewer.pal(1, 'Set1')
red=col[1]
blue=col[2]
green=col[3]
yellow=col[6]
# relationship between vars
model=dagify(hui~t+aoo+sex+yoo+trt,
             aoo~sex,
             trt~aoo+yoo+sex,
             coords = coord_dag,
             outcome='hui',
             exposure='trt',
             labels = c("aoo"="aoo", 
                        "yoo"="yoo",
                        "hui"="hui",
                        "trt"="trt",
                        "sex" = "sex",
                        "t"="t")) %>% 
  tidy_dagitty() %>%
  mutate(Variable = case_when(name == 'trt' ~ 'Exposure',
                              name == 'hui' ~ 'Outcome',
                              T~ "Observed"))
# plot DAG
ggdag(model, text = T, text_size = 3, text_col = 'white', edge_type = "link_arc",node_size=12)+
  #geom_dag_point(color = "grey", show.legend = T)+
 # ggdag_paths(model,shadow = T)+
  geom_dag_point(aes(colour = Variable))+
  geom_dag_label_repel(aes(label = label), nudge_x = 2,nudge_y = 5) +
  scale_color_manual(values=c(blue,green,red,'grey'))+
  theme_dag()+
  guides(color = guide_legend(override.aes = list(size = 2)))


## -------------------------------------------------------------------------------------------------------
set.seed(3111963)
# define N
n=1e4
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
# see distribution
hist(as.numeric(onsetdatev),
     breaks = 100,
     main = "Time before end date",
     xlab = "Date of onset (numeric)")
abline(v=as.numeric(as.Date('2021-05-31')),col="red",lty=2,lwd=2)


## -------------------------------------------------------------------------------------------------------
set.seed(181299)
covs=data.frame(id=1:n,
              aoo=runif(n,0,8),
              event=1,
              dateonset=as.Date(onsetdatev),
              u=runif(n,0,1),
              e=rnorm(n,0.5,0.02),
              sex=sample(c('male','female'),n,replace=T))


## -------------------------------------------------------------------------------------------------------
set.seed(181299)
t=seq(0,2190,by=182.625)
id=1:n
# expand time
tid=expand.grid(id=id,time=t)
# merge with covariates
df=merge(covs,tid,by='id') %>% 
  arrange(id,time) %>% 
  mutate(terr=ifelse(time==0,0,runif(n,0,+90)),
         date=dateonset+time+terr,
         yoo=year(dateonset),
         yoo_c=yoo-mean(yoo))


## -------------------------------------------------------------------------------------------------------
b0=0.95-(0.018*mean(df$sex=="male")+0.042*mean(df$aoo)+0.0002*mean(df$yoo))


## -------------------------------------------------------------------------------------------------------
set.seed(181299)
df=df %>% 
  group_by(id) %>% 
  mutate(row=row_number(),
         u0 = rnorm(1,0,0.07),         
         epsilon = rnorm(n(),0,0.08))


## -------------------------------------------------------------------------------------------------------
df=df %>% 
  mutate(timen=as.numeric(date-dateonset)/365.25,
         hui_untrt = (b0+u0) + 0.03*aoo + 0.01*(sex=="male") + 0.0002*yoo -0.1*timen+epsilon)


## -------------------------------------------------------------------------------------------------------
# random draw
ids=sample(unique(df$id),15)
# plot those random individuals
df %>% 
  filter(id %in% ids) %>% 
  ggplot(.,aes(x=timen,y=hui_untrt,group=id))+
  geom_point()+
  geom_line()+
  labs(x='Disease duration',y='HUI (quadratic)')+
  geom_hline(yintercept=-0.5,lty=2)+
  theme_bw()


## -------------------------------------------------------------------------------------------------------
# random draw
ids=sample(unique(df$id),15)
# plot those random individuals
df %>% 
  filter(id %in% ids) %>% 
  ggplot(.,aes(x=date,y=hui_untrt,group=id))+
  geom_point()+
  geom_line()+
  labs(x='Disease duration',y='HUI (linear)')+
  geom_hline(yintercept=-0.5,lty=2)+
  theme_bw()


## -------------------------------------------------------------------------------------------------------
# simulate probability of treatment
df=df %>% 
  group_by(id) %>% 
  mutate(trialdate=as.Date(round(runif(1,15000,17000))),
         trt_yoo=ifelse(dateonset>=trialdate,1,0),
         ptrt_conf=1/(1+exp(-(+2 - 0.55*aoo + 2*I(sex=='male')))),
         trt_conf=rbinom(1,1,ptrt_conf),
         trt=ifelse(trt_conf == 1 & trt_yoo == 1,1,0))
# see distribution of covariates
hist(df$ptrt_conf,col='white',
     main='Distribution of P(T|X)',
     xlab='Probability')
range(df$ptrt_conf)


## -------------------------------------------------------------------------------------------------------
df=df %>% 
  group_by(id) %>% 
  mutate(trt_overall=ifelse(sum(trt)>=1,1,0))


## -------------------------------------------------------------------------------------------------------
ids=df %>%
  group_by(id) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  mutate(flag=ifelse(trt_overall==0 & trt_yoo==1,1,0)) %>% 
  filter(flag==0) %>% 
  select(id)
# now select
df=df %>% 
  filter(id %in% ids$id)


## -------------------------------------------------------------------------------------------------------
# create baseline
baseline=df %>% 
  group_by(id) %>% 
  filter(row_number()==1) %>% 
  ungroup()
# smd
smd(baseline$aoo,baseline$trt_overall)
smd(baseline$sex,baseline$trt_overall)


## -------------------------------------------------------------------------------------------------------
df=df %>% 
  arrange(desc(trt_overall),id,date) %>% 
  group_by(id) %>%
  mutate(lastvisit=max(date),
         starttrial=as.Date(ifelse(any(trt == 1),
                              min(date[trt == 1], na.rm = TRUE),
                              as.Date(NA)))) %>%
  ungroup()


## -------------------------------------------------------------------------------------------------------
df=df %>% 
  arrange(id, date) %>% 
  group_by(id)                       


## -------------------------------------------------------------------------------------------------------
# loop to redefine dates
df=df %>% 
  group_by(id) %>% 
  mutate(date=if_else(trt==1, starttrial+90*(cumsum(trt)-1), date)) %>% 
  ungroup()


## -------------------------------------------------------------------------------------------------------
df=df %>% 
  mutate(timen=as.numeric(date-dateonset)/365.25,
         timen_trt=ifelse(trt==0,
                          0,
                          as.numeric(date-starttrial)/365.25),
         hui_trt=hui_untrt+0*trt+0.05*trt*timen_trt) %>% 
  group_by(id,trt) %>% 
  filter(!(trt==1 & !(row_number() %in% c(1,2,3,5,7,9))))


## -------------------------------------------------------------------------------------------------------
set.seed(181299)
df=df %>% 
  mutate(pmiss=ifelse(trt==1 | time==0 ,0,0.4),
         miss=rbinom(n(),1,pmiss)) %>% 
  group_by(id) %>% 
  mutate(flag=ifelse(trt_overall==1 & sum(trt)<6,1,0))%>%
  ungroup() %>% 
  filter(flag==0) %>% 
  filter(timen_trt<=2)


## -------------------------------------------------------------------------------------------------------
# plot
df %>% 
  #filter(year(date)>2013 & trt_overall==0) %>% 
  ggplot(aes(x = date, y = hui_trt,col=factor(trt),group=id)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values = c('#FF3300','#0099FF'))+
  labs(x = 'Disease duration', y = 'HUI (linear)') +
  geom_hline(yintercept = -0.5, lty = 2) +
  theme_bw()


## -------------------------------------------------------------------------------------------------------
# sample one id
singleid=sample(unique(df$id[which(df$trt == 1)]),1)
# plot
df %>% 
  filter(id %in% singleid) %>%
  filter(miss==0) %>% 
  ggplot() +
  geom_point(aes(x = date, y = hui_trt), col = '#0099FF') +
  geom_line(aes(x = date, y = hui_trt), col = '#0099FF') +
  geom_point(aes(x = date, y = hui_untrt), col = '#FF3300') +
  geom_line(aes(x = date, y = hui_untrt), col = '#FF3300') +
  labs(x = 'Disease duration', y = 'HUI (linear)') +
  geom_hline(yintercept = -0.5, lty = 2) +
  theme_bw()


## -------------------------------------------------------------------------------------------------------
writexl::write_xlsx(df,"C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_hui_like_without_dur.xlsx")


## -------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\without duration\\continuous\\dgm\\linear\\shape_linear.qmd")

