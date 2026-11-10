## ----------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggh4x)


## ----------------------------------------------------------------------------------------------------------------------------------------------
results=readRDS("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\main_analysis\\without duration\\time_to_event\\results_population_rr_single_boot_bias_0.25_calibrate.rds")


## ----------------------------------------------------------------------------------------------------------------------------------------------
results_filtered=results %>% 
  filter(is.na(message) & result<100) %>% 
  mutate(truth=0.25,
         se=1,
         lower=1.96*se-result,
         upper=1.96*se+result,
         method2=paste(method,as.character(ratio),sep=""),
         z=abs((truth-result)/se))


## ----------------------------------------------------------------------------------------------------------------------------------------------
psm_p_greedy_caliper_011=results %>% 
  filter(method=="psm_p_greedy_caliper_01" & ratio == 1)


## ----------------------------------------------------------------------------------------------------------------------------------------------
exact=results %>% 
  filter(grepl("exact",method) & grepl("optimal",method))


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics=results_filtered %>% 
  mutate(count=ifelse(result>lower & result<upper,1,0)) %>%
  group_by(method2,ntrt,nctrl,truth,ratio) %>%
  summarise(n_reps=n(),
            bias=mean(result - truth),
            mcse_bias= sd(result) / sqrt(n_reps),
            mse=mean((truth - result)^2),
            mcse_mse=sqrt(sum(((truth - result)^2 - mse)^2) / (n_reps*(n_reps-1))),
            coverage=mean(count),
            mcse_cov=sqrt(coverage*(1-coverage)/n_reps),
            .groups = "drop")


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics=metrics %>%
  left_join(results_filtered %>% group_by(method2, nctrl,ntrt,ratio) %>%
              summarise(empirical_se = sd(result), model_se = mean(se),
                        .groups="drop"),
    by = c("method2","nctrl","ntrt")) %>%
  mutate(se_ratio = model_se / empirical_se) 


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------
# # reshape SMDs into long format
# smd_long=results_filtered %>%
#   select(method, ntrt, nctrl, smdageb, smdsexb, smdagea, smdsexa,ratio) %>%
#   pivot_longer(cols = c(smdageb, smdsexb, smdagea, smdsexa),
#                names_to = "var", values_to = "smd") %>%
#   mutate(covariate = case_when(grepl("age", var) ~ "Age of onset",
#                            grepl("sex", var) ~ "Sex"),
#          timing = case_when(grepl("b$", var) ~ "Before matching",
#                         grepl("a$", var) ~ "After matching")) %>%
#   select(-var)
# # summarise mean |SMD| across replicates
# smd_summary = smd_long %>%
#   group_by(method, ntrt, nctrl,ratio, covariate, timing) %>%
#   summarise(mean_abs_smd = mean(abs(smd), na.rm=TRUE),
#             se_smd = sd(abs(smd), na.rm=TRUE)/sqrt(n()),
#             .groups="drop")


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
  ggplot(aes(x=ntrt, y=bias, col=method2)) +
  geom_point(position=position_dodge(width=2)) +
  # geom_errorbar(aes(ymin=bias-1.96*mcse_bias, 
  #                   ymax=bias+1.96*mcse_bias),
  #               position=position_dodge(width=2),
  #               width=2) +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~nctrl) +
  #lims(y=c(0,0.6))+
  #ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Sample size of single treated arm', y='Bias', col='Matching method') +
  theme_bw() + theme(legend.position='bottom')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
   filter(method2 %in% c("rawNA","gcompNA",
                         "iptw, stabilised, truncated 0.01NA","iptw, stabilised, truncated 0.1NA",
                         "iptw, stabilisedNA","iptw, unstabilisedNA")) %>%
  #filter(bias<3) %>% 
  #filter(!(method2 %in% c("exact_optimal","psm_p_greedy_noc"))) %>% 
  ggplot(aes(x=ntrt, y=bias, col=method2)) +
  geom_point(position=position_dodge(width=2)) +
  # geom_errorbar(aes(ymin=bias-1.96*mcse_bias, 
  #                   ymax=bias+1.96*mcse_bias),
  #               position=position_dodge(width=2),
  #               width=2) +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~nctrl) +
  #lims(y=c(0,0.6))+
  #ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Sample size of single treated arm', y='Bias', col='Matching method') +
  theme_bw() + theme(legend.position='bottom')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
   filter(!(method2 %in% c("rawNA","gcompNA",
                           "iptw, stabilised, truncated 0.01NA","iptw, stabilised, truncated 0.1NA","iptw, stabilisedNA","iptw, unstabilisedNA"))) %>%
  mutate(distance=case_when(grepl("psm_p",method2)~"probability",
                            grepl("psm_m",method2)~"mahalanobis",
                            grepl("exact",method2)~"exact",
                            T~"none"),
         search=case_when(grepl("greedy",method2)~"greedy",
                            grepl("optimal",method2)~"optimal",
                            T~"none"),
         caliper=case_when(grepl("caliper",method2)~"yes",
                           grepl("noc",method2)~"no",
                           T~"none")) %>% 
  #filter(bias<3) %>% 
  #filter(!(method2 %in% c("exact_optimal","psm_p_greedy_noc"))) %>% 
  ggplot(aes(x=ntrt, y=bias, col=method2)) +
  geom_point(position=position_dodge(width=2)) +
  # geom_errorbar(aes(ymin=bias-1.96*mcse_bias, 
  #                   ymax=bias+1.96*mcse_bias),
  #               position=position_dodge(width=2),
  #               width=2) +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~nctrl+distance) +
  #lims(y=c(0,0.6))+
  #ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Sample size of single treated arm', y='Bias', col='Matching method') +
  theme_bw() + theme(legend.position='bottom')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
   filter(method2 %in% c("psm_p_optimal1","psm_p_optimal2","psm_p_optimal3",
                         "psm_p_greedy_noc1","psm_p_greedy_noc2","psm_p_greedy_noc3",
                         "psm_p_greedy_caliper_011","psm_p_greedy_caliper_012","psm_p_greedy_caliper_013",
                         "psm_p_greedy_caliper_0051","psm_p_greedy_caliper_0052","psm_p_greedy_caliper_0053",
                         "rawNA")) %>%
  mutate(distance=case_when(grepl("psm_p",method2)~"probability",
                            grepl("psm_m",method2)~"mahalanobis",
                            grepl("exact",method2)~"exact",
                            T~"none"),
         search=case_when(grepl("greedy",method2)~"greedy",
                          grepl("optimal",method2)~"optimal",
                          T~"none"),
         caliper=case_when(grepl("caliper",method2)~"yes",
                           grepl("noc",method2)~"no",
                           T~"none")) %>%  
  ggplot(aes(x=ntrt, y=bias, col=method2)) +
  geom_point(position=position_dodge(width=2)) +
  # geom_errorbar(aes(ymin=bias-1.96*mcse_bias, 
  #                   ymax=bias+1.96*mcse_bias),
  #               position=position_dodge(width=2),
  #               width=2) +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~nctrl) +
  #lims(y=c(0,0.6))+
  #ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Sample size of single treated arm', y='Bias', col='Matching method') +
  theme_bw() + theme(legend.position='bottom')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
   filter(method2 %in% c("psm_m_optimal1","psm_m_optimal2","psm_m_optimal3",
                         "psm_m_greedy_noc1","psm_m_greedy_noc2","psm_m_greedy_noc3",
                         "psm_m_greedy_caliper_011","psm_m_greedy_caliper_012","psm_m_greedy_caliper_013",
                         "psm_m_greedy_caliper_0051","psm_m_greedy_caliper_0052","psm_m_greedy_caliper_0053",
                         "rawNA")) %>%
  #filter(bias<3) %>% 
  #filter(!(method2 %in% c("exact_optimal","psm_p_greedy_noc"))) %>% 
  ggplot(aes(x=ntrt, y=bias, col=method2)) +
  geom_point(position=position_dodge(width=2)) +
  # geom_errorbar(aes(ymin=bias-1.96*mcse_bias, 
  #                   ymax=bias+1.96*mcse_bias),
  #               position=position_dodge(width=2),
  #               width=2) +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~nctrl) +
  #lims(y=c(0,0.6))+
  #ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Sample size of single treated arm', y='Bias', col='Matching method') +
  theme_bw() + theme(legend.position='bottom')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
  ggplot(data=.,aes(x=ntrt_pool,y=mse,group=method,col=method))+
  geom_point(position=position_dodge(width=2))+
  geom_errorbar(aes(ymin=mse-1.96*mcse_mse, 
                    ymax=mse+1.96*mcse_mse), 
                width=2,
                position=position_dodge(width=2))+
  #scale_x_continuous(breaks=c(50,75,100))+
  labs(x='Sample size of single treated arm',
       y='Mean Squared Error',
       col='Matching method')+
  facet_wrap(~nctrl_pool)+
  #ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  #scale_color_manual(values=c('#FF0000','#0066FF'))+
  theme_bw()+
  guides(fill='none')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>%
  ggplot(aes(x=empirical_se, y=model_se, col=method)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, lty=2) +
  ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Empirical SE', y='Model-based SE') +
  theme_bw()+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
metrics %>% 
  ggplot(data=.,aes(x=ntrt_pool,y=coverage,group=method,col=method))+
  geom_point(position=position_dodge(width=2))+
  geom_errorbar(aes(ymin=coverage-1.96*mcse_cov, 
                    ymax=coverage+1.96*mcse_cov), width=2,
                position=position_dodge(width=2))+
  #facet_grid(wrap2~wrap1)+
  facet_wrap(~nctrl_pool)+
  scale_y_continuous(breaks=seq(0.75,0.95,0.05))+
  labs(x='Sample size of single treated arm',
       y='Coverage',
       col='Matching method')+
  #ggh4x::facet_grid2(wrap2~wrap1, scales = "free", independent = "x")+
  #scale_color_manual(values=c('#FF0000','#0066FF'))+
  theme_bw()+
  guides(fill='none')+
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom')


## ----------------------------------------------------------------------------------------------------------------------------------------------
results %>% 
  ggplot(.,aes(x=method,y=beta))+
  geom_boxplot()+
  ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  theme_bw()+
  geom_hline(yintercept=log(0.6),lty=2)+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(face = "bold"))


## ----------------------------------------------------------------------------------------------------------------------------------------------
results %>% 
  ggplot(.,aes(x=beta,group=method,fill=method))+
  geom_density(alpha=0.7,col='black')+
  ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "x")+
  labs(x='Beta coefficient',y='Count')+
  theme_bw()+
  guides(fill='none')+
  geom_vline(xintercept=log(0.6),lty=2)+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(face = "bold"))


## ----------------------------------------------------------------------------------------------------------------------------------------------
results %>% 
  filter(method=="iptw") %>% 
  ggplot(.,aes(x=beta,y=se))+
  geom_point(alpha=0.7)+
  #facet_wrap(~method)+
  ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "all")+
  geom_vline(xintercept=log(0.6),lty=2 )+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(face = "bold"))


## ----------------------------------------------------------------------------------------------------------------------------------------------
results %>%
  filter(method=="iptw") %>% 
  mutate(count=ifelse(dfbetas>lower & dfbetas<upper,1,0)) %>% 
  ggplot(.,aes(y=z,x=beta,col=factor(count))) +
  #facet_grid(dfnctrl~dfntrt)+
  geom_pointrange(aes(xmin=lower,xmax=upper))+
  ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "all")+
  #facet_wrap(~method)+
  labs(col='Contains')+
  scale_color_manual(values=c('#FF0000','#0066FF'))+
  geom_vline(xintercept=log(0.6),lty=2)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(face = "bold"))


## ----------------------------------------------------------------------------------------------------------------------------------------------
smd_summary %>%
  filter(covariate!="Sex") %>% 
  mutate(timing = factor(timing, levels=c("Before matching","After matching"))) %>%
  ggplot(aes(x=timing, y=mean_abs_smd, 
             col=method, group=method)) +
  geom_point(position=position_dodge(width=0.3)) +
  geom_errorbar(aes(ymin=mean_abs_smd-1.96*se_smd, 
                    ymax=mean_abs_smd+1.96*se_smd),
                width=0.15, position=position_dodge(width=0.3)) +
  geom_hline(yintercept=0.1, lty=2, col="grey40") + 
  ggh4x::facet_grid2(nctrl_pool~ntrt_pool, scales = "free", independent = "all")+
  #facet_grid(covariate ~ ntrt_pool + nctrl_pool, labeller=label_both) +
  labs(x=NULL, y="Mean |SMD|", col="Method") +
  theme_bw() +
  theme(strip.background=element_rect(fill="white"), 
        strip.text=element_text(face="bold"),
        legend.position="bottom")


## ----------------------------------------------------------------------------------------------------------------------------------------------
smd_long %>%
  ggplot(aes(x=smd, fill=method)) +
  geom_density(alpha=0.6) +
  geom_vline(xintercept=0, lty=2) +
  ggh4x::facet_grid2(covariate~timing, scales = "free", independent = "all")+
  labs(x="SMD", y="Density", fill="Method") +
  theme_bw() +
  theme(strip.background=element_rect(fill="white"), 
        strip.text=element_text(face="bold"),
        legend.position="bottom")


## ----------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\without duration\\time_to_event\\risk\\reporting\\report_facet_bias_0.25_calibrate.qmd")

