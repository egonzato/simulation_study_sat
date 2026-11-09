# libraries
pacman::p_load(dagitty, ggdag, ggplot2,dplyr)
## survival
# set relationships
dag = dagitty("dag {
  sex -> aoo
  sex -> trt
  sex -> S
  aoo -> trt
  aoo -> S
  trt -> S
  yoo -> trt
  yoo -> S
}")
# coordinates
coordinates(dag) = list(x=c(sex = 0.10, aoo = 0.40, trt = 0.25, yoo = 0.75, S = 0.90),
                        y=c(sex = 0.75, aoo = 0.50, trt = 0.90, yoo = 0.90, S = 0.55))
# labels (must match node name order)
labs = c(sex = "sex", aoo = "aoo", trt = "trt", yoo = "yoo", S = "S")
# turn to tidy object
dag_tidy = tidy_dagitty(dag)
# plot
p=ggplot(dag_tidy,aes(x=x,y=y,xend=xend,yend=yend))+
  geom_dag_edges(edge_color="black",arrow_directed=grid::arrow(length=unit(0.2,"cm")))+
  geom_dag_point(aes(color=name),size=14)+
  geom_dag_label(aes(label=name),fill="white",color="black",size=5,nudge_y=0.05)+
  scale_color_manual(name=NULL,
                     #values=c(confounders="#4daf4a",treatment="#377eb8",outcome="#e41a1c"),
                     values = c(S="#e41a1c",
                                trt="#377eb8",
                                yoo="#4daf4a",
                                sex="#4daf4a",
                                aoo="#4daf4a"))+
  theme_dag()+
  guides(col='none')+
  theme(legend.position="right")
# save plot
ggsave("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Figures\\paper\\dag_surv.png",
       p,
       width = 6, height = 6, units = "in",
       dpi = 300)
## continuous
# set relationships
dag = dagitty("dag {
  sex -> aoo
  sex -> trt
  sex -> HUI
  aoo -> trt
  aoo -> HUI
  trt -> HUI
  yoo -> trt
  yoo -> HUI
  t -> HUI
}")
# coordinates
coordinates(dag) = list(x=c(sex = 0.10, aoo = 0.40, trt = 0.25, yoo = 0.75,t = 0.75,HUI = 0.90),
                        y = c(sex = 0.75, aoo = 0.50, trt = 0.90, yoo = 0.90,t = 0.25,HUI = 0.55))
# labels (must match node name order)
labs = c(sex = "sex", aoo = "aoo", trt = "trt", yoo = "yoo", t = "t", HUI = "HUI")
# define labels
dag_tidy=tidy_dagitty(dag) %>%
  mutate(group = case_when(name %in% c("sex", "aoo", "yoo") ~ "confounders",
                           name == "trt" ~ "treatment",
                           name == "HUI" ~ "outcome",
                           TRUE ~ "other"))
# plot dag
p=ggplot(dag_tidy,aes(x=x,y=y,xend=xend,yend=yend))+
  geom_dag_edges(edge_color="black",arrow_directed=grid::arrow(length=unit(0.2,"cm")))+
  geom_dag_point(aes(color=name),size=14)+
  geom_dag_label(aes(label=name),fill="white",color="black",size=5,nudge_y=0.05)+
  scale_color_manual(name=NULL,
                     #values=c(confounders="#4daf4a",treatment="#377eb8",outcome="#e41a1c"),
                     values = c(HUI="#e41a1c",
                                trt="#377eb8",
                                t="#4daf4a",
                                yoo="#4daf4a",
                                sex="#4daf4a",
                                aoo="#4daf4a"))+
  theme_dag()+
  guides(col='none')+
  theme(legend.position="right")
# save plot
ggsave("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Figures\\paper\\dag_cont.png",
       p,
       width = 6, height = 6, units = "in",
       dpi = 300)
