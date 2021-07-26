library("RColorBrewer")
library("ggplot2")
library("cowplot")

source("preprocess_data.R")



# Model fitting ----------------------------------------------------------------
#
source("FitUVSDT.R")
# models <- c("GaussianEVSDT","GaussianUVSDT","GumbelEVSDT")
#
# for (model in models){
#
# model <- "ExGaussNormEVSDT"
#
# for(subjid in unique(testphase$id)){
# # for(experiment in c("SB2021_e1","SB2021_e2")){
#   fullsubj <- NULL
#   for(cond in c("A","B","C","D")){
#
#
#
#     data <- testphase %>% filter(id == subjid) %>% filter(condition==cond)
#
#
#
#   fit <- FitSDT(data = data, model = model, rep=20, freqdat = F) %>% mutate(condition = cond)
#
#     fullsubj <- fullsubj %>% bind_rows(fit)
#   }
#
#
#   saveRDS(fullsubj,file = paste0("Fits2/fit_seperateconditions_",model,"_",subjid,".rds"))
#
# }
#
# }


models <- c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

fits <- load_files("Fits2/","fit")

bestfit <- fits %>%
  #mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1)

#
#
# # Stability fit -------------------
#
# numbermodels <- 5
#
# minDelta <- fits %>%
#   filter(model %in% models) %>%
#   mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
#   mutate(AIC = 2 * npar + 2 * objective) %>%
#   group_by(model,id,condition) %>%
#   arrange(objective)  %>%
#   group_by(model,exp,id,condition) %>%
#   mutate(mindiffLL = outer(objective,objective, `-`)[2,1],
#          bestLL = objective[[1]],
#          diffDev = mindiffLL * 2,
#          bestDev = bestLL * 2) %>%
#   select(model,exp,id,condition,diffDev,bestDev) %>%
#   distinct() %>%
#   arrange(id) %>%
#   ungroup() %>%
#   mutate(diffObj = diffDev/bestDev * 100,
#          #trials = ntrials$n,
#          id = id) %>%
#   mutate(model = factor(model,
#                         levels =models)) %>%
#   mutate(modelnum = as.numeric(model)) %>%
#
#   mutate(modelposition = model) %>%
#   mutate(modelposition = factor(modelposition,
#                                 levels = models)) %>%
#   mutate(position = as.numeric(modelposition))
#
#
#
# # Deviance in Points figure:
#
# ## Text: data points per model not plotted in the graph
#
# PtsAbove10 <- minDelta %>% mutate(above = ifelse(diffDev > 10, 1, 0)) %>%
#   group_by(exp,condition,model,modelnum,position) %>%
#   summarize(notplotted = sum(above)/length(above) * 100) %>%
#   arrange(exp,condition,model) %>% filter(!notplotted %in% c(0,100) )
#
#
# ## Text: proportion of data sets underneath a threshold per model
#
# thresholdmodel <- minDelta %>%
#   group_by(exp,condition,model,modelnum,position) %>%
#   # in deviance terms...
#   summarize(belowDev05 = length(diffDev[diffDev < 0.5]),
#             belowDev1 = length(diffDev[diffDev < 1]),
#             belowDev2 = length(diffDev[diffDev < 2]),
#             belowDev10 = length(diffDev[diffDev < 10])) %>%
#   gather(key = "crit",value="value",-exp,-model,-condition,-modelnum,-position) %>%
#   group_by(exp,condition,model,modelnum,position,crit) %>%
#   summarize(propbelow = value/64) %>%
#   # y axis graph placement
#   mutate(yplace = case_when(crit == "belowDev05" ~ 0.25,
#                             crit == "belowDev1" ~ 0.75,
#                             crit == "belowDev2" ~ 1.5,
#                             crit == "belowDev10" ~ 6))
#
#
# ## Text: proportion of data sets underneath a treshold for any of the models
#
# thresholdtotal <- minDelta %>% ungroup() %>%
#   mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
#          belowDev1 = ifelse(diffDev < 1,1,0),
#          belowDev2 = ifelse(diffDev < 2,1,0),
#          belowDev10 = ifelse(diffDev < 10,1,0)) %>%
#   group_by(exp,id,condition) %>%
#   summarize(modelsbelowDev05 = sum(belowDev05),
#             modelsbelowDev1 = sum(belowDev1),
#             modelsbelowDev2 = sum(belowDev2),
#             modelsbelowDev10 = sum(belowDev10)) %>%
#   gather(key = "crit",value="value",-id,-exp,-condition) %>%
#   group_by(exp,condition,crit) %>%
#   summarize(propbelow = length(value[value == 3])/length(value)) %>%
#   mutate(yplace = case_when(crit == "modelsbelowDev05" ~ 0.25,
#                             crit == "modelsbelowDev1" ~ 0.75,
#                             crit == "modelsbelowDev2" ~ 1.5,
#                             crit == "modelsbelowDev10" ~ 6))
#
#
# thresholdtotal2 <- minDelta %>% ungroup() %>%
#   mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
#          belowDev1 = ifelse(diffDev < 1,1,0),
#          belowDev2 = ifelse(diffDev < 2,1,0),
#          belowDev10 = ifelse(diffDev < 10,1,0)) %>%
#   group_by(exp,id) %>%
#   summarize(modelsbelowDev05 = sum(belowDev05),
#             modelsbelowDev1 = sum(belowDev1),
#             modelsbelowDev2 = sum(belowDev2),
#             modelsbelowDev10 = sum(belowDev10)) %>%
#   gather(key = "crit",value="value",-id,-exp) %>%
#   group_by(exp,crit) %>%
#   summarize(propbelow = length(value[value == 3*4])/length(value)) %>%
#   mutate(yplace = case_when(crit == "modelsbelowDev05" ~ 0.25,
#                             crit == "modelsbelowDev1" ~ 0.75,
#                             crit == "modelsbelowDev2" ~ 1.5,
#                             crit == "modelsbelowDev10" ~ 6))
#
# # TotalPtsAbove10 <- (1 - thresholdtotal %>%
# #                       filter(crit == "modelsbelowDev10") %>%
# #                       .$propbelow) * 100
#
# TotalPtsAbove10 <- thresholdtotal %>%
#   filter(crit == "modelsbelowDev10") %>%
#   mutate(probpercent = (1 - propbelow) * 100)
#
# ggplotColours <- function(n = 10, h = c(0, 360) + 15){
#   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
#   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }
#
# colorin <- ggplotColours(n = 10)
#
# plotcolors <- colorin[unique(thresholdmodel$modelnum)]
#
# datastabpoints <- ggplot(minDelta %>% filter(diffDev <= 10),
#                          aes(x = modelnum,y = diffDev )) +
#
#   geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 0,ymax=0.5),
#             fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 0.5,ymax=1),
#             fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 1,ymax=2),
#             fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 2,ymax=10),
#             fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 10,ymax=10.5),
#             fill = "transparent",colour="grey")+
#
#   geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
#                 ymin = 0,ymax=0.5),fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
#                 ymin = 0.5,ymax=1),fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
#                 ymin = 1,ymax=2),fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
#                 ymin = 2,ymax=10),fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
#                 ymin = 10,ymax=10.5),fill = "transparent",colour="grey")+
#   geom_rect(aes(xmin = 0.75,xmax = max(position) + 1.5,
#                 ymin = 10.4,ymax=10.5),fill = "white",colour="white")+
#   facet_grid(exp~condition)+
#
#   geom_abline(aes(intercept=0.5,slope = 0),linetype="dotted",colour="grey") +
#   geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
#   geom_abline(aes(intercept=1,slope = 0),linetype="dotted",colour="grey") +
#   geom_abline(aes(intercept=2,slope = 0),linetype="dotted",colour="grey") +
#   geom_abline(aes(intercept=10,slope = 0),linetype="dotted",colour="grey") +
#
#   geom_text(data = thresholdmodel,
#             mapping = aes(x = position + 0.2, y = yplace,
#                           color = factor(modelnum),
#                           label = sub("^0+", "", paste(round(propbelow,2)))),
#             size = 3.5) +
#   geom_text(data = thresholdtotal,
#             mapping = aes(x = numbermodels + 1.25, y = yplace,
#                           label = sub("^0+", "", paste(round(propbelow,2)))),
#             color = "#A9A9A9",size = 3.5) +
#
#   geom_jitter(data = minDelta %>% filter(diffDev <10), aes(x= position - 0.2,
#                                                            color = factor(modelnum)),
#               size = 2,alpha=0.6,
#               width=0.15) +
#   geom_jitter(data = minDelta %>% filter(diffDev <10), aes(x = max(position) + 0.8,
#                                                            y = diffDev),
#               color = "#A9A9A9",size = 2,alpha=0.6,width=0.15) +
#   #scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
#   geom_text(data = PtsAbove10,
#             mapping = aes(x = position - 0.3, y = 10.5,
#                           label = paste0(round(notplotted),"%")),size = 2.5) +
#   geom_text(data = TotalPtsAbove10,
#             mapping = aes(x = 3 + .75, y = 10.5,
#                           label = ifelse(probpercent > 0,paste0(round(probpercent),"%"),NA)),
#                           size = 2.5) +
#
#
#   scale_y_continuous(name = expression(paste(Delta,"Deviance (in points)")),
#                      limits = c(0,10.5),
#                      breaks = c(c(0,0.5,1),seq(2,10,2))) +
#   scale_x_continuous(name = "Model", limits = c(.5,3 + 1.5),
#                      breaks = c(c(1:3),3 + 1),
#                      labels = c("EVSDT","UVSDT","Gumbel","Total")
#   ) +
#
#   theme(axis.text.x = element_text(size=11),#,angle=90,vjust=0.5
#         axis.text.y = element_text(size=12),
#         axis.title.y = element_text(size = 12),
#         axis.title.x = element_blank(),
#         legend.title = element_blank(),
#         axis.line = element_line(),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         strip.background = element_rect(fill = "transparent"),
#         strip.text = element_text(size = 12),
#         legend.position = "none",
#         legend.key = element_rect(fill = "transparent",colour="transparent"))
#

# Analyse fits -----------------------------------------------------------------

#
# exclLL <- unique(fits %>%
#   filter(model %in% models) %>%
#   mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
#   mutate(AIC = 2 * npar + 2 * objective) %>%
#   group_by(model,id,condition) %>%
#   arrange(objective)  %>%
#   group_by(model,exp,id,condition) %>%
#   dplyr::slice(1:2) %>%
#   mutate(mindiffLL =diff(objective)) %>%
#   select(exp,id,condition,model,objective,mindiffLL) %>%
#   dplyr::slice(1) %>%
#   filter(mindiffLL > 0.5) %>%
#   .$id)



  AICmeans <- bestfit %>%
    filter(model %in% c("ExGaussNormEVSDT","GumbelEVSDT","GaussianUVSDT")) %>%
    group_by(exp,condition,model) %>%
    summarize(means = mean(AIC)) %>%
    arrange(exp,condition,means)

  full<-NULL
  for(experiment in c("SB2021_e1","SB2021_e2")){
    for(cond in c("A","B","C","D")){

AICmeans1 <- AICmeans %>% filter(condition == cond) %>% filter(exp == experiment)

    #sum all aics
    AICtot <- bestfit %>%
      filter(condition == cond) %>% filter(exp == experiment) %>%
      filter(model %in% c("ExGaussNormEVSDT","GumbelEVSDT","GaussianUVSDT")) %>%
      mutate(deltaAIC = AIC - AICmeans1 %>% dplyr::slice(1) %>% .$means) %>%
      group_by(exp,condition,model) %>%
      mutate(model = factor(model,levels = c("ExGaussNormEVSDT","GumbelEVSDT","GaussianUVSDT"),
                            labels = c("ExGauss","Gumbel","UVSDT"))) %>%

      summarize(mdeltaAIC = mean(deltaAIC)) %>%
      mutate(condition = factor(condition,
                                levels = c("A","B","C","D"),
                                labels = c("High strength\nHigh variability",
                                           "High strength\nLow variability",
                                           "Low strength\nHigh variability",
                                           "Low strength\nLow variability"))) %>%
      mutate(exp = factor(exp,levels=c("SB2021_e1","SB2021_e2"),
                          labels=c("S&B (2021)\nExp 1","S&B (2021)\nExp 2")))


    full <- full %>% bind_rows(AICtot)

    }
  }



  indstrict_PW <- bestfit %>%
    filter(model %in% c("ExGaussNormEVSDT","GumbelEVSDT","GaussianUVSDT")) %>%
    group_by(exp,id,condition) %>%

    mutate(delAIC = AIC - min(AIC)) %>%
    #filter(id == subjid) %>% .$delAIC
    filter(delAIC == 0) %>% # only 1 model can win
    group_by(exp,condition,model) %>%
    summarize(numBest = length(delAIC)) %>%
    mutate(model = factor(model,levels = c("ExGaussNormEVSDT","GumbelEVSDT","GaussianUVSDT"),
                          labels = c("ExGauss","Gumbel","UVSDT"))) %>%


    mutate(condition = factor(condition,
                              levels = c("A","B","C","D"),
                              labels = c("High strength\nHigh variability\n",
                                         "High strength\nLow variability\n",
                                         "Low strength\nHigh variability\n",
                                         "Low strength\nLow variability\n"))) %>%
    mutate(exp = factor(exp,levels=c("SB2021_e1","SB2021_e2"),
                        labels=c("S&B (2021)\nExp 1","S&B (2021)\nExp 2"))) %>%
    complete(model,nesting(exp,condition), fill = list(0)) %>%
    group_by(exp,condition,model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    group_by(exp,condition) %>%
    mutate(wu =wu/sum(wu)) %>%
    arrange(exp,condition,model)


  ylim.prim <- c(0,5)
  ylim.sec <- c(0,1)

  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- b*(ylim.prim[1] - ylim.sec[1])


  devtools::install_github("tidyverse/ggplot2")
  library(ggplot2)

 p<-ggplot(full,
                           aes(y = mdeltaAIC,x = model )) +
    #geom_rect(aes(xmin=0, xmax=24, ymin=0, ymax=1), fill = "grey") +
    facet_grid(exp~condition)+
    geom_bar(stat = "identity",fill="#666666",width=1) +
    coord_flip(ylim = c(0,5.5)) +
    geom_point(aes(y = a +indstrict_PW$wu*b),
               color = "black",size=5,shape=21,fill="red") +


    #scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
    # scale_fill_gradientn(name = expression(paste(Delta)),
    #                      colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
    #                      limits=c(0,50))+
   # geom_text(aes(label=round(labeldelta), y = 52), hjust = 0.5,size=3)+
    scale_y_continuous(name = expression(paste(Delta," AIC")),
                       expand = c(0,0),
                       breaks = seq(0,5,1),
                       #labels = c(1:5),
                       sec.axis = sec_axis(~ (. - a)/b,
                                           name = "% Winning (AIC)",
                                           labels = c("0",".25",".50",".75","1"))) +

    scale_x_discrete(name = "Model") +

    theme(axis.text.x = element_text(size=12),
          axis.text.x.bottom = element_text(colour="#666666"),
          axis.ticks.x.bottom = element_line(color="#666666"),
          axis.title.x.bottom = element_text(colour="#666666"),
          axis.text.x.top= element_text(colour = "red"),
          axis.title.x.top = element_text(color= "red",vjust=1),
          axis.ticks.x.top = element_line(color="red"),
          axis.text.y.left = element_text(size=12,hjust=0),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill = "transparent"),
          strip.background.x = element_rect(fill = "white"),
          strip.background.y = element_rect(fill = "white"),
          strip.placement = "outside",
          strip.text = element_text(size=12))

 gp <- ggplotGrob(p)

 # gp$layout #helps you to understand the gtable object
 # gtable_show_layout(ggplotGrob(p)) #helps you to understand the gtable object


 gp$layout[gp$layout[['name']] == 'xlab-t' ,][['t']] <- 7#unique(t_strip)
 gp$layout[gp$layout[['name']] == 'xlab-t' ,][['b']] <- 7#unique(t_strip)

 gp$layout[grepl('strip', gp$layout[['name']]),][['t']] <- c(5,5,5,5,9,11)
 gp$layout[grepl('strip', gp$layout[['name']]),][['b']] <- c(5,5,5,5,9,11)


 library(grid)
 grid.newpage()
 empAIC <- ggdraw(gp)

 ggsave(paste0("Figures/SB2021Results.png"), units="cm", width=35, height=15, dpi=600)

#
# testdat <- bestfit %>% dplyr::select(exp,condition,model,objective,AIC) %>%
#   filter(model %in% c("GaussianUVSDT","ExGaussNormEVSDT","GumbelEVSDT")) %>%
#   mutate(objective = 2 * objective)
#
#
# maketable <- testdat %>%
#     pivot_longer(cols = c("objective","AIC"),names_to = "penalty",values_to="value") %>%
#     group_by(exp,condition,penalty,id) %>%
#     mutate(winning = ifelse(value == min(value),1,0)) %>%
#     group_by(exp,condition,penalty,model) %>%
#     summarize(num = sum(winning),
#               meanval = mean(value)) %>%
#     group_by(exp,condition,penalty) %>%
#     mutate(winprop = num/sum(num))
#
#
# winningprop <- maketable %>%
#   mutate(model = case_when(model == "GumbelEVSDT" ~ "sayA",
#                            model == "ExGaussNormEVSDT" ~ "sayB",
#                            TRUE~"sayC")) %>% ungroup() %>% dplyr::select(-num,-meanval) %>%
#   rename("value" = winprop) %>% mutate(type="Winprop")
#
# meanvalues <- maketable %>%
#   mutate(model = case_when(model == "GumbelEVSDT" ~ "sayA",
#                            model == "ExGaussNormEVSDT" ~ "sayB",
#                            TRUE~"sayC")) %>% ungroup() %>% dplyr::select(-num,-winprop) %>%
#   rename("value" = meanval) %>%
#   group_by(exp,condition,penalty) %>%
#   mutate(value = value - min(value)) %>%
#   mutate(type="Ave")
#
# total <- bind_rows(winningprop,meanvalues) %>%
#   mutate(penalty = factor(penalty,levels=c("objective","AIC"),
#                           labels=c("-2LL","AIC")))
#
#
# Exp1 <- ggplot(total %>% filter(exp == "SB2021_e1"),aes(y = type,x = model)) +
#   geom_tile(color="white",fill="white") +
#
#   geom_text(aes(label=ifelse(type=="Winprop" & value != 0,
#                              stringr::str_remove(round(value,2), "^0+"),
#                             paste(round(value,2)))),size=4)+
#   scale_x_discrete(position = "top",name="Recovered (% winning)",labels=c("Gumbel","EGNorm","UVSDT"))+
#   scale_y_discrete(position = "left",
#                    name="Generating",labels = c(expression(paste(Delta,frac(1,N),Sigma)),"% Winning"))+
#
#   facet_grid(penalty~condition)+
#   theme(axis.text.y = element_text(size=10),
#         axis.text.x = element_text(size=10),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = "transparent",colour="black"),
#         strip.placement = "outside",
#         legend.position = "none",
#         strip.background = element_rect(fill = "transparent",
#                                         colour="transparent"),
#         strip.text = element_text(size = 12))
#
# Exp2 <- ggplot(total %>% filter(exp == "SB2021_e2"),aes(y = type,x = model)) +
#   geom_tile(color="white",fill="white") +
#
#   geom_text(aes(label=ifelse(type=="Winprop" & value != 0,
#                              stringr::str_remove(round(value,2), "^0+"),
#                              paste(round(value,2)))),size=4)+
#   scale_x_discrete(position = "top",name="Recovered (% winning)",labels=c("Gumbel","EGNorm","UVSDT"))+
#   scale_y_discrete(position = "left",
#                    name="Generating",labels = c(expression(paste(Delta,frac(1,N),Sigma)),"% Winning"))+
#
#   facet_grid(penalty~condition)+
#   theme(axis.text.y = element_text(size=10),
#         axis.text.x = element_text(size=10),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = "transparent",colour="black"),
#         strip.placement = "outside",
#         legend.position = "none",
#         strip.background = element_rect(fill = "transparent",
#                                         colour="transparent"),
#         strip.text = element_text(size = 12))
#
# plot_grid(Exp1,Exp2,nrow=1,labels=c("Exp 1","Exp 2"))
# ggsave(paste0("EmpiricalData.png"), units="cm", width=45, height=10, dpi=600)
#
# differences <- bestfit %>%
#   filter(model %in% models_charvector) %>%
#   mutate(model = factor(model,levels = models)) %>%
#   group_by(exp,id,condition) %>%
#   mutate(delAIC = AIC - min(AIC)) %>%
#   mutate(colordelAIC = case_when(delAIC == 0 ~ 0,
#                                  delAIC > 0 & delAIC <= 1 ~ 1,
#                                  delAIC > 1 & delAIC <= 2 ~ 2,
#                                  delAIC > 2 & delAIC <= 5 ~ 3,
#                                  delAIC > 5 & delAIC <= 10 ~ 4,
#                                  delAIC > 10 & delAIC <= 25 ~ 5,
#                                  TRUE ~ 6)) %>%
#   mutate(condition = factor(condition,
#                             levels = c("A","B","C","D"),
#                             labels = c("High strength\nHigh variability",
#                                        "High strength\nLow variability",
#                                        "Low strength\nHigh variability",
#                                        "Low strength\nLow variability")))
#
# ggplot(differences %>% filter(exp == "SB2021_e2"),aes(y = delAIC,x=id)) +
#   geom_bar(stat = "identity")+
#   facet_grid(model~condition,scales="free_x")
#
#
# IndividualDelta <- ggplot(data = differences,
#                aes(y = model,x = id)) +
#   geom_tile(data = differences, aes(fill = factor(colordelAIC))) +
#
#
#   #geom_text(aes(label = round(delAIC)),size = 3.5) +
#   scale_fill_manual(values = c("#000000","#bd0026","#f03b20","#fd8d3c","#feb24c","#fed976","#ffffb2"),
#                     labels = c("= 0","> 0 & <= 1","> 1 & <= 2","> 2 & <= 5","> 5 & <= 10","> 10 & <= 25","> 25"),
#                     name = "Delta AIC")+
#   scale_x_discrete(position = "top",name = "id") +
#
#
#   facet_wrap(exp~condition,scales="free_x",nrow=2) +
#   theme(axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         axis.text.y = element_text(size=12),
#         axis.ticks = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = "transparent",colour="black"),
#         strip.placement = "outside",
#         #legend.position = "none",
#         strip.background = element_rect(fill = "transparent",
#                                         colour="transparent"),
#         strip.text = element_text(size = 12))
#
#
#
# plot_grid(datastabpoints,DeltaAIC,IndividualDelta,nrow=3,rel_heights = c(0.65,0.5,0.5),
#           labels = c("A","B","C"))
#
# ggsave("Test2.png", units="cm", width=30, height=30, dpi=600)
#
# # Investigate fits in detail ---------------------------------------------------
# allfits <- bestfit %>%
#   filter(model %in% models) %>%
#   group_by(exp,id,condition) %>%
#   mutate(delAIC = AIC - min(AIC)) %>%
#   mutate()
#   filter(condition == "D") %>%
#   filter(exp == "SB2021_e1")
#
#
# winning <- allfits %>%
#   filter(delAIC == 0)
#
#
# GumbelWin <- winning %>% filter(model == "GumbelEVSDT" && delAIC == 0) %>% .$id
# EVSDTWin <- winning %>% filter(model == "GaussianEVSDT" && delAIC == 0) %>% .$id
# UVSDTWin <- winning %>% filter(model == "GaussianUVSDT" && delAIC == 0) %>%
#   relocate(sigo, .after = muo) %>% .$id
#
#
# checkwhatishappening <- allfits %>%
#   mutate(whowins = case_when(id %in% GumbelWin ~ "Gumbel",
#                              id %in% EVSDTWin ~ "EVSDT",
#                              id %in% UVSDTWin ~ "UVSDT",
#                              TRUE ~ "Else")) %>%
#   relocate(sigo, .after = muo) %>%
#   mutate(unstable = ifelse(id %in% exclLL,"Unstable","Stable")) %>%
#   mutate(model = factor(model,levels = models,labels = c("EVSDT","UVSDT","Gumbel"))) %>%
#   mutate(whowins = factor(whowins,levels = c("EVSDT","UVSDT","Gumbel"))) %>%
#   mutate(muo = ifelse(model == "Gumbel",-muo,muo))
#
# NsperWin <- checkwhatishappening %>%
#   group_by(whowins) %>%
#   summarize(nsum = length(delAIC)/3) %>%
#   mutate(model = NA) %>%
#   mutate(model = factor(model,levels = models,labels = c("EVSDT","UVSDT","Gumbel"))) %>%
#   mutate(whowins = factor(whowins,levels = c("EVSDT","UVSDT","Gumbel")))
#
# sigold <- ggplot(checkwhatishappening,aes(x = whowins,y = sigo)) +
#   geom_hline(yintercept = 0)+
#   geom_hline(yintercept = 1,linetype = "dashed")+
#   geom_boxplot(outlier.shape = NA,alpha=0.2,data =checkwhatishappening,
#                aes(x = whowins,y = sigo),fill="lightgrey")+
#   geom_jitter(width=0.3,data =checkwhatishappening,
#               aes(x = whowins,y = sigo,color=model))+
#   scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
#   scale_y_continuous(limits = c(0,12),name = expression(paste(sigma[old], " estimate of best-fitting UVSDT model")))+
#   geom_text(data = NsperWin, aes(x = whowins,y = 12,label=paste0("N = ",nsum)))+
#   scale_x_discrete(name = "Winning Model") +
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12,hjust=0),
#         panel.background = element_rect(fill = "white"),
#         legend.position = "none",
#         panel.border = element_rect(colour = "black", fill = "transparent"),
#         strip.background = element_rect(fill = "white"),
#         strip.placement = "outside",
#         strip.text = element_text(size=12))
#
#
# muold <- ggplot(checkwhatishappening,aes(x = whowins,y = muo)) +
#   geom_hline(yintercept = 0)+
#   geom_boxplot(outlier.shape = NA,alpha=0.2,data =checkwhatishappening,
#                aes(x = whowins,y = muo,fill=model))+
#   scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
#   geom_point(aes(x = whowins,y = abs(muo),color=model),
#              position=position_jitterdodge(dodge.width = 0.8,
#                                            jitter.width = 0.3))+
#   scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
#   geom_text(data = NsperWin, aes(x = whowins,y = 20,label=paste0("N = ",nsum)))+
#   scale_y_continuous(name = expression(paste(mu[old]," estimate")))+
#   scale_x_discrete(name = "Winning Model") +
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12,hjust=0),
#         panel.background = element_rect(fill = "white"),
#         legend.position = "none",
#         panel.border = element_rect(colour = "black", fill = "transparent"),
#         strip.background = element_rect(fill = "white"),
#         strip.placement = "outside",
#         strip.text = element_text(size=12))
#
#
#
#
# delAICplot <- ggplot(checkwhatishappening,aes(x = whowins,y = delAIC)) +
#
#   geom_boxplot(data = checkwhatishappening,
#                aes(x = whowins,y = delAIC,fill=model),alpha=0.2,outlier.shape=NA,color="black")+
#   scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
#   geom_point(data = checkwhatishappening,aes(x = whowins,y = delAIC,color=model),
#              position=position_jitterdodge(dodge.width = 0.8,
#                                            jitter.width = 0.3))+
#   scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
#   geom_text(data = NsperWin, aes(x = whowins,y = 20,label=paste0("N = ",nsum)))+
#   geom_hline(yintercept = 0)+
#   scale_y_continuous(name = expression(paste(Delta, " AIC to winning model")),limits = c(0,20))+
#   scale_x_discrete(name = "Winning Model") +
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12,hjust=0),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(colour = "black", fill = "transparent"),
#         strip.background = element_rect(fill = "white"),
#         strip.placement = "outside",
#         strip.text = element_text(size=12))
#
# exp1condD <- plot_grid(sigold,muold,delAICplot,nrow=1,rel_widths=c(0.5,1,1.2))
#
# plot_grid(exp1condA,
#           exp1condB,
#           exp1condC,
#           exp1condD,nrow=4,labels = c("A","B","C","D"))
#
# ggsave("Exp1_correctzeros.png", units="cm", width=35, height=40, dpi=600)
#
#
# # Parameter estimates ANOVA /exclusions basis of UVSDT ---------------------
#
# # Correlation of muo and sigo estimates in UVSDT
#
# excludeExtreme <- unique(bestfit %>%
#   filter(model == "GaussianUVSDT") %>% ungroup() %>%
#   filter(muo >= mean(muo) + 3 * sd(muo) | sigo >= mean(sigo) + 3 * sd(sigo)) %>%
#   .$id)
#
# bestfit %>%
#   filter(model == "GaussianUVSDT") %>%
#   filter(! id %in% excludeExtreme) %>%
#   group_by(exp) %>%
#   summarize(corrmusig = cor.test(muo,sigo)$estimate[[1]],
#             pcorrmusig = cor.test(muo,sigo)$p.value[[1]])
#
#
# # ANOVA on UVSDT estimates
#
# library(afex)
#
# foranova <- bestfit %>%
#   filter(model == "GaussianUVSDT") %>%
#   filter(! id %in% excludeExtreme) %>%
#   mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
#          variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
#   group_by(exp)
#
# aov_car(sigo ~ Error(id/(strength*variability)),
#                     data = foranova %>% filter(exp == "SB2021_e1") ,factorize = FALSE)
#
# aov_car(muo ~ Error(id/(strength*variability)),
#         data = foranova %>% filter(exp == "SB2021_e1") ,factorize = FALSE)
#
#
#
# foranova <- bestfit %>%
#   filter(model == "GumbelEVSDT") %>%
#   filter(! id %in% excludeExtreme) %>%
#   mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
#          variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
#   group_by(exp)
#
# aov_car(muo ~ Error(id/(strength*variability)),
#         data = foranova %>% filter(exp == "SB2021_e1") ,factorize = FALSE)
#
#
# meanparameters <- bestfit %>%
#   filter(! id %in% excludeExtreme) %>%
#   mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
#          variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
#   group_by(exp,condition,model) %>%
#   mutate(muo = ifelse(model == "GumbelEVSDT",-muo,muo)) %>% #mirror Gumbel
#   summarize(meanmuo = mean(muo),
#             sdmuo = sd(muo),
#             meansigo = mean(sigo),
#             sdsigo = sd(sigo))
#
#
