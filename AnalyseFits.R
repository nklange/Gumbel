library("RColorBrewer")
library("ggplot2")
library("cowplot")

source("preprocess_data.R")

# Model fitting ----------------------------------------------------------------

source("FitUVSDT.R")
models <- c("GaussianEVSDT","GaussianUVSDT","GumbelEVSDT")

for (model in models){

  model <- "ExGaussNormEVSDT"

  for(subjid in unique(testphase$id)){
  # for(experiment in c("SB2021_e1","SB2021_e2")){
    fullsubj <- NULL
    for(cond in c("A","B","C","D")){


      data <- testphase %>% filter(id == subjid) %>% filter(condition==cond)



      fit <- FitSDT(data = data, model = model, rep=20, freqdat = F) %>% mutate(condition = cond)

      fullsubj <- fullsubj %>% bind_rows(fit)
    }


    saveRDS(fullsubj,file = paste0("FitsExGauss/fit_seperateconditions_",model,"_",subjid,".rds"))

  }

}


models <- c("GaussianEVSDT","GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

fits <- load_files("Fits2/","fit")

bestfit <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1)



# Stability fit -------------------

numbermodels <- 5

minDelta <- fits %>%
  filter(model %in% models) %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective)  %>%
  group_by(model,exp,id,condition) %>%
  mutate(mindiffLL = outer(objective,objective, `-`)[2,1],
         bestLL = objective[[1]],
         diffDev = mindiffLL * 2,
         bestDev = bestLL * 2) %>%
  select(model,exp,id,condition,diffDev,bestDev) %>%
  distinct() %>%
  arrange(id) %>%
  ungroup() %>%
  mutate(diffObj = diffDev/bestDev * 100,
         #trials = ntrials$n,
         id = id) %>%
  mutate(model = factor(model,
                        levels =models)) %>%
  mutate(modelnum = as.numeric(model)) %>%

  mutate(modelposition = model) %>%
  mutate(modelposition = factor(modelposition,
                                levels = models)) %>%
  mutate(position = as.numeric(modelposition))



# Deviance in Points figure:

## Text: data points per model not plotted in the graph

PtsAbove10 <- minDelta %>% mutate(above = ifelse(diffDev > 10, 1, 0)) %>%
  group_by(exp,condition,model,modelnum,position) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(exp,condition,model) %>% filter(!notplotted %in% c(0,100) )


## Text: proportion of data sets underneath a threshold per model

thresholdmodel <- minDelta %>%
  group_by(exp,condition,model,modelnum,position) %>%
  # in deviance terms...
  summarize(belowDev05 = length(diffDev[diffDev < 0.5]),
            belowDev1 = length(diffDev[diffDev < 1]),
            belowDev2 = length(diffDev[diffDev < 2]),
            belowDev10 = length(diffDev[diffDev < 10])) %>%
  gather(key = "crit",value="value",-exp,-model,-condition,-modelnum,-position) %>%
  group_by(exp,condition,model,modelnum,position,crit) %>%
  summarize(propbelow = value/64) %>%
  # y axis graph placement
  mutate(yplace = case_when(crit == "belowDev05" ~ 0.25,
                            crit == "belowDev1" ~ 0.75,
                            crit == "belowDev2" ~ 1.5,
                            crit == "belowDev10" ~ 6))


## Text: proportion of data sets underneath a treshold for any of the models

thresholdtotal <- minDelta %>% ungroup() %>%
  mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
         belowDev1 = ifelse(diffDev < 1,1,0),
         belowDev2 = ifelse(diffDev < 2,1,0),
         belowDev10 = ifelse(diffDev < 10,1,0)) %>%
  group_by(exp,id,condition) %>%
  summarize(modelsbelowDev05 = sum(belowDev05),
            modelsbelowDev1 = sum(belowDev1),
            modelsbelowDev2 = sum(belowDev2),
            modelsbelowDev10 = sum(belowDev10)) %>%
  gather(key = "crit",value="value",-id,-exp,-condition) %>%
  group_by(exp,condition,crit) %>%
  summarize(propbelow = length(value[value == 3])/length(value)) %>%
  mutate(yplace = case_when(crit == "modelsbelowDev05" ~ 0.25,
                            crit == "modelsbelowDev1" ~ 0.75,
                            crit == "modelsbelowDev2" ~ 1.5,
                            crit == "modelsbelowDev10" ~ 6))


thresholdtotal2 <- minDelta %>% ungroup() %>%
  mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
         belowDev1 = ifelse(diffDev < 1,1,0),
         belowDev2 = ifelse(diffDev < 2,1,0),
         belowDev10 = ifelse(diffDev < 10,1,0)) %>%
  group_by(exp,id) %>%
  summarize(modelsbelowDev05 = sum(belowDev05),
            modelsbelowDev1 = sum(belowDev1),
            modelsbelowDev2 = sum(belowDev2),
            modelsbelowDev10 = sum(belowDev10)) %>%
  gather(key = "crit",value="value",-id,-exp) %>%
  group_by(exp,crit) %>%
  summarize(propbelow = length(value[value == 3*4])/length(value)) %>%
  mutate(yplace = case_when(crit == "modelsbelowDev05" ~ 0.25,
                            crit == "modelsbelowDev1" ~ 0.75,
                            crit == "modelsbelowDev2" ~ 1.5,
                            crit == "modelsbelowDev10" ~ 6))

# TotalPtsAbove10 <- (1 - thresholdtotal %>%
#                       filter(crit == "modelsbelowDev10") %>%
#                       .$propbelow) * 100

TotalPtsAbove10 <- thresholdtotal %>%
  filter(crit == "modelsbelowDev10") %>%
  mutate(probpercent = (1 - propbelow) * 100)

ggplotColours <- function(n = 10, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

colorin <- ggplotColours(n = 10)

plotcolors <- colorin[unique(thresholdmodel$modelnum)]

datastabpoints <- ggplot(minDelta %>% filter(diffDev <= 10),
                         aes(x = modelnum,y = diffDev )) +

  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 0,ymax=0.5),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 0.5,ymax=1),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 1,ymax=2),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 2,ymax=10),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 10,ymax=10.5),
            fill = "transparent",colour="grey")+

  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 0,ymax=0.5),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 0.5,ymax=1),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 1,ymax=2),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 2,ymax=10),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 10,ymax=10.5),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = 0.75,xmax = max(position) + 1.5,
                ymin = 10.4,ymax=10.5),fill = "white",colour="white")+
  facet_grid(exp~condition)+

  geom_abline(aes(intercept=0.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=2,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=10,slope = 0),linetype="dotted",colour="grey") +

  geom_text(data = thresholdmodel,
            mapping = aes(x = position + 0.2, y = yplace,
                          color = factor(modelnum),
                          label = sub("^0+", "", paste(round(propbelow,2)))),
            size = 3.5) +
  geom_text(data = thresholdtotal,
            mapping = aes(x = numbermodels + 1.25, y = yplace,
                          label = sub("^0+", "", paste(round(propbelow,2)))),
            color = "#A9A9A9",size = 3.5) +

  geom_jitter(data = minDelta %>% filter(diffDev <10), aes(x= position - 0.2,
                                                           color = factor(modelnum)),
              size = 2,alpha=0.6,
              width=0.15) +
  geom_jitter(data = minDelta %>% filter(diffDev <10), aes(x = max(position) + 0.8,
                                                           y = diffDev),
              color = "#A9A9A9",size = 2,alpha=0.6,width=0.15) +
  #scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
  geom_text(data = PtsAbove10,
            mapping = aes(x = position - 0.3, y = 10.5,
                          label = paste0(round(notplotted),"%")),size = 2.5) +
  geom_text(data = TotalPtsAbove10,
            mapping = aes(x = 3 + .75, y = 10.5,
                          label = ifelse(probpercent > 0,paste0(round(probpercent),"%"),NA)),
                          size = 2.5) +


  scale_y_continuous(name = expression(paste(Delta,"Deviance (in points)")),
                     limits = c(0,10.5),
                     breaks = c(c(0,0.5,1),seq(2,10,2))) +
  scale_x_continuous(name = "Model", limits = c(.5,3 + 1.5),
                     breaks = c(c(1:3),3 + 1),
                     labels = c("EVSDT","UVSDT","Gumbel","Total")
  ) +

  theme(axis.text.x = element_text(size=11),#,angle=90,vjust=0.5
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent",colour="transparent"))


# Analyse fits -----------------------------------------------------------------


exclLL <- unique(fits %>%
  filter(model %in% models) %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective)  %>%
  group_by(model,exp,id,condition) %>%
  dplyr::slice(1:2) %>%
  mutate(mindiffLL =diff(objective)) %>%
  select(exp,id,condition,model,objective,mindiffLL) %>%
  dplyr::slice(1) %>%
  filter(mindiffLL > 0.5) %>%
  .$id)


propbestfit <- function(models_charvector,bestfit){

  AICmeans <- bestfit %>%
    filter(model %in% models_charvector) %>%
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
      filter(model %in% models_charvector) %>%
      mutate(deltaAIC = AIC - AICmeans1 %>% dplyr::slice(1) %>% .$means) %>%
      group_by(exp,condition,model) %>%
      mutate(model = factor(model,levels = models)) %>%

      summarize(mdeltaAIC = mean(deltaAIC)) %>%
      mutate(condition = factor(condition,
                                levels = c("A","B","C","D"),
                                labels = c("High strength\nHigh variability",
                                           "High strength\nLow variability",
                                           "Low strength\nHigh variability",
                                           "Low strength\nLow variability")))

    full <- full %>% bind_rows(AICtot)

    }
  }



  indstrict_PW <- bestfit %>%
    filter(model %in% models_charvector) %>%
    group_by(exp,id,condition) %>%

    mutate(delAIC = AIC - min(AIC)) %>%
    #filter(id == subjid) %>% .$delAIC
    filter(delAIC == 0) %>% # only 1 model can win
    group_by(exp,condition,model) %>%
    summarize(numBest = length(delAIC)) %>%
    mutate(model = factor(model,levels = models)) %>%


    mutate(condition = factor(condition,
                              levels = c("A","B","C","D"),
                              labels = c("High strength\nHigh variability",
                                         "High strength\nLow variability",
                                         "Low strength\nHigh variability",
                                         "Low strength\nLow variability"))) %>%
    mutate(exp = factor(exp)) %>%
    complete(model,nesting(exp,condition), fill = list(0)) %>%
    group_by(exp,condition,model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    group_by(exp,condition) %>%
    mutate(wu =wu/sum(wu)) %>%
    arrange(exp,condition,model)


  ylim.prim <- c(0,11)
  ylim.sec <- c(0,1)

  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- b*(ylim.prim[1] - ylim.sec[1])



  plot <- ggplot(full,
                           aes(y = mdeltaAIC,x = model )) +
    #geom_rect(aes(xmin=0, xmax=24, ymin=0, ymax=1), fill = "grey") +
    facet_grid(exp~condition)+
    geom_bar(stat = "identity",fill="lightgrey") +
    coord_flip(ylim = c(0,11)) +
    geom_point(aes(y = a +indstrict_PW$wu*b),
               color = "black",size=2,shape=21,fill="black") +


    #scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
    # scale_fill_gradientn(name = expression(paste(Delta)),
    #                      colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
    #                      limits=c(0,50))+
   # geom_text(aes(label=round(labeldelta), y = 52), hjust = 0.5,size=3)+
    scale_y_continuous(name = expression(paste(Delta," AIC")),
                       breaks = seq(0,10,2),
                       labels = c(seq(0,8,2),"..."),
                       sec.axis = sec_axis(~ (. - a)/b,
                                           name = "Proportion datasets best fit",
                                           labels = c("0",".25",".50",".75","1"))) +

    scale_x_discrete(name = "Model") +

    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12,hjust=0),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill = "transparent"),
          strip.background = element_rect(fill = "white"),
          strip.placement = "outside",
          strip.text = element_text(size=12))


return(plot)
}


models_charvector <- c("GaussianEVSDT","GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

models <- models_charvector
DeltaAICexclLL <- propbestfit(models_charvector,
                        bestfit %>% filter(!id  %in% exclLL))

DeltaAIC <- propbestfit(models_charvector,
                              bestfit)


differences <- bestfit %>%
  filter(model %in% models_charvector) %>%
  mutate(model = factor(model,levels = models)) %>%
  group_by(exp,id,condition) %>%
  mutate(delAIC = AIC - min(AIC)) %>%
  mutate(colordelAIC = case_when(delAIC == 0 ~ 0,
                                 delAIC > 0 & delAIC <= 1 ~ 1,
                                 delAIC > 1 & delAIC <= 2 ~ 2,
                                 delAIC > 2 & delAIC <= 5 ~ 3,
                                 delAIC > 5 & delAIC <= 10 ~ 4,
                                 delAIC > 10 & delAIC <= 25 ~ 5,
                                 TRUE ~ 6)) %>%
  mutate(condition = factor(condition,
                            levels = c("A","B","C","D"),
                            labels = c("High strength\nHigh variability",
                                       "High strength\nLow variability",
                                       "Low strength\nHigh variability",
                                       "Low strength\nLow variability")))

ggplot(differences %>% filter(exp == "SB2021_e2"),aes(y = delAIC,x=id)) +
  geom_bar(stat = "identity")+
  facet_grid(model~condition,scales="free_x")


IndividualDelta <- ggplot(data = differences,
               aes(y = model,x = id)) +
  geom_tile(data = differences, aes(fill = factor(colordelAIC))) +


  #geom_text(aes(label = round(delAIC)),size = 3.5) +
  scale_fill_manual(values = c("#000000","#bd0026","#f03b20","#fd8d3c","#feb24c","#fed976","#ffffb2"),
                    labels = c("= 0","> 0 & <= 1","> 1 & <= 2","> 2 & <= 5","> 5 & <= 10","> 10 & <= 25","> 25"),
                    name = "Delta AIC")+
  scale_x_discrete(position = "top",name = "id") +


  facet_wrap(exp~condition,scales="free_x",nrow=2) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour="black"),
        strip.placement = "outside",
        #legend.position = "none",
        strip.background = element_rect(fill = "transparent",
                                        colour="transparent"),
        strip.text = element_text(size = 12))



plot_grid(datastabpoints,DeltaAIC,IndividualDelta,nrow=3,rel_heights = c(0.65,0.5,0.5),
          labels = c("A","B","C"))

ggsave("Test2.png", units="cm", width=30, height=30, dpi=600)

# Investigate fits in detail ---------------------------------------------------
allfits <- bestfit %>%
  filter(model %in% models) %>%
  group_by(exp,id,condition) %>%
  mutate(delAIC = AIC - min(AIC)) %>%
  filter(condition == "D") %>%
  filter(exp == "SB2021_e1")


winning <- allfits %>%
  filter(delAIC == 0)


GumbelWin <- winning %>% filter(model == "GumbelEVSDT" && delAIC == 0) %>% .$id
EVSDTWin <- winning %>% filter(model == "GaussianEVSDT" && delAIC == 0) %>% .$id
UVSDTWin <- winning %>% filter(model == "GaussianUVSDT" && delAIC == 0) %>%
  relocate(sigo, .after = muo) %>% .$id


checkwhatishappening <- allfits %>%
  mutate(whowins = case_when(id %in% GumbelWin ~ "Gumbel",
                             id %in% EVSDTWin ~ "EVSDT",
                             id %in% UVSDTWin ~ "UVSDT",
                             TRUE ~ "Else")) %>%
  relocate(sigo, .after = muo) %>%
  mutate(unstable = ifelse(id %in% exclLL,"Unstable","Stable")) %>%
  mutate(model = factor(model,levels = models,labels = c("EVSDT","UVSDT","Gumbel"))) %>%
  mutate(whowins = factor(whowins,levels = c("EVSDT","UVSDT","Gumbel"))) %>%
  mutate(muo = ifelse(model == "Gumbel",-muo,muo))

NsperWin <- checkwhatishappening %>%
  group_by(whowins) %>%
  summarize(nsum = length(delAIC)/3) %>%
  mutate(model = NA) %>%
  mutate(model = factor(model,levels = models,labels = c("EVSDT","UVSDT","Gumbel"))) %>%
  mutate(whowins = factor(whowins,levels = c("EVSDT","UVSDT","Gumbel")))

sigold <- ggplot(checkwhatishappening,aes(x = whowins,y = sigo)) +
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 1,linetype = "dashed")+
  geom_boxplot(outlier.shape = NA,alpha=0.2,data =checkwhatishappening,
               aes(x = whowins,y = sigo),fill="lightgrey")+
  geom_jitter(width=0.3,data =checkwhatishappening,
              aes(x = whowins,y = sigo,color=model))+
  scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
  scale_y_continuous(limits = c(0,12),name = expression(paste(sigma[old], " estimate of best-fitting UVSDT model")))+
  geom_text(data = NsperWin, aes(x = whowins,y = 12,label=paste0("N = ",nsum)))+
  scale_x_discrete(name = "Winning Model") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12,hjust=0),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


muold <- ggplot(checkwhatishappening,aes(x = whowins,y = muo)) +
  geom_hline(yintercept = 0)+
  geom_boxplot(outlier.shape = NA,alpha=0.2,data =checkwhatishappening,
               aes(x = whowins,y = muo,fill=model))+
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
  geom_point(aes(x = whowins,y = abs(muo),color=model),
             position=position_jitterdodge(dodge.width = 0.8,
                                           jitter.width = 0.3))+
  scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
  geom_text(data = NsperWin, aes(x = whowins,y = 20,label=paste0("N = ",nsum)))+
  scale_y_continuous(name = expression(paste(mu[old]," estimate")))+
  scale_x_discrete(name = "Winning Model") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12,hjust=0),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))




delAICplot <- ggplot(checkwhatishappening,aes(x = whowins,y = delAIC)) +

  geom_boxplot(data = checkwhatishappening,
               aes(x = whowins,y = delAIC,fill=model),alpha=0.2,outlier.shape=NA,color="black")+
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
  geom_point(data = checkwhatishappening,aes(x = whowins,y = delAIC,color=model),
             position=position_jitterdodge(dodge.width = 0.8,
                                           jitter.width = 0.3))+
  scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73"))+
  geom_text(data = NsperWin, aes(x = whowins,y = 20,label=paste0("N = ",nsum)))+
  geom_hline(yintercept = 0)+
  scale_y_continuous(name = expression(paste(Delta, " AIC to winning model")),limits = c(0,20))+
  scale_x_discrete(name = "Winning Model") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12,hjust=0),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))

exp1condD <- plot_grid(sigold,muold,delAICplot,nrow=1,rel_widths=c(0.5,1,1.2))

plot_grid(exp1condA,
          exp1condB,
          exp1condC,
          exp1condD,nrow=4,labels = c("A","B","C","D"))

ggsave("Exp1_correctzeros.png", units="cm", width=35, height=40, dpi=600)


# Parameter estimates ANOVA /exclusions basis of UVSDT ---------------------

# Correlation of muo and sigo estimates in UVSDT

excludeExtreme <- unique(bestfit %>%
  filter(model == "GaussianUVSDT") %>% ungroup() %>%
  filter(muo >= mean(muo) + 3 * sd(muo) | sigo >= mean(sigo) + 3 * sd(sigo)) %>%
  .$id)

bestfit %>%
  filter(model == "GaussianUVSDT") %>%
  filter(! id %in% excludeExtreme) %>%
  group_by(exp) %>%
  summarize(corrmusig = cor.test(muo,sigo)$estimate[[1]],
            pcorrmusig = cor.test(muo,sigo)$p.value[[1]])


# ANOVA on UVSDT estimates

library(afex)

foranova <- bestfit %>%
  filter(model == "GaussianUVSDT") %>%
  filter(! id %in% excludeExtreme) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp)

aov_car(sigo ~ Error(id/(strength*variability)),
                    data = foranova %>% filter(exp == "SB2021_e1") ,factorize = FALSE)

aov_car(muo ~ Error(id/(strength*variability)),
        data = foranova %>% filter(exp == "SB2021_e1") ,factorize = FALSE)



foranova <- bestfit %>%
  filter(model == "GumbelEVSDT") %>%
  filter(! id %in% excludeExtreme) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp)

aov_car(muo ~ Error(id/(strength*variability)),
        data = foranova %>% filter(exp == "SB2021_e1") ,factorize = FALSE)


meanparameters <- bestfit %>%
  filter(! id %in% excludeExtreme) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  mutate(muo = ifelse(model == "GumbelEVSDT",-muo,muo)) %>% #mirror Gumbel
  summarize(meanmuo = mean(muo),
            sdmuo = sd(muo),
            meansigo = mean(sigo),
            sdsigo = sd(sigo))
