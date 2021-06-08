source("preprocess_data.R")
source("FitUVSDT.R")


# Specific residuals -----------------------------------------------------------
# models <- "ExGaussNormEVSDT"
#
# fullsubj <- NULL
# for (subjid in unique(testphase$id)){
#
#   for(cond in c("A","B","C","D")){
#
#
#     data <- testphase %>%  filter(id == subjid) %>%
#       filter(condition == cond)
#
#     dp <- prep_data(data,freqdat=FALSE)
#
#     ROCs <- NULL
#     partot <- NULL
#
#
#
#     for (modelname in models){
#
#       par <- bestfit %>%
#         filter(model == modelname) %>%
#         filter(id == subjid) %>%
#         filter(condition == cond) %>%
#         ungroup() %>%
#         select(names(get_start_par(modelname)))
#
#       makePredictions <- PredictSDT(data,modelname,as.numeric(par))
#
#       predicted <- tibble(
#         exp = dp$exp,
#         id = dp$id,
#         condition = cond,
#         model = modelname,
#         pred = "predicted",
#         confidence = rep(dp$confidence,2),
#         oldnew = c(rep("New",length(dp$confidence)),rep("Old",length(dp$confidence))),
#         freq = makePredictions$predicted
#       )
#
#       ROCs <- ROCs %>% bind_rows(predicted)
#
#     }
#
#
#
#     observed <-  tibble(
#       exp = dp$exp,
#       id = dp$id,
#       condition = cond,
#       model = "Data",
#       pred = "observed",
#       confidence = rep(dp$confidence,2),
#       oldnew = c(rep("New",length(dp$confidence)),rep("Old",length(dp$confidence))),
#       freq =c(dp$datalist$New,dp$datalist$Old)
#     )
#
#     raw <- bind_rows(ROCs,observed) %>%
#       mutate(confidence = factor(confidence,levels=c(6:1))) %>%
#       group_by(exp,id,condition,model,pred,oldnew) %>%
#       arrange(exp,id,condition,model,pred,oldnew,confidence) %>%
#       mutate(roc = (cumsum(freq)/sum(freq))) %>%
#       mutate(zroc =  qnorm(roc)) %>%
#       mutate(zroc = ifelse(!is.finite(zroc),NA,zroc)) %>%
#       mutate(condid = paste0(id,"_",condition))
#
#     fullsubj <- fullsubj %>% bind_rows(raw)
#
#   }
#
# }

#saveRDS(fullsubj,file = "SB2021_individualpredictions_ExGaussNorm.rds")

models <- c("GaussianEVSDT","GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

fullsubj <- bind_rows(readRDS("SB2021_individualpredictions.rds"),
                      readRDS("SB2021_individualpredictions_ExGaussNorm.rds") %>%
                        filter(!pred=="observed"))

alls <- fullsubj %>%
  mutate(condid = paste0(id,"_",condition)) %>%
  group_by(exp,id,condid,condition,model,oldnew) %>%
  mutate(responseprop = freq/sum(freq))

allresiduals <- NULL
for(subjid in unique(alls$condid)){



  datonly <- alls %>% filter(model == "Data") %>% ungroup() %>%

    filter(condid == subjid)

  for (modelname in models) {

    modelsprop <- alls %>% filter(model == modelname) %>%

      filter(condid == subjid) %>%
      bind_cols(datprop = datonly$responseprop,
                datfreq = datonly$freq) %>%
      mutate(residuals = responseprop - datprop)

    allresiduals <- allresiduals %>% bind_rows(modelsprop)

  }
}

allresiduals <- allresiduals %>%
  mutate(confidence = factor(confidence,levels = c(1:6))) %>%
  mutate(model = factor(model,levels = c(models),labels = c("EVSDT","UVSDT","Gumbel","ExGaussNorm")))



allfits <- bestfit %>%
  filter(model %in% models) %>%
  group_by(exp,id,condition) %>%
  mutate(delAIC = AIC - min(AIC)) %>%
  mutate(condid = paste0(id,"_",condition))


winning <- allfits %>%
  filter(delAIC == 0)


GumbelWin <- winning %>% filter(model == "GumbelEVSDT") %>% .$condid
EVSDTWin <- winning %>% filter(model == "GaussianEVSDT") %>% .$condid
UVSDTWin <- winning %>% filter(model == "GaussianUVSDT") %>%
  relocate(sigo, .after = muo) %>% .$condid
ExGaussNormWin <- winning %>% filter(model == "ExGaussNormEVSDT") %>% .$condid


checkwhatishappening <- allfits %>%   filter(delAIC == 0) %>%
  mutate(whowins = case_when(id %in% GumbelWin ~ "Gumbel",
                             id %in% EVSDTWin ~ "EVSDT",
                             id %in% UVSDTWin ~ "UVSDT",
                             id %in% ExGaussNormWin ~ "ExGaussNorm",
                             TRUE ~ "Else")) %>%
  relocate(sigo, .after = muo) %>%
  mutate(model = factor(model,levels = models,labels = c("EVSDT","UVSDT","Gumbel","ExGaussNorm"))) %>%
  mutate(whowins = factor(whowins,levels = c("EVSDT","UVSDT","Gumbel","ExGaussNorm")))

calcresiduals <- allresiduals %>%
  mutate(whowins = case_when(condid %in% GumbelWin ~ "Gumbel",
                             condid %in% EVSDTWin ~ "EVSDT",
                             condid %in% UVSDTWin ~ "UVSDT",
                             condid %in% ExGaussNormWin ~ "ExGaussNorm",
                             TRUE ~ "Else")) %>%
  mutate(whowins = factor(whowins,levels =c("EVSDT","UVSDT","Gumbel","ExGaussNorm"),
                          labels = c("EVSDT wins","UVSDT wins","Gumbel wins","ExGaussNorm wins"))) %>%
  group_by(exp,condition,whowins,model,oldnew,confidence) %>%
  summarize(mresid = mean(residuals),
            nresid = length(residuals),
            ciresid = qnorm(.975) * sd(residuals)/sqrt(length(residuals)))

Nresids <- calcresiduals %>% ungroup() %>% select(exp,whowins,condition,oldnew,nresid,model) %>%
  filter(model == "EVSDT") %>% filter(oldnew == "New") %>%  distinct()


experiment <- "SB2021_e1"
whowinsmodel <- "UVSDT wins"

plotresiduals <- function(experiment,whowinsmodel){

  if(whowinsmodel == "Gumbel wins"){
    legendpos <- "right"
  } else {
    legendpos <- "none"
  }

  ggplot(data = calcresiduals %>% filter(exp == experiment) %>%
           filter(whowins == whowinsmodel),
         aes(x = confidence,y = mresid,color=model,group=model,fill=model,shape=model)) +

    scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73","#CC79A7"))+
    scale_fill_manual(values = c("#56B4E9","#E69F00", "#009E73","#CC79A7"))+
    scale_shape_manual(values = c(22,21,24,23))+
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim=c(-.05,.05))+
    scale_y_continuous(name="Residuals")+
    geom_line(data = calcresiduals %>% filter(exp == experiment) %>%
                filter(whowins == whowinsmodel),
              aes(x = confidence,y = mresid,color=model,group=model),
              size=1,position = position_dodge(.5))+
    geom_errorbar(data = calcresiduals %>% filter(exp == experiment) %>%
                    filter(whowins == whowinsmodel),
                  aes(x = confidence,ymin = mresid - ciresid,
                      ymax=mresid + ciresid,color=model),
                  position = position_dodge(.5),width=.1)+
    geom_point(data = calcresiduals %>% filter(exp == experiment) %>%
                 filter(whowins == whowinsmodel),
               aes(x = confidence,y = mresid,fill=model,shape=model),
               position = position_dodge(.5),color="black")+
    geom_text(data = Nresids %>% filter(exp == experiment) %>%
                filter(whowins == whowinsmodel),
              aes(x = 3.5, y = .05,label=paste0("N = ",nresid)),color="black")+
    facet_grid(condition~oldnew) +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(size=12),
          #axis.ticks = element_blank(),
          legend.position = legendpos,
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = "transparent",colour="black"),
          strip.placement = "outside",
          #legend.position = "none",
          strip.background = element_rect(fill = "transparent",
                                          colour="transparent"),
          strip.text = element_text(size = 12))
}

Exp1resids <- plot_grid(plotresiduals("SB2021_e1","EVSDT wins"),
                        plotresiduals("SB2021_e1","UVSDT wins"),
                        plotresiduals("SB2021_e1","ExGaussNorm wins"),
                        plotresiduals("SB2021_e1","Gumbel wins"),

                        labels = levels(Nresids$whowins),ncol=4,
                        rel_widths=c(1,1,1,1.1),scale=.95)
Exp2resids <- plot_grid(plotresiduals("SB2021_e2","EVSDT wins"),
                        plotresiduals("SB2021_e2","UVSDT wins"),
                        plotresiduals("SB2021_e1","ExGaussNorm wins"),
                        plotresiduals("SB2021_e2","Gumbel wins"),

                        labels = levels(Nresids$whowins),ncol=4,
                        rel_widths=c(1,1,1,1.1),scale=.95)

plot_grid(Exp1resids,Exp2resids,labels=c("A","B"),nrow=2)

ggsave("Test.png", units="cm", width=35, height=40, dpi=600)

