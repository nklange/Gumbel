# in UVSDT model, correlation of mu and sigma
# generate data from other models


library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(readr)
library(tibble)
library(ggplot2)
library(stringr)
library(cowplot)
source("FitUVSDT.R")


# UVSDT from Gumbel-generated datasets ------------------------------------------------






GetGsquared <- function(LL,observedData){


  temp <- c(observedData[1:6] * log(observedData[1:6]/sum(observedData[1:6])),
            observedData[7:12] * log(observedData[7:12]/sum(observedData[7:12])))
  temp[observedData == 0] <- 0

  GSq = 2*(LL - -sum(temp) ) # Calculate G^2 (better approximator to Chi-sq dist than Pearson Chi-Sq)
  # f$df     <- length(observedData)-length(names(get_start_par(model)))-1            # df model (L&F,p.63) n_observations - n_freeparas - 1
  #f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq



  return(GSq)
}

# MimicryFits <- NULL
# for (gen in c("genGumbelEVSDT","genGaussianUVSDT","genExGaussNormEVSDT")){
#
#
#   fits <- readRDS(paste0("SimulationFits/",gen,"_mimicry_bestfits.rds")) %>%
#     arrange(id,model) %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen)))
#
#   data <- readRDS(paste0("SimulationData/simulate_",gen,"_data_mimicry.rds")) %>%
#     arrange(id)
#
#   LL <- fits %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen))) %>%
#     group_by(genmodel,id,model) %>%
#     group_nest(keep=T,.key="fit")
#
#   observedData <- data %>% group_by(genmodel,id) %>%
#     dplyr::slice(rep(1:n(),each=3))
#
#   getGsquaredvals <- LL %>% bind_cols(observedData %>% ungroup() %>%  dplyr::select(data)) %>%
#     bind_cols(observedData %>% ungroup() %>%  dplyr::select(genparameters)) %>%
#     mutate(Gsq = map2(.x = fit,.y = data, .f = ~GetGsquared(.x$objective,.y$Freq))) %>%
#     mutate(objective = map(.x = fit,.f = ~mean(.x$objective))) %>%
#     mutate(Dev = map(.x = fit,.f = ~mean(2*.x$objective))) %>%
#     mutate(AIC = map(.x = fit,.f = ~mean(.x$AIC)))
#
#
#   MimicryFits <- MimicryFits %>% bind_rows(getGsquaredvals)
#
# }


MimicryFits <- NULL
for (gen in c("genGumbelEVSDT","genGaussianUVSDT","genExGaussNormEVSDT")){

  fits <- readRDS("LargeNSimulation/LargeN_mimicry_bestfits.rds") %>%
    filter(genmodel == substr(gen,start=4,stop=nchar(gen))) %>% arrange(id,model)

  # fits <- readRDS(paste0("SimulationFits/",gen,"_mimicry_bestfits.rds")) %>%
  # arrange(id,model) %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen)))

  data <- readRDS(paste0("LargeNSimulation/simulate_",gen,"_data_LN.rds")) %>%
    arrange(id)

  # data <- readRDS(paste0("SimulationData/simulate_",gen,"_data_mimicry.rds")) %>%
  #   arrange(id)


  LL <- fits %>% #mutate(genmodel = substr(gen,start=4,stop=nchar(gen))) %>%
    group_by(genmodel,id,model) %>%
    group_nest(keep=T,.key="fit")

  # LL <- fits %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen))) %>%
  #   group_by(genmodel,id,model) %>%
  #   group_nest(keep=T,.key="fit")

  observedData <- data %>% group_by(genmodel,id)  %>% group_nest(keep=T)%>%
    dplyr::slice(rep(1:n(),each=3))

  parameters <- readRDS(paste0("LargeNSimulation/simulate_",gen,"_parametervalues_LN.rds")) %>%
    pivot_wider(values_from="value",names_from="parameter") %>%  group_by(genmodel,parid)  %>%
    group_nest(keep=T)%>%dplyr::slice(rep(1:n(),each=3)) %>% rename("genparameters" = data)


  # observedData <- data %>% group_by(genmodel,id) %>%
  #   dplyr::slice(rep(1:n(),each=3))

  getGsquaredvals <- LL %>% bind_cols(observedData %>% ungroup() %>%  dplyr::select(data)) %>%
    bind_cols(parameters %>% ungroup() %>% dplyr::select(genparameters)) %>%
    mutate(Gsq = map2(.x = fit,.y = data, .f = ~GetGsquared(.x$objective,.y$Freq))) %>%
    mutate(objective = map(.x = fit,.f = ~mean(.x$objective))) %>%
    mutate(Dev = map(.x = fit,.f = ~mean(2*.x$objective))) %>%
    mutate(AIC = map(.x = fit,.f = ~mean(.x$AIC)))


  MimicryFits <- MimicryFits %>% bind_rows(getGsquaredvals)

}

testcorr <- MimicryFits %>% .$fit %>% bind_rows()


load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

fits <- load_files("Fits2/","fit")

bestfitDat <- fits %>%
  #mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1) %>% filter(model=="GaussianUVSDT") %>%
  mutate(genmodel = "SB2021")
lim <- 5
cortext <- tibble(genmodel=c("SB2021","ExGaussNormEVSDT","GumbelEVSDT","GaussianUVSDT"),
                  corrlim = c(cor(bestfitDat %>% filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$muo,
                               bestfitDat %>% filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$sigo),
                    cor(testcorr %>% filter(genmodel == "ExGaussNormEVSDT") %>%
                          filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$muo,
                               testcorr %>% filter(genmodel == "ExGaussNormEVSDT") %>%
                          filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$sigo),
                           cor(testcorr %>% filter(genmodel == "GumbelEVSDT") %>%
                                 filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$muo,
                               testcorr %>% filter(genmodel == "GumbelEVSDT") %>%
                                 filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$sigo),
                           cor(testcorr %>% filter(genmodel == "GaussianUVSDT") %>%
                                 filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$muo,
                               testcorr %>% filter(genmodel == "GaussianUVSDT") %>%
                                 filter(muo < lim & sigo < lim) %>%
                                 filter(model == "GaussianUVSDT") %>% .$sigo)),
                  corr = c(cor(bestfitDat %>%
                                    filter(model == "GaussianUVSDT") %>% .$muo,
                                  bestfitDat %>%
                                    filter(model == "GaussianUVSDT") %>% .$sigo),
                              cor(testcorr %>% filter(genmodel == "ExGaussNormEVSDT") %>%

                                    filter(model == "GaussianUVSDT") %>% .$muo,
                                  testcorr %>% filter(genmodel == "ExGaussNormEVSDT") %>%

                                    filter(model == "GaussianUVSDT") %>% .$sigo),
                              cor(testcorr %>% filter(genmodel == "GumbelEVSDT") %>%

                                    filter(model == "GaussianUVSDT") %>% .$muo,
                                  testcorr %>% filter(genmodel == "GumbelEVSDT") %>%

                                    filter(model == "GaussianUVSDT") %>% .$sigo),
                              cor(testcorr %>% filter(genmodel == "GaussianUVSDT") %>%

                                    filter(model == "GaussianUVSDT") %>% .$muo,
                                  testcorr %>% filter(genmodel == "GaussianUVSDT") %>%

                                    filter(model == "GaussianUVSDT") %>% .$sigo))) %>%
  mutate(genmodel = factor(genmodel, levels=c("SB2021","GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                           labels=c("SB2021","Gen: UVSDT","Gen: Gumbel","Gen: ExGauss"))) %>%
  mutate(designatealpha = NA)

testcorr <- testcorr %>% mutate(genmodel = factor(genmodel, levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                                      labels=c("Gen: UVSDT","Gen: Gumbel","Gen: ExGauss")))





correlationUVSDT <- bind_rows(testcorr,bestfitDat) %>%
  mutate(genmodel = factor(genmodel, levels=c("SB2021","Gen: UVSDT","Gen: Gumbel","Gen: ExGauss"))) %>%
  mutate(designatealpha = ifelse(genmodel == "SB2021","high","low"))

uvsdtmuosigo <- ggplot(correlationUVSDT %>% filter(model == "GaussianUVSDT"),
                       aes(x=muo,y=sigo,color=genmodel,alpha=designatealpha,size=designatealpha))+
  # annotate("text", x = 0, y = 4, aes(label = lm_eqn_poly1(testcorr %>% filter(model == "GaussianUVSDT"))),
  #          parse = TRUE, color="black",hjust = 0)+
  geom_point() +
  facet_grid(.~genmodel)+
  scale_alpha_discrete(range=c(.1,.01))+
  scale_size_discrete(range=c(2,1))+
  geom_text(data = cortext,aes(x = 5,y=4.5,label=paste0("graph: r = ", stringr::str_remove(round(corrlim,2), "^0+")),
                               color=genmodel),hjust = 1)+
  geom_text(data = cortext,aes(x = 5,y=4,label=paste0("total: r = ", stringr::str_remove(round(corr,2), "^0+")),
                               color=genmodel),hjust = 1)+
  #geom_smooth()+
  scale_color_manual(values = c("#616161","#E69F00","#009E73","#0072B2"))+
  scale_x_continuous(name=expression(paste(mu[s], "(estimated by UVSDT)")),limits=c(-0.5,5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste(sigma[s]," (estimated by UVSDT)")),limits=c(-0.5,5))+

  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


parametervals <- MimicryFits %>% .$genparameters %>% bind_rows() %>% ungroup() %>% distinct() %>%
  filter(genmodel == "GumbelEVSDT")
bestfit <- testcorr %>% filter(genmodel=="Gen: Gumbel")

gumbelest <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

# gumbeluvsdtest <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
#   arrange(id) %>% .$muo
#
# gumbelfit <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
#   arrange(id) %>% .$message

uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

uvsdtsigma <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,sigo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$sigo

# gumbeluvsdtsigma <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,betao) %>% mutate(id = str_remove(id,"_set1")) %>%
#   arrange(id) %>% mutate(sigo = pi/sqrt(6)*betao) %>% .$sigo



compareboths <- parametervals %>% select(parid,muo) %>% arrange(parid) %>% mutate(muo = -muo)

modelorder <- c("Gumbel","UVSDT")
recoveringmu <- bind_rows(compareboths,
                          compareboths) %>%
  mutate(muo_estimate = c(-gumbelest,uvsdtest),
         #   message = c(gumbelfit,uvsdtfit),
         model = c(
                   rep("Gumbel",10000),
                   rep("UVSDT",10000))) %>%
  mutate(model = factor(model,levels=modelorder))

recoveringsigma <- bind_rows(compareboths,
                             compareboths) %>% mutate(sigo_estimate = c(
                                                                                     rep(pi/sqrt(6),10000),
                                                                                     uvsdtsigma),
                                                                   #   message = c(gumbelfit,uvsdtfit),
                                                                   model = c(
                                                                             rep("Gumbel",10000),
                                                                             rep("UVSDT",10000))
                             ) %>%
  mutate(model = factor(model,levels=modelorder))



muoestimateall <- ggplot(recoveringmu %>% filter(model != "Gumbel (free sigma[o])") ,
                         aes(y=muo_estimate,
                             x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.01) +
  geom_abline(slope=1,intercept=0)+
  scale_color_manual(name = "Fitted model",values = c("#009E73","#E69F00"),
                     labels=c(#expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("Gumbel")),
                              expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[s] - mu[n],"| (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated |", mu[s] - mu[n],"|")),limits=c(-0.5,5),breaks = c(0:5))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +

  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "transparent",color="transparent"),
    legend.box.background = element_rect(fill = "transparent",color="transparent"),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

sigoestimateall <- ggplot(recoveringsigma %>% filter(model != "Gumbel (free sigma[o])"),
                          aes(y=sigo_estimate,
                              x=muo))+
  geom_line(data =  recoveringsigma %>% filter(model == "Gumbel"),
            aes(y=sigo_estimate, x=muo,linetype=model),color="#009E73",size=1)+
  geom_point(data = recoveringsigma %>% filter(model == "UVSDT"), aes(y=sigo_estimate,
                                                                       x=muo,color=model),alpha=0.01) +
  scale_color_manual(name = "Fitted model",values = c("#E69F00"),
                     labels=c(#expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("UVSDT"))))+

  scale_linetype_manual(name = "Generating model",values="dashed")+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)),
         linetype = guide_legend(override.aes = list(size = 1,
                                                     alpha = 1) ) ) +
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[s] - mu[n],"| (Gumbel, generating value)")),limits=c(0,4),
                     breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", sigma[s])),limits=c(0,4),breaks = c(0:4))+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.background = element_rect(fill = "transparent",color="transparent"),
    legend.box.background = element_rect(fill = "transparent",color="transparent"),
    legend.position = c(0.5,0.8),
    #legend.text.align = 0,
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))



LargeN_genGumbel <- plot_grid(muoestimateall,sigoestimateall,rel_widths = c(1,1),nrow=1)


# Gen ExGaussNorm --------------------------------------------------------------

parametervals <- MimicryFits %>% .$genparameters %>% bind_rows() %>% ungroup() %>% distinct() %>%
  filter(genmodel == "ExGaussNormEVSDT")
bestfit <- testcorr %>% filter(genmodel=="Gen: ExGauss")

exgest <- bestfit %>% filter(model == "ExGaussNormEVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo



exgaussbeta <- bestfit %>% filter(model == "ExGaussNormEVSDT") %>% select(id,betao) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$betao

lambda <- 1/exgaussbeta
exgausssigma <- sqrt(1 + (1/(lambda^2)))



uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

uvsdtsigma <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,sigo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$sigo



compareboths <- parametervals %>% select(parid,muo) %>% arrange(parid)

modelorder <- c("ExGauss","UVSDT")

recoveringmu <- bind_rows(compareboths,compareboths) %>%
  mutate(muo_estimate = c(exgest,uvsdtest),
         #   message = c(gumbelfit,uvsdtfit),
         model = c(
                   rep("ExGauss",10000),
                   rep("UVSDT",10000))) %>%
  mutate(model = factor(model,levels=modelorder))

modelorder <- c("ExGauss","UVSDT")

recoveringsigma <- bind_rows(compareboths,compareboths) %>% mutate(sigo_estimate = c(
                                                                                     exgausssigma,
                                                                        uvsdtsigma),
                                                      #   message = c(gumbelfit,uvsdtfit),
                                                      model = c(
                                                                rep("ExGauss",10000),
                                                                rep("UVSDT",10000))
                             ) %>%
  mutate(model = factor(model,levels=modelorder))



muoestimateall <- ggplot(recoveringmu ,
                         aes(y=muo_estimate,
                             x=muo,color=model))+
  geom_abline(slope=1,intercept=0)+

  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.01) +
  scale_color_manual(name = "Fitted model",values = c("#0072B2","#E69F00"),
                     labels=c(
                              expression(paste("ExGauss")),
                              expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[s] - mu[n],"| (ExGauss, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated |", mu[s] - mu[n],"|")),limits=c(-0.5,5),breaks = c(0:5))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +

  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "transparent",color="transparent"),
    legend.box.background = element_rect(fill = "transparent",color="transparent"),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

sigoestimateall<-ggplot(recoveringsigma,
       aes(y=sigo_estimate,
           x=muo))+
  geom_point(data = recoveringsigma %>% filter(model != "sigma[G] (ExGauss)"), aes(y=sigo_estimate,
                                                                       x=muo,color=model),alpha=0.01) +
  scale_color_manual(name = "Fitted model",values = c("#0072B2","#E69F00"),
                     labels=c(expression(paste("ExGauss")),
                              expression(paste("UVSDT"))))+
  # geom_line(data =  recoveringsigma %>% filter(model == "sigma[G] (ExGauss)"),
  #           aes(y=sigo_estimate, x=muo,linetype=model),color="#CC79A7",size=1)+
  scale_linetype_manual(name = "Generating model",values="dashed")+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)),
         linetype = guide_legend(override.aes = list(size = 1,
                                                     alpha = 1) ) ) +
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[s] - mu[n],"| (ExGauss, generating value)")),limits=c(0,5.5),
                     breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", sigma[s])),limits=c(0,8),breaks = c(0:8))+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),
    legend.background = element_rect(fill = "transparent",color="transparent"),
    legend.box.background = element_rect(fill = "transparent",color="transparent"),
    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.8),
    legend.text.align = 0,
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


LargeN_genExG <- plot_grid(muoestimateall,sigoestimateall,rel_widths = c(1,1),nrow=1)

correlationsplot1 <- plot_grid(LargeN_genGumbel,LargeN_genExG,ncol=2,labels=c("A","B"))
correlationsplot2 <- plot_grid(correlationsplot1,uvsdtmuosigo,nrow=2,labels=c("","C"))

ggsave(paste0("Figures/UVSDTcorrelations_LN.png"), units="cm", width=35, height=15, dpi=600)

