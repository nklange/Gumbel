

# Model mimicry ----------------------------------------------------------------


GetGsquared <- function(LL,observedData){


  temp <- c(observedData[1:6] * log(observedData[1:6]/sum(observedData[1:6])),
            observedData[7:12] * log(observedData[7:12]/sum(observedData[7:12])))
  temp[observedData == 0] <- 0

  GSq = 2*(LL - -sum(temp) ) # Calculate G^2 (better approximator to Chi-sq dist than Pearson Chi-Sq)
 # f$df     <- length(observedData)-length(names(get_start_par(model)))-1            # df model (L&F,p.63) n_observations - n_freeparas - 1
  #f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq



  return(GSq)
}


fits <- readRDS(paste0("SimulationData/simulate_","genGumbelEVSDT","_data_","LOWN",".rds"))
unique(fits$id)
# LARGEN

genmodels <- c("genGumbelEVSDT","genGaussianUVSDT","genExGaussNormEVSDT")

listids <- unique(data$id)


MimicryFits <- NULL
for (sim in c("LARGEN","LOWN","KEEPCRITS")){
for (gen in genmodels){


fits <- readRDS(paste0("SimulationFits/",gen,"_",sim,"_bestfits.rds")) %>%
  arrange(id,model) %>%
  mutate(simulationtype = sim)

selectedid <- unique(fits %>% .$id)

fits <- fits %>% filter(id %in% selectedid) %>% arrange(id,model)

data <- readRDS(paste0("SimulationData/simulate_",gen,"_data_",sim,".rds")) %>%
  filter(id %in% selectedid) %>% arrange(id)%>%
  mutate(simulationtype = sim)


LL <- fits %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen))) %>%
  group_by(simulationtype,sig,genmodel,id,model) %>%
  group_nest(keep=T,.key="fit")

observedData <- data %>% group_by(simulationtype,sig,genmodel,id) %>% group_nest(keep=T,.key="data") %>%
 dplyr::slice(rep(1:n(),each=4))

getGsquaredvals <- LL %>% bind_cols(observedData %>% select(data)) %>%
  mutate(Gsq = map2(.x = fit,.y = data, .f = ~GetGsquared(.x$objective,.y$Freq))) %>%
  mutate(Dev = map(.x = fit,.f = ~mean(2*.x$objective))) %>%
  mutate(AIC = map(.x = fit,.f = ~mean(.x$AIC))) %>%
  filter(!model == "GumbelUVSDT")

MimicryFits <- MimicryFits %>% bind_rows(getGsquaredvals)

}
}

Bestfits <- MimicryFits %>% select(simulationtype,genmodel,id,model,Gsq,Dev,AIC) %>%
  unnest(cols=c(Gsq,Dev,AIC)) %>%
  group_by(simulationtype,genmodel,id) %>%
  mutate(deltaDev = Dev - min(Dev),
         minDev = min(Dev),
         deltaAIC = AIC - min(AIC)) %>%
  mutate(percentdeltaDev = deltaDev/minDev * 100) %>%
  mutate(winningAIC = ifelse( deltaAIC== 0,1,0),
         winningDev = ifelse( deltaDev== 0,1,0))


Results <- Bestfits %>% group_by(simulationtype,genmodel,model) %>%
  summarize(medianGsq = median(Gsq),
            meanGsq = mean(Gsq),
            sdGsq = sd(Gsq),
            quant25Gsq = quantile(Gsq)[[2]],
            quant75Gsq = quantile(Gsq)[[4]],
            medianDev = median(deltaDev),
            medianpercentDev = median(percentdeltaDev),
            winAIC = sum(winningAIC),
            winDev = sum(winningDev)) %>%
  group_by(simulationtype,genmodel) %>%
  mutate(winpropAIC = winAIC/sum(winAIC),
         winpropDev = winDev/sum(winDev))



test<-Results %>% select(simulationtype,genmodel,model,winpropAIC,medianGsq,medianpercentDev)

