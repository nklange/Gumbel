
# Simulation models ------------------------------------------------------------

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

fits <- load_files("Fits2/","fit")

genmodel <- "ExGaussNormEVSDT"
bestfit <- fits %>%
  # mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1)


data <- bestfit %>%
  filter(model %in% genmodel) %>% group_by(condition) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  select(c(model,id,condition,muo,betao,c1,dc1,dc2,dc3,dc4)) %>% ungroup() %>%
  pivot_longer(!c(exp,model,id,condition),names_to = "parameter",values_to="value") %>%
  group_by(condition,parameter) %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(parameter,value) %>%
  drop_na() %>%
  select(-c(exp,model,id)) %>%
  select(c(condition,names(get_start_par(genmodel))))

for(cond in c("C")){
# for each model: identify mean/sd of paramaters and correlations between them

datas <- data %>% filter(condition == cond) %>% ungroup() %>% select(-condition)
correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
sigmas <- sds %*% t(sds) * correlations

# Use multivariate to generate parameters on basis of multivariate distribution
# Sample 1000 possible sets of generating parameters per model

# Genpar: mvn from par estimates ------------------------------------------------

set.seed(1)
parameters <- tmvtnorm::rtmvnorm(n=2000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,0,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)
saveRDS(parameters,
        file=paste0("SimulationData/simulate_genGaussianUVSDT_parametervalues_cond",cond,".rds"))
}

outall <- NULL

for(cond in c("A","B","C","D")){
  # for each model: identify mean/sd of paramaters and correlations between them


  datas <- data %>% filter(condition == cond) %>% ungroup() %>% select(-condition)
  correlations <- cor(datas)
  means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
  sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
  sigmas <- sds %*% t(sds) * correlations


  out <- bind_rows(as_tibble(t(means)) %>% set_colnames(names(get_start_par(genmodel))) %>% mutate(type="M"),
  as_tibble(t(sds)) %>% set_colnames(names(get_start_par(genmodel))) %>% mutate(type="SD")) %>%
    mutate(condition = cond) %>%
    mutate(model = "ExGaussNormEVSDT")

  # Use multivariate to generate parameters on basis of multivariate distribution
  # Sample 1000 possible sets of generating parameters per model

outall <- outall %>% bind_rows(out)

}

outall <- bind_rows(outallUVSDT,outallGumbel,outallEGNorm)
saveRDS(outall,file="mimicry_meansforMVNgeneration.rds")

# Simulate data -----------------------------------------------------------------




genmodel <- "GaussianUVSDT"

for(cond in c("C")){

  simulation <- NULL
  parameters <- readRDS(file=paste0("SimulationData/simulate_genGaussianUVSDT_parametervalues_cond",cond,".rds"))

for(i in c(1:dim(parameters)[[1]])){

    id <- paste0("gen",genmodel,"cond_",cond,"_par",i)

    pars <- as.numeric(parameters[i,])


    sim <- PredictSDT(data = NULL, model = genmodel,par = pars,
                      itemspertype = c(60,60))

    partibble <- as_tibble(t(pars)) %>% set_colnames(names(get_start_par(genmodel))) %>%
      mutate(id = id, genmodel=genmodel,condition=cond) %>% group_by(condition,id,genmodel) %>% group_nest(keep=T)

   simtibble <- sim %>% mutate(id = id,
                               genmodel = genmodel,
                               condition = cond) %>%
     group_by(condition,id,genmodel) %>%
     group_nest(keep=T) %>%
     mutate(genparameters = partibble$data)

   simulation <- simulation %>% bind_rows(simtibble)



}


saveRDS(simulation,file=paste0("SimulationData/simulate_genGaussianUVSDT_data_mimicry_cond",cond,".rds"))

}

# Basic mimicry ----------------------------------------------------------------


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

