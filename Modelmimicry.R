
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
  filter(model %in% genmodel) %>% #group_by(condition) %>%
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

for(cond in c("A","B","C","D")){
# for each model: identify mean/sd of paramaters and correlations between them

datas <- data %>%
 # filter(condition == cond) %>%
  ungroup() %>% select(-condition)
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
        file=paste0("SimulationData/simulate_genExGaussNormEVSDT_parametervalues.rds"))
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




genmodel <- "ExGaussNormEVSDT"

for(cond in c("C")){

  simulation <- NULL
  parameters <- readRDS(file=paste0("SimulationData/simulate_genExGaussNormEVSDT_parametervalues.rds"))

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


saveRDS(simulation,file=paste0("SimulationData/simulate_genExGaussNormEVSDT_data_mimicry_cond",cond,".rds"))

}
