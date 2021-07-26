  library("RColorBrewer")
  library("ggplot2")
  library("cowplot")
  library("stringr")


  source("FitUVSDT.R")



# check out a reasonable range of parameter estimates

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


datas <- bestfit %>%
  filter(model %in% genmodel) %>% group_by(condition) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  select(c(model,id,condition,muo,betao,c1,dc1,dc2,dc3,dc4)) %>% ungroup() %>%
  pivot_longer(!c(exp,model,id,condition),names_to = "parameter",values_to="value") %>%
  group_by(parameter) %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(parameter,value) %>%
  drop_na() %>%
  select(-c(exp,model,id,condition)) %>%
  select(names(get_start_par(genmodel)))

# for each model: identify mean/sd of paramaters and correlations between them

correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
sigmas <- sds %*% t(sds) * correlations

# Use multivariate to generate parameters on basis of multivariate distribution
# Sample 1000 possible sets of generating parameters per model

# Genpar: mvn from par estimates ------------------------------------------------

set.seed(1)
parameters <- tmvtnorm::rtmvnorm(n=10000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,0,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)
saveRDS(parameters,file="LargeNSimulation/simulate_genExGaussNormEVSDT_parametervalues.rds")

# Simulate data ----------------------------------------------------------------

genmodel <- "GaussianUVSDT"
simulation <- NULL
parametertibble <-NULL
parametersext <- readRDS(paste0("LargeNSimulation/simulate_gen",genmodel,"_parametervalues.rds"))


for(i in c(1:dim(parametersext)[[1]])){

# for(j in c(2:51)){ # for each set of values, make 50 small-N datasets

id <- paste0("gen",genmodel,"_par",i)
#id <- paste0(parametersext[i,]$parid,"_set1")
parid <-paste0("gen",genmodel,"_par",i)# parametersext[i,]$parid#
pars <- as.numeric(parametersext[i,])#parameters[i,]


sim <- PredictSDT(data = NULL, model = genmodel,par = pars,
           itemspertype = c(60,60))

# itemspertype choice as a LargeN (50000) and standard-experiment (100)

simtibble <- sim %>% mutate(id = id,
                            parid = parid,
               genmodel = genmodel)

partibble <- tibble(value = pars,
                    parameter = names(get_start_par(genmodel))) %>%
  mutate(parid = parid,
         genmodel = genmodel)

  # partibble <- pars %>%
  #   set_colnames(names(get_start_par(genmodel))) %>%
  # mutate(parid = parid,
  #        genmodel = genmodel)


simulation <- simulation %>% bind_rows(simtibble)
parametertibble <- parametertibble %>% bind_rows(partibble)



}

saveRDS(simulation,file=paste0("LargeNSimulation/simulate_gen",genmodel,"_data_LN.rds"))
saveRDS(parametertibble,file=paste0("LargeNSimulation/simulate_gen",genmodel,"_parametervalues_LN.rds"))
