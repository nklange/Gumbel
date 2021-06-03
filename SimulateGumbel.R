library("RColorBrewer")
library("ggplot2")
library("cowplot")

source("preprocess_data.R")
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


bestfit <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1)


datas <- bestfit %>%
  filter(model %in% c("GumbelEVSDT")) %>% group_by(condition) %>%
  filter(muo > mean(muo) - 3 * sd(muo)) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  select(c(model,id,condition,muo,c1,dc1,dc2,dc3,dc4)) %>% ungroup() %>%
  pivot_longer(!c(exp,model,id,condition),names_to = "parameter",values_to="value") %>%
  group_by(parameter) %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(parameter,value) %>%
  drop_na() %>%
  select(-c(exp,model,id,condition)) %>%
  select(names(get_start_par("GumbelEVSDT")))

# for each model: identify mean/sd of paramaters and correlations between them

correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
sigmas <- sds %*% t(sds) * correlations

# Use multivariate to generate parameters on basis of multivariate distribution
# Sample 1000 possible sets of generating parameters per model

set.seed(1)
parameters <- tmvtnorm::rtmvnorm(n=1000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)


# Simulate from Gumbel and fit with UVSDT

genmodel <- "GumbelEVSDT"

simulation <- NULL
parametertibble <- NULL

for(i in c(1:1000)){

  for(j in c(2:101)){

id <- paste0("gen",genmodel,"_par",i,"_set",j)
parid <- paste0("gen",genmodel,"_par",i)
pars <- parameters[i,]


sim <- PredictSDT(data = NULL, model = "GumbelEVSDT",par = pars,
           itemspertype = c(200,200))

# itemspertype choice as high enough where parameter recovery is pretty spot-on.

simtibble <- sim %>% mutate(id = id,
                            parid = parid,
               genmodel = genmodel)

partibble <- tibble(value = pars,
                    parameter = names(get_start_par(genmodel))) %>%
  mutate(parid = parid,
         genmodel = genmodel)

simulation <- simulation %>% bind_rows(simtibble)
parametertibble <- parametertibble %>% bind_rows(partibble)

}

}

saveRDS(simulation,file="simulate_genGumbelEVSDT_data_trials400.rds")
#saveRDS(parametertibble,file="simulate_genGumbelEVSDT_parametervalues.rds")
library("doParallel")
library("foreach")

doParallel::registerDoParallel(cores=16)
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
model <- "GumbelEVSDT"
model <- "GaussianUVSDT"

res <- foreach(subjid = unique(simulation %>% .$id)
               ,.combine = 'rbind'
               ,.options.multicore=mcoptions) %dopar% {

                 fitted <- FitSDT(simulation %>% filter(id=="genGumbelEVSDT_par1_set1"),
                                  model,rep=20,freqdat=T)

        saveRDS(fitted, file = paste0("FitsSimulation/", subjid, "_", model, ".rds"))
               }






