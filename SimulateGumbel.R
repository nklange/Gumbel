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


bestfit %>%
  filter(model %in% c("GumbelEVSDT")) %>% group_by(condition) %>%
  filter(muo > mean(muo) - 3 * sd(muo)) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  select(c(model,id,condition,muo,c1,dc1,dc2,dc3,dc4)) %>% ungroup() %>%
  pivot_longer(!c(exp,model,id,condition),names_to = "parameter",values_to="value") %>%
  group_by(parameter) %>%
  summarize(mpar = mean(value),
            sdpar = sd(value))

# make generating parameters sampling from multivariatenormal



# Simulate from Gumbel and fit with UVSDT

PredictSDT(data = NULL, model = "GumbelEVSDT",par = c(-1.12,-1.74,1.20,.7,.45,.53),
           itemspertype = c(100,100) )
