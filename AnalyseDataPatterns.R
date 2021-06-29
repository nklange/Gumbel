source("preprocess_data.R")

# Load fits

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
  dplyr::slice(1) %>%
  filter(model %in%  c("GaussianUVSDT","GaussianEVSDT","GumbelEVSDT")) %>%
  group_by(exp,id,condition) %>%
  mutate(delAIC = AIC - min(AIC))


# Visualize data patterns

frequencies <- NULL
for(subj in unique(testphase$id)){
  for(cond in unique(testphase$condition)){


    dp <- prep_data(testphase %>% filter(id == subj) %>% filter(condition == cond),freqdat=F)

    ind <- tibble(rating = as.character(rep(dp$confidence,2)),
                  oldnew = rep(c("New","Old"),each= 6),
                  freq = c(dp$datalist$New,dp$datalist$Old),
                  id = dp$id,
                  exp = dp$exp,
                  condition = cond) %>%
      mutate(freq = ifelse(freq < 1, 0, freq))

    frequencies <- frequencies %>% bind_rows(ind)

  }
}

incompletes <- frequencies %>% mutate(freq0 = ifelse(freq == 0,0,1))

# find misfits

colorthings <- bestfit %>% mutate(freq0 = case_when(
  model == "GumbelEVSDT" & delAIC > 2 ~ 3,
  model == "GaussianEVSDT" & delAIC > 2 ~ 4,
  model == "GaussianUVSDT" & delAIC > 2 ~ 5,
  TRUE ~ 2
)) %>% mutate(rating = case_when(model == "GumbelEVSDT" ~ "1",
                                 model == "GaussianEVSDT" ~ "2",
                                 model == "GaussianUVSDT" ~ "3")) %>%

  select(exp,id,condition,rating,freq0) %>%
  mutate(oldnew = "Z")


# find unstable



responses <- bind_rows(incompletes,colorthings)

ggplot(data = responses %>% filter(exp == "SB2021_e2") %>% filter(condition == "D"),
       aes(y = oldnew,x = rating)) +
  geom_tile(aes(fill = factor(freq0)),color="black") +


  #geom_text(aes(label = round(delAIC)),size = 3.5) +
  scale_fill_manual(values = c("#ffff00","#000000","#ffffff","#ff0000","#00ff00","#0000ff"),
                    labels = c("no responses","responses",NA,"Gumbel > 2","GaussianEVSDT > 2",
                               "GaussianUVSDT > 2"))+
  scale_x_discrete(position = "top",name = "id") +
  scale_y_discrete(labels = c("New","Old","Delta AIC"))+

  facet_wrap(.~id) +
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

bestfit %>% filter(id == "SB2021_e2_9")
