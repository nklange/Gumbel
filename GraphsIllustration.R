source("preprocess_data.R")
source("FitUVSDT.R")
library(cowplot)
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



# density graphs, based on median par estimates from fits to Spanton & Berry

generatedensities <- function(model,par){

  nitems <- 1e6
if(model == "GaussianUVSDT"){
      strength_old <- rnorm(nitems,mean = par[[1]],sd = par[[2]])
      strength_new <- rnorm(nitems,mean = 0, sd = 1)
  #     crits <- c(-Inf,cumsum(as.numeric(par[c(3:7)])),Inf)

   } else if(model == "GumbelEVSDT"){

       # switch the sign on muo to get equivalent results to the pgumbel fitting
       strength_old <- ordinal::rgumbel(nitems,location = -par[[1]],scale = 1,max=F)
       strength_new <- ordinal::rgumbel(nitems,location = 0,scale = 1,max=F)
  #    crits <- c(-Inf,cumsum(as.numeric(par[c(2:6)])),Inf)

   } else if(model == "ExGaussNormEVSDT"){

     strength_old <- brms::rexgaussian(nitems,mu = par[[1]],sigma = 1,beta=par[[2]])
     strength_new <- rnorm(nitems,mean = 0, sd = 1)

   } else if(model=="GaussianEVSDT"){

       strength_old <- rnorm(nitems,mean = par[[1]],sd = 1)
       strength_new <- rnorm(nitems,mean = 0, sd = 1)
   }

  return(list(old = strength_old,new=strength_new))
}

pars <- bestfit %>% ungroup() %>% filter(condition == "A") %>%
  select(model,muo,sigo,betao,c1,dc1,dc2,dc3,dc4) %>%
  group_by(model) %>%
  summarise_all(median,na.rm=T) %>% ungroup()


densityplot <- function(modelname){

par <- pars %>% filter(model==modelname) %>%
  select(names(get_start_par(modelname)))

crits <- par[(length(par)-4):length(par)]

genden <- generatedensities(modelname,par)


dats <- tibble(type = rep(c("old","new"),each=1e6),
               value = c(genden$old,genden$new))


critsprop <- cumsum(as.numeric(crits))
muo <- par[[1]]

if(modelname == "GumbelEVSDT"){
  #critsprop <- -critsprop
  muo <- -muo

}



plot <- ggplot(dats,aes(x = value,color=type))+
  geom_density(size=1) +
  coord_cartesian(xlim = c(-6,6))+
  scale_color_manual(labels=c("New","Studied"),values=c("#0072B2", "#D55E00"))+
  scale_x_continuous(breaks=c(0),labels=c("0"))+
  geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.2,0.5),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))

plot
}

plot_grid(densityplot("GaussianEVSDT"),densityplot("GaussianUVSDT"),
          densityplot("GumbelEVSDT"),densityplot("ExGaussNormEVSDT"),nrow=2,
          labels="AUTO",scale=.9)

#based on median par estimates from fits to Spanton & Berry
#ggsave("IllustrateModels.png", units="cm", width=30, height=20, dpi=600)
