source("preprocess_data.R")
source("FitUVSDT.R")
library(cowplot)
models <- c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

fits <- load_files("Fits2/","fit")

bestfit <- fits %>%
  #mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
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

pars <- bestfit %>% ungroup() %>% filter(condition == "D") %>%
  select(model,muo,sigo,betao,c1,dc1,dc2,dc3,dc4) %>%
  group_by(model) %>%
  summarise_all(median,na.rm=T) %>% ungroup() %>%
  mutate(muo = muo * 2) %>%
  mutate(betao = betao * 1.5)





densityplot <- function(modelname){

par <- pars %>% filter(model==modelname) %>%
  select(names(get_start_par(modelname)))

crits <- par[(length(par)-4):length(par)]

genden <- generatedensities(modelname,par)


dats <- tibble(type = rep(c("old","new"),each=1e6),
               value = c(genden$old,genden$new)) %>%
  mutate(model = modelname)

}

dats <- bind_rows(densityplot("GaussianUVSDT"),
densityplot("GumbelEVSDT"),densityplot("ExGaussNormEVSDT")) %>%
  mutate(model = factor(model, levels = c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                        labels = c("UVSDT","Gumbel","ExGauss")))

#c("#E69F00", "#009E73", "#56B4E9")
Ill_UVSDT <- ggplot(dats %>% filter(model=="UVSDT"),aes(x = value,fill=type))+
  geom_density(size=1,alpha=0.7) +
  coord_cartesian(xlim = c(-6,8),ylim=c(0,0.4))+
  scale_fill_manual(values=c("#8c6100","#E69F00"),labels=c("New","Studied"))+
  scale_x_continuous(breaks=c(0),labels=c("0"))+
  annotate("text",x=2,y=0.38,hjust=0,size=3.5,label=expression(paste("Free: ", mu[s], ", ",sigma[s])))+
  annotate("text",x=2,y=0.4,hjust=0,size=3.5,label=expression(paste("Fixed: ", mu[n]," = 0, ", sigma[n]," = 1")))+
  # geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  facet_grid(.~model)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.2,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))

Ill_Gum <- ggplot(dats %>% filter(model=="Gumbel"),aes(x = value,fill=type))+
  geom_density(size=1,alpha=0.7) +
  coord_cartesian(xlim = c(-6,8),ylim=c(0,0.4))+
  scale_fill_manual(values=c("#005941","#009E73"),labels=c("New","Studied"))+
  scale_x_continuous(breaks=c(0),labels=c("0"))+
  annotate("text",x=2,y=0.38,hjust=0,size=3.5,label=expression(paste("Free: ", mu[s])))+
  annotate("text",x=2,y=0.4,hjust=0,size=3.5,label=expression(paste("Fixed: ", mu[n]," = 0, ", beta[n]," = ",beta[s]," = 1")))+
  # geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  facet_grid(.~model)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.2,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


Ill_EG <- ggplot(dats %>% filter(model=="ExGauss"),aes(x = value,fill=type))+
  geom_density(size=1,alpha=0.7) +
  coord_cartesian(xlim = c(-6,8),ylim=c(0,0.4))+
  scale_fill_manual(values=c("#2b5e7a","#56B4E9"),labels=c("New","Studied"))+
  scale_x_continuous(breaks=c(0),labels=c("0"))+
  annotate("text",x=2,y=0.38,hjust=0,size=3.5,label=expression(paste("Free: ", mu[s[G]], ", ",beta[s]," = 1/",lambda[s])))+
  annotate("text",x=2,y=0.4,hjust=0,size=3.5,label=expression(paste("Fixed: ", mu[n]," = 0, ", sigma[n]," = ",sigma[s[G]]," = 1")))+
  # geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  facet_grid(.~model)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.2,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))

#based on median par estimates from fits to Spanton & Berry

plot_grid(Ill_UVSDT,Ill_Gum,Ill_EG,nrow=1)


ggsave("Figures/IllustrateModels.png", units="cm", width=30, height=10, dpi=600)


print(citation(package="brms"), bibtex=TRUE)
