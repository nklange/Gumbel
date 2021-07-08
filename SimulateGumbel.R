library("RColorBrewer")
library("ggplot2")
library("cowplot")
library("stringr")


source("FitUVSDT.R")

calcDArms <- function(meanold,meannew,sdold,sdnew){

  (meanold - meannew)/sqrt(mean(c(sdold,sdnew))^2)

}

# check out a reasonable range of parameter estimates

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

fits <- load_files("Fits2/","fit")

genmodel <- "GumbelEVSDT"
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
  select(c(model,id,condition,muo,c1,dc1,dc2,dc3,dc4)) %>% ungroup() %>%
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
parameters <- tmvtnorm::rtmvnorm(n=2000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)
saveRDS(parameters,file="SimulationData/simulate_genGumbelEVSDT_parametervalues.rds")

# Genpar: ext mvn from par estimates --------------------------------------------
#rejectionsampling: for range of generating mu
# set up mu a bit to have a higher chance of high mu, increase sds for all parameters
# keep correlations
# then sample groups across range of 0 - 5 mu

correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
#means[1] <- means[1] + 1.5
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1]) * 2

sigmas <- sds %*% t(sds) * correlations


parameters2 <- tmvtnorm::rtmvnorm(n=100000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,0,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)

parametersext <- bind_rows(as_tibble(parameters2) %>% filter(V1 > 4.5 & V1 < 5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 4 & V1 < 4.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 3.5 & V1 < 4) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 3 & V1 < 3.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 2.5 & V1 < 3) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 2 & V1 < 2.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 1.5 & V1 < 2) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 1 & V1 < 1.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 0.5 & V1 < 1) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 > 0 & V1 < 0.5) %>% dplyr::slice(1:200))

saveRDS(parametersext,file="simulate_genExGaussNormEVSDT_parametervalues_mvnext_median.rds")

# Genpar: uniform --------------------------------------------------------------

nmatch <- 1000 #same number of datasets as mvtnorm total

parametersunif <- tibble(muo = runif(nmatch,min=-5,max=.Machine$double.eps),
       c1 = runif(nmatch,min=-3,max=1),
       dc1 = runif(nmatch,min=.01,2),
       dc2 = runif(nmatch,min=.01,2),
       dc3 = runif(nmatch,min=.01,2),
       dc4 = runif(nmatch,min=.01,2))

# Genpar: keep crit ------------------------------------------------------------

extparvalues<-readRDS(file="SimulationData/simulate_genExGaussNormEVSDT_parametervalues_mvnext.rds")

set.seed(1)
somecritsamples <- sample(1:2000,50)
keepcrits <- extparvalues[somecritsamples,] %>%
  set_colnames(names(get_start_par(genmodel))) %>%
  mutate(critnum = c(1:50)) %>%
  select(-muo) %>%
  dplyr::slice(rep(1:n(),each=26)) %>%
  mutate(muo = rep(seq(0,5,0.2),50),
         munum = rep(1:26,50)) %>%
  mutate(parid = paste0("genExGaussNormVSDT_extcrit",critnum,"_mu",munum) ) %>%
  select(-critnum,-munum)

saveRDS(keepcrits,file="SimulationData/simulate_genExGaussNorm_parametervalues_keepcrits.rds")

# Simulate data ----------------------------------------------------------------

simulation <- NULL
parametersext <- extparvalues#readRDS("SimulationData/simulate_genExGaussNorm_parametervalues_keepcrits.rds") %>%
 # relocate(muo)


for(i in c(1:dim(parametersext)[[1]])){

  for(j in c(2:51)){ # for each set of values, make 50 small-N datasets

id <- paste0("gen",genmodel,"_extpar",i,"_set",j)
#id <- paste0(parametersext[i,]$parid,"_set1")
parid <-paste0("gen",genmodel,"_extpar",i)# parametersext[i,]$parid#
pars <- as.numeric(parametersext[i,])#parameters[i,]


sim <- PredictSDT(data = NULL, model = "ExGaussNormEVSDT",par = pars,
           itemspertype = c(100,100))

# itemspertype choice as a LargeN (50000) and standard-experiment (100)

simtibble <- sim %>% mutate(id = id,
                            parid = parid,
               genmodel = genmodel)

# partibble <- tibble(value = pars,
#                     parameter = names(get_start_par(genmodel))) %>%
#   mutate(parid = parid,
#          genmodel = genmodel)

  # partibble <- pars %>%
  #   set_colnames(names(get_start_par(genmodel))) %>%
  # mutate(parid = parid,
  #        genmodel = genmodel)


simulation <- simulation %>% bind_rows(simtibble)
#parametertibble <- parametertibble %>% bind_rows(partibble)


}
}

saveRDS(simulation,file="simulate_genExGaussNormEVSDT_data_LOWN.rds")
#saveRDS(parametertibble,file="simulate_genGumbelEVSDT_parametervalues_mvnext.rds")

# Fitting routine for server --------------------------------------------------
# library("doParallel")
# library("foreach")



# doParallel::registerDoParallel(cores=16)
# mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
# model <- "GumbelEVSDT"
# model <- "GaussianUVSDT"
#
# res <- foreach(subjid = unique(simulation %>% .$id)
#                ,.combine = 'rbind'
#                ,.options.multicore=mcoptions) %dopar% {
#
#                  fitted <- FitSDT(simulation %>% filter(id=="genGumbelEVSDT_par1_set1"),
#                                   model,rep=20,freqdat=T)
#
#         saveRDS(fitted, file = paste0("FitsSimulation/", subjid, "_", model, ".rds"))
#                }


# Analysis ---------------------------------------------------------


# Large N - uniform #############################################################
# uniformly-dist. gen parameters
#


fits <- load_files("FitsSimulationunif/","gen")


bestfit <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  dplyr::slice(1)


# Analyse correlation r(mu_o, sig_o) in UVSDT fits of Gumbel-generated data

cor.test(bestfit %>% filter(model == "GaussianUVSDT") %>% .$muo,
         bestfit %>% filter(model == "GaussianUVSDT") %>% .$sigo)


linmodel0 <- lm(sigo ~ 1,data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel1 <- lm(sigo ~ poly(muo,1),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel2 <- lm(sigo ~ poly(muo,2),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel3 <- lm(sigo ~ poly(muo,3),data = bestfit %>% filter(model == "GaussianUVSDT"))
summary(linmodel)



anova(linmodel0, linmodel1)
anova(linmodel1, linmodel2)
anova(linmodel2, linmodel3)

logLik(linmodel1)
logLik(linmodel2)
logLik(linmodel3)


uvsdtmuosigo <- ggplot(bestfit %>% filter(model == "GaussianUVSDT"),aes(x=muo,y=sigo))+
  annotate("text", x = 0, y = 4, label = lm_eqn_poly1(bestfit %>% filter(model == "GaussianUVSDT")),
           parse = TRUE, color="black",hjust = 0)+
  annotate("text", x = 0, y = 3.8, label = lm_eqn_poly2(bestfit %>% filter(model == "GaussianUVSDT")),
           parse = TRUE, color="blue",hjust = 0)+
  annotate("text", x = 0, y = 3.6, label = lm_eqn_poly3(bestfit %>% filter(model == "GaussianUVSDT")),
           parse = TRUE, color="green",hjust = 0)+
  geom_point(alpha=0.1,color="red",size=2) +
  geom_smooth(method=lm,formula=y~x,color="black",size=1, se=T)+
  geom_smooth(method=lm,formula=y~poly(x,2),color="blue",size=1, se=T)+
  geom_smooth(method=lm,formula=y~poly(x,3),color="green",size=1, se=T)+

  scale_x_continuous(name=expression(paste(mu[o]," (estimated by UVSDT)")),limits=c(-0.5,10),breaks = c(0:10))+
  scale_y_continuous(name=expression(paste(sigma[o]," (estimated by UVSDT)")),limits=c(0,4))+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    #legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


parametervals <- readRDS(file="simulate_genGumbelEVSDT_parametervalues_unif.rds")


gumbelest <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

gumbelfit <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$message

uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

uvsdtfit <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$message


compareboths <- parametervals %>% select(parid,muo) %>% arrange(parid) %>% mutate(muo = -muo)

recoveringmu <- bind_rows(compareboths,compareboths) %>% mutate(muo_estimate = c(-gumbelest,uvsdtest),
                                                message = c(gumbelfit,uvsdtfit),
                                                model = c(rep("Gumbel",1000),rep("UVSDT",1000)))

relativeconverg_gen <- intersect(recoveringmu %>% filter(model== "Gumbel") %>% filter(message == "relative convergence (4)") %>%
  .$parid,recoveringmu %>% filter(model== "UVSDT") %>% filter(message == "relative convergence (4)") %>%
  .$parid)

muoestimateall <- ggplot(recoveringmu ,
       aes(y=muo_estimate,
           x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#009E73","#E69F00"))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(mu[o]," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", mu[o])),limits=c(-0.5,10),breaks = c(0:10))+
  annotate("text",x = 0,y=10,label=paste("1000/1000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.5),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


muoestimateconverge <- ggplot(recoveringmu %>% filter(parid %in% relativeconverg_gen),
       aes(y=muo_estimate,
           x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model", values = c("#009E73","#E69F00"))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(mu[o]," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", mu[o])),limits=c(-0.5,10),breaks = c(0:10))+
  annotate("text",x = 0,y=10,label=paste(length(relativeconverg_gen),"/1000 data sets - relative convergence (4)"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.5),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

uniformgenparameters <- plot_grid(uvsdtmuosigo,muoestimateall,muoestimateconverge,rel_widths = c(1,0.5,0.5),nrow=1)


# Large N - mvn #####################################################
# generated from reasonable Gumbel parameters (based on fits to Spanton & Berry)
# generated from multivariate normal of parameter space
# multivariate sample, no rejection sampling for range
#

fits <- load_files("FitsSimulation/","gen")

bestfit <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  dplyr::slice(1)

# Analyse correlation r(mu_o, sig_o) in UVSDT fits of Gumbel-generated data

cor.test(bestfit %>% filter(model == "GaussianUVSDT") %>% .$muo,
         bestfit %>% filter(model == "GaussianUVSDT") %>% .$sigo)


linmodel0 <- lm(sigo ~ 1,data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel1 <- lm(sigo ~ poly(muo,1),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel2 <- lm(sigo ~ poly(muo,2),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel3 <- lm(sigo ~ poly(muo,3),data = bestfit %>% filter(model == "GaussianUVSDT"))
summary(linmodel)



anova(linmodel0, linmodel1)
anova(linmodel1, linmodel2)
anova(linmodel2, linmodel3)

logLik(linmodel1)
logLik(linmodel2)
logLik(linmodel3)


uvsdtmuosigo <- ggplot(bestfit %>% filter(model == "GaussianUVSDT"),aes(x=muo,y=sigo))+
  annotate("text", x = 0, y = 4, label = lm_eqn_poly1(bestfit %>% filter(model == "GaussianUVSDT")),
           parse = TRUE, color="black",hjust = 0)+
  annotate("text", x = 0, y = 3.8, label = lm_eqn_poly2(bestfit %>% filter(model == "GaussianUVSDT")),
           parse = TRUE, color="blue",hjust = 0)+
  annotate("text", x = 0, y = 3.6, label = lm_eqn_poly3(bestfit %>% filter(model == "GaussianUVSDT")),
           parse = TRUE, color="green",hjust = 0)+
  geom_point(alpha=0.2,color="red",size=2) +
  geom_smooth(method=lm,formula=y~x,color="black",size=1, se=T)+
  geom_smooth(method=lm,formula=y~x + I(x^2),color="blue",size=1, se=T)+
  geom_smooth(method=lm,formula=y~x + I(x^2)  + I(x^3),color="green",size=1, se=T)+

  scale_x_continuous(name=expression(paste(mu[o]," (estimated by UVSDT)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste(sigma[o]," (estimated by UVSDT)")),limits=c(0,4))+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    #legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


parametervals <- readRDS(file="simulate_genGumbelEVSDT_parametervalues.rds")


gumbelest <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

gumbelfit <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$message

uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

uvsdtfit <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$message


compareboths <- parametervals %>% filter(parameter == "muo")%>% arrange(parid) %>% mutate(muo = -value) %>% select(parid,muo)

recoveringmu <- bind_rows(compareboths,compareboths) %>% mutate(muo_estimate = c(-gumbelest,uvsdtest),
                                                                message = c(gumbelfit,uvsdtfit),
                                                                model = c(rep("Gumbel",1000),rep("UVSDT",1000)))

relativeconverg_gen <- intersect(recoveringmu %>% filter(model== "Gumbel") %>% filter(message == "relative convergence (4)") %>%
                                   .$parid,recoveringmu %>% filter(model== "UVSDT") %>% filter(message == "relative convergence (4)") %>%
                                   .$parid)

muoestimateall <- ggplot(recoveringmu ,
                         aes(y=muo_estimate,
                             x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#009E73","#E69F00"))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(mu[o]," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", mu[o])),limits=c(-0.5,6),breaks = c(0:6))+
  annotate("text",x = 0,y=6,label=paste("1000/1000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.5),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


muoestimateconverge <- ggplot(recoveringmu %>% filter(parid %in% relativeconverg_gen),
                              aes(y=muo_estimate,
                                  x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model", values = c("#009E73","#E69F00"))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(mu[o]," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", mu[o])),limits=c(-0.5,6),breaks = c(0:6))+
  annotate("text",x = 0,y=6,label=paste(length(relativeconverg_gen),"/1000 data sets - relative convergence (4)"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.2,0.5),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

mvnstandardparameters <- plot_grid(uvsdtmuosigo,muoestimateall,muoestimateconverge,rel_widths = c(1,0.5,0.5),nrow=1)



#Large N - Ext Mvn #############################################################
# generated from reasonable Gumbel parameters (based on fits to Spanton & Berry)
# generated from multivariate normal of parameter space
# multivariate sample BUTsome changes to values + sds, correlation the same
# rejection sampling for range
# 10 bins, each with 200 data sets, to cover mu_gen of 0 - 5
# LARGE N
#



lm_eqn_poly3 <- function(df){
  m <- lm(sigo ~ muo + I(muo^2) + I(muo^3), df);
  eq <- substitute(sigma[o] == a + b %.% d*"'"* + c %.% d*"'"^2* + e %.% d*"'"^3*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        e = format(unname(coef(m)[4]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn_poly2 <- function(df){
  m <- lm(sigo ~ muo + I(muo^2), df);
  eq <- substitute(sigma[o] == a + b %.% d*"'"* + c %.% d*"'"^2*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn_poly1 <- function(df){
  m <- lm(sigo ~ muo, df);
  eq <- substitute(sigma[o] == a + b %.% d*"'"*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

bestfit <- readRDS("SimulationFits/genGumbelEVSDT_LARGEN_bestfits.rds")


# bestfit <- bind_rows(fits,fits2) %>%
#   mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
#   mutate(AIC = 2 * npar + 2 * objective) %>%
#   group_by(model,id) %>%
#   arrange(objective) %>%
#   dplyr::slice(1)



# Analyse correlation r(mu_o, sig_o) in UVSDT fits of Gumbel-generated data

cor.test(bestfit %>% filter(model == "GaussianUVSDT") %>% .$muo,
         bestfit %>% filter(model == "GaussianUVSDT") %>% .$sigo)


linmodel0 <- lm(sigo ~ 1,data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel1 <- lm(sigo ~ poly(muo,1),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel2 <- lm(sigo ~ poly(muo,2),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel3 <- lm(sigo ~ poly(muo,3),data = bestfit %>% filter(model == "GaussianUVSDT"))
summary(linmodel)



anova(linmodel0, linmodel1)
anova(linmodel1, linmodel2)
anova(linmodel2, linmodel3)

logLik(linmodel1)
logLik(linmodel2)
logLik(linmodel3)


uvsdtmuosigo <- ggplot(bestfit %>% filter(model == "GaussianUVSDT"),aes(x=muo,y=sigo))+
  # annotate("text", x = 0, y = 4, label = lm_eqn_poly1(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="black",hjust = 0)+
  # annotate("text", x = 0, y = 3.8, label = lm_eqn_poly2(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="blue",hjust = 0)+
  # annotate("text", x = 0, y = 3.6, label = lm_eqn_poly3(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="green",hjust = 0)+
  geom_point(alpha=0.2,color="#E69F00",size=2) +
  # geom_smooth(method=lm,formula=y~x,color="black",size=1, se=T)+
  # geom_smooth(method=lm,formula=y~x + I(x^2),color="blue",size=1, se=T)+
  # geom_smooth(method=lm,formula=y~x + I(x^2)  + I(x^3),color="green",size=1, se=T)+

  scale_x_continuous(name=expression(paste(mu[o], "(estimated by UVSDT)")),limits=c(-0.5,15),breaks = c(0:15))+
  scale_y_continuous(name=expression(paste(sigma[o]," (estimated by UVSDT)")),limits=c(0,8))+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    #legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


parametervals <- readRDS(file="SimulationData/simulate_genGumbelEVSDT_parametervalues_mvnext.rds")




exgaussbeta <- bestfit %>% filter(model == "ExGaussNormEVSDT") %>% select(id,betao) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$betao

exgaussest <- bestfit %>% filter(model == "ExGaussNormEVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo + exgaussbeta

exgausssigma <- sqrt(1 + 1/((1/exgaussbeta)^2))
gumbelest <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

gumbeluvsdtest <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

gumbelfit <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$message

uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$muo

uvsdtsigma <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,sigo) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$sigo

gumbeluvsdtsigma <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,betao) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% mutate(sigo = pi/sqrt(6)*betao) %>% .$sigo


uvsdtfit <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,message) %>% mutate(id = str_remove(id,"_set1")) %>%
  arrange(id) %>% .$message


uvsdtfitdprime <- calcDArms(uvsdtest,0,uvsdtsigma,1)
exgaussdprime <- calcDArms(exgaussest,0,exgausssigma,1)

compareboths <- parametervals %>% select(parid,muo) %>% arrange(parid) %>% mutate(muo = -muo)

modelorder <- c("Gumbel (free sigma[o])","Gumbel","UVSDT","ExGauss-Norm EVSDT")
recoveringmu <- bind_rows(compareboths,compareboths,
                          compareboths,compareboths) %>%
  mutate(muo_estimate = c(-gumbeluvsdtest,-gumbelest,uvsdtest,exgaussest),
                                                             #   message = c(gumbelfit,uvsdtfit),
                                                                model = c(rep("Gumbel (free sigma[o])",2000),
                                                                          rep("Gumbel",2000),
                                                                          rep("UVSDT",2000),
                                                                          rep("ExGauss-Norm EVSDT",2000))) %>%
  mutate(model = factor(model,levels=modelorder))

recoveringsigma <- bind_rows(compareboths,compareboths,
                          compareboths,compareboths) %>% mutate(sigo_estimate = c(gumbeluvsdtsigma,
                                                                    rep(pi/sqrt(6),2000),
                                                                    uvsdtsigma,
                                                                    exgausssigma),
                                                   #   message = c(gumbelfit,uvsdtfit),
                                                   model = c(rep("Gumbel (free sigma[o])",2000),
                                                             rep("Gumbel",2000),
                                                             rep("UVSDT",2000),
                                                             rep("ExGauss-Norm EVSDT",2000))
                                                   ) %>%
  mutate(model = factor(model,levels=modelorder))

recoveringdprime <- bind_rows(compareboths,compareboths) %>% mutate(dprime_estimate = c(uvsdtfitdprime,exgaussdprime),
                                                       #   message = c(gumbelfit,uvsdtfit),
                                                       model = c(rep("UVSDT",2000),
                                                                 rep("ExGauss-Norm EVSDT",2000))) %>%
  mutate(model = factor(model,levels=modelorder))


muoestimateall <- ggplot(recoveringmu ,
                         aes(y=muo_estimate,
                             x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00","#000000"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("Gumbel")),
                              expression(paste("UVSDT")),
                              expression(paste("ExGauss-Norm"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[o] - mu[n],"| (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated |", mu[o] - mu[n],"|")),limits=c(-0.5,15),breaks = c(0:15))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +

  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,

    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

ggplot(recoveringdprime ,
       aes(y=dprime_estimate,
           x=muo,color=model))+
  geom_point(aes(y=dprime_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#E69F00","#000000"),
                     labels=c(
                              expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[o] - mu[n],"| (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", d[a])),limits=c(-0.5,15),breaks = c(0:15))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +

  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,

    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

sigoestimateall <- ggplot(recoveringsigma,
                         aes(y=sigo_estimate,
                             x=muo))+
  geom_point(data = recoveringsigma %>% filter(model != "Gumbel"), aes(y=sigo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#E69F00","#000000"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("UVSDT")),
                              expression(paste("ExGauss-Norm"))))+
  geom_line(data =  recoveringsigma %>% filter(model == "Gumbel"),
            aes(y=sigo_estimate, x=muo,linetype=model),color="#009E73",size=1)+
  scale_linetype_manual(name = "Generating model",values="dashed")+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)),
         linetype = guide_legend(override.aes = list(size = 1,
                                                  alpha = 1) ) ) +
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste("|",mu[o] - mu[n],"| (Gumbel, generating value)")),limits=c(0,5.5),
                     breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", sigma[o])),limits=c(0,8),breaks = c(0:8))+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.1,0.7),
    legend.text.align = 0,
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))



LargeN <- plot_grid(muoestimateall,sigoestimateall,uvsdtmuosigo,rel_widths = c(1,1,1),nrow=1,
                    labels=c("1)","2)","3)"),scale=.95, label_x = 0, label_y = 0.9)


#ggsave("LargeN_Gumbelsimulation.png", units="cm", width=30, height=10, dpi=600)

# Low N- Ext Mvn ###############################################################
# generated from reasonable Gumbel parameters (based on fits to Spanton & Berry)
# generated from multivariate normal of parameter space
# multivariate sample BUTsome changes to values + sds, correlation the same
# rejection sampling for range
# 10 bins, each with 200 data sets, to cover mu_gen of 0 - 5
# LOW N
#


# Choose limited set of data points to display in plots

# # 2000 sets of parameters. Choose arbitrary 400 to avoid busy plot
#
set.seed(1)
#criticalpar <- paste0(rep("extpar",200),sample(1:2000,200))
criticalpar <- paste0(rep("extpar",2000),c(1:2000))

#
# # generated/fitted 50 sets per set of parameters
# # only show 10 sets per set to avoid busy plot
#
#
# criticalsets <- paste0(rep("set",10),sample(2:51,10))
#
# critsets <- paste0("genGumbelEVSDT_",rep(criticalpar,each=10),"_",criticalsets)
set.seed(1)
criticalsets <- paste0(rep("set",10),sample(2:51,10))

critsets <- paste0("genGumbelEVSDT_",rep(criticalpar,each=10),"_",criticalsets)




setdiff(critsets,bestfit %>% filter(model == "ExGaussNormEVSDT") %>% .$id)

# Analyse correlation r(mu_o, sig_o) in UVSDT fits of Gumbel-generated data

cor.test(bestfit %>% filter(model == "GaussianUVSDT") %>% .$muo,
         bestfit %>% filter(model == "GaussianUVSDT") %>% .$sigo)


linmodel0 <- lm(sigo ~ 1,data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel1 <- lm(sigo ~ poly(muo,1),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel2 <- lm(sigo ~ poly(muo,2),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel3 <- lm(sigo ~ poly(muo,3),data = bestfit %>% filter(model == "GaussianUVSDT"))
summary(linmodel)



anova(linmodel0, linmodel1)
anova(linmodel1, linmodel2)
anova(linmodel2, linmodel3)

logLik(linmodel1)
logLik(linmodel2)
logLik(linmodel3)



parametervals <- readRDS(file="simulate_genGumbelEVSDT_parametervalues_mvnext.rds") %>%
  filter(parid %in% paste0("genGumbelEVSDT_",criticalpar) )

ids <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo)  %>%
  arrange(id) %>% .$id

gumbelest <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo)  %>%
  arrange(id) %>% .$muo


gumbeluvsdtest <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,muo) %>%
  arrange(id) %>% .$muo

gumbelfit <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,message) %>%
  arrange(id) %>% .$message

gumbeluvsdtfit <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,message) %>%
  arrange(id) %>% .$message


uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>%
  arrange(id) %>% .$muo

uvsdtsigma <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,sigo) %>%
  arrange(id) %>% .$sigo

gumbeluvsdtsigma <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,betao) %>%
  arrange(id) %>% mutate(sigo = pi/sqrt(6)*betao) %>% .$sigo


uvsdtfit <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,message) %>%
  arrange(id) %>% .$message


compareboths <- parametervals %>% select(parid,muo) %>% arrange(parid) %>% mutate(muo = -muo) %>%
  dplyr::slice(rep(1:n(),each=10))

recoveringmu <- bind_rows(compareboths,compareboths,
                          compareboths) %>% mutate(muo_estimate = c(-gumbeluvsdtest,-gumbelest,uvsdtest),
                                                   message = c(gumbeluvsdtfit, gumbelfit,uvsdtfit),
                                                   id = c(ids,ids,ids),
                                                   model = c(rep("Gumbel (free sigma[o])",20000),
                                                             rep("Gumbel EVSDT",20000),
                                                             rep("UVSDT",20000)))

recoveringsigma <- bind_rows(compareboths,compareboths,
                             compareboths) %>% mutate(sigo_estimate = c(gumbeluvsdtsigma,
                                                                        rep(pi/sqrt(6),20000),
                                                                        uvsdtsigma),
                                                      message = c(gumbeluvsdtfit, gumbelfit,uvsdtfit),
                                                      id = c(ids,ids,ids),
                                                      model = c(rep("Gumbel (free sigma[o])",20000),
                                                                rep("Gumbel",20000),
                                                                rep("UVSDT",20000)))

sect1 <- intersect(recoveringmu %>% filter(model== "Gumbel EVSDT") %>% filter(message == "relative convergence (4)") %>%
                     .$id,
                   recoveringmu %>% filter(model== "Gumbel (free sigma[o])") %>% filter(message == "relative convergence (4)") %>%
                     .$id)
relativeconverg_gen <- intersect(sect1,
                                 recoveringmu %>% filter(model== "UVSDT") %>% filter(message == "relative convergence (4)") %>%
                                   .$id)
muoestimateall <- ggplot(recoveringmu %>% filter(id %in% relativeconverg_gen),
                         aes(y=muo_estimate,
                             x=muo,color=model))+
  # stat_binhex(aes(y=muo_estimate,
  #                           x=muo,fill=model),alpha=0.3) +
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.1) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("Gumbel")),
                              expression(paste("UVSDT"))))+
  # scale_fill_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00"),
  #                    labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
  #                             expression(paste("Gumbel")),
  #                             expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(-d,"'"," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", d,"'")),limits=c(-0.5,15),breaks = c(0:15))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +

  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,
    legend.background = element_rect(fill = "transparent"),

    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


ggplot(recoveringmu %>% filter(id %in% relativeconverg_gen),
       aes(y=muo_estimate,
           x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("Gumbel")),
                              expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(d,"'"," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", d,"'")),limits=c(-0.5,15),breaks = c(0:15))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +
  facet_grid(.~model)+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,

    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

sigoestimateall <- ggplot(recoveringsigma %>% filter(id %in% relativeconverg_gen),
                          aes(y=sigo_estimate,
                              x=muo))+
  geom_point(data = recoveringsigma %>% filter(model != "Gumbel") %>%
               filter(id %in% relativeconverg_gen), aes(y=sigo_estimate,
                                                                       x=muo,color=model),alpha=0.1) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),expression(paste("UVSDT"))))+
  geom_line(data =  recoveringsigma %>% filter(model == "Gumbel") %>% filter(id %in% relativeconverg_gen),
            aes(y=sigo_estimate, x=muo,linetype=model),color="#009E73",size=1)+
  scale_linetype_manual(name = "Generating model",values="dashed")+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)),
         linetype = guide_legend(override.aes = list(size = 1,
                                                     alpha = 1) ) ) +
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(-d,"'"," (Gumbel, generating value)")),limits=c(0,5.5),
                     breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", sigma[o])),limits=c(0,8),breaks = c(0:8))+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),

    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

ggplot(recoveringsigma %>% filter(id %in% relativeconverg_gen),
       aes(y=sigo_estimate,
           x=muo))+
  geom_point(data = recoveringsigma %>%
               filter(model != "Gumbel")%>%
               filter(id %in% relativeconverg_gen), aes(y=sigo_estimate,
                                                                       x=muo,color=model),alpha=0.1) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),expression(paste("UVSDT"))))+
  geom_line(data =  recoveringsigma %>% filter(model == "Gumbel") %>% filter(id %in% relativeconverg_gen),
            aes(y=sigo_estimate, x=muo,linetype=model),color="#009E73",size=1)+
  scale_linetype_manual(name = "Generating model",values="dashed")+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1)),
         linetype = guide_legend(override.aes = list(size = 1,
                                                     alpha = 1) ) ) +
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  facet_grid(.~model)+
  scale_x_continuous(name=expression(paste(-d,"'"," (Gumbel, generating value)")),limits=c(0,5.5),
                     breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", sigma[o])),limits=c(0,8),breaks = c(0:8))+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),

    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


uvsdtmuosigo <- ggplot(bestfit %>% filter(model == "GaussianUVSDT") %>% filter(id %in% relativeconverg_gen),aes(x=muo,y=sigo))+
  # annotate("text", x = 0, y = 4, label = lm_eqn_poly1(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="black",hjust = 0)+
  # annotate("text", x = 0, y = 3.8, label = lm_eqn_poly2(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="blue",hjust = 0)+
  # annotate("text", x = 0, y = 3.6, label = lm_eqn_poly3(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="green",hjust = 0)+
  geom_point(alpha=0.1,color="#E69F00",size=2) +
  # geom_smooth(method=lm,formula=y~x,color="black",size=1, se=T)+
  # geom_smooth(method=lm,formula=y~x + I(x^2),color="blue",size=1, se=T)+
  # geom_smooth(method=lm,formula=y~x + I(x^2)  + I(x^3),color="green",size=1, se=T)+

  scale_x_continuous(name=expression(paste("d' (estimated by UVSDT)")),limits=c(-0.5,15),breaks = c(0:15))+
  scale_y_continuous(name=expression(paste(sigma[o]," (estimated by UVSDT)")),limits=c(0,8))+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    #legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


LOWN <- plot_grid(muoestimateall,sigoestimateall,uvsdtmuosigo,rel_widths = c(1,1,1),nrow=1,
          labels=c("1)","2)","3)"),scale=.95, label_x = 0, label_y = 0.9)

# Large N - keep crit ##########################################################
# vary mu-only (0 to 5), while keeping crits identical
# based on multivariate-extended generated crits above
# LARGE N



fits <- load_files("FitSimulation_keepcrit_largeN/","gen")
parametervals <- readRDS(file="simulate_genGumbelEVSDT_parametervalues_keepcrits.rds")

critids <-rep(str_sort(rep(paste0(rep("genGumbelEVSDT_extcrit",50),c(1:50)),each=26)),4)


bestfit <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  dplyr::slice(1) %>%
  arrange(model,id)

bestfit2 <- readRDS("genGumbel_ExGauss_bestfit_KEEPCRITS.rds")

best <- bind_rows(bestfit,bestfit2) %>%  bind_cols(critid = critids)

best <- readRDS("genUVSDT_bestfit_KEEPCRITS.rds")
saveRDS(best,file="SimulationFits/genGaussianUVSDT_KEEPCRITS_bestfits.rds")

# Analyse correlation r(mu_o, sig_o) in UVSDT fits of Gumbel-generated data

cor.test(bestfit %>% filter(model == "GaussianUVSDT") %>% .$muo,
         bestfit %>% filter(model == "GaussianUVSDT") %>% .$sigo)


linmodel0 <- lm(sigo ~ 1,data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel1 <- lm(sigo ~ poly(muo,1),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel2 <- lm(sigo ~ poly(muo,2),data = bestfit %>% filter(model == "GaussianUVSDT"))
linmodel3 <- lm(sigo ~ poly(muo,3),data = bestfit %>% filter(model == "GaussianUVSDT"))
summary(linmodel)



anova(linmodel0, linmodel1)
anova(linmodel1, linmodel2)
anova(linmodel2, linmodel3)

logLik(linmodel1)
logLik(linmodel2)
logLik(linmodel3)





ids <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo)  %>%
  arrange(id) %>% .$id

gumbelest <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,muo)  %>%
  arrange(id) %>% .$muo


gumbeluvsdtest <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,muo) %>%
  arrange(id) %>% .$muo

gumbelfit <- bestfit %>% filter(model == "GumbelEVSDT") %>% select(id,message) %>%
  arrange(id) %>% .$message

gumbeluvsdtfit <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,message) %>%
  arrange(id) %>% .$message


uvsdtest <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,muo) %>%
  arrange(id) %>% .$muo

uvsdtsigma <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,sigo) %>%
  arrange(id) %>% .$sigo

gumbeluvsdtsigma <- bestfit %>% filter(model == "GumbelUVSDT") %>% select(id,betao) %>%
  arrange(id) %>% mutate(sigo = pi/sqrt(6)*betao) %>% .$sigo


uvsdtfit <- bestfit %>% filter(model == "GaussianUVSDT") %>% select(id,message) %>%
  arrange(id) %>% .$message


compareboths <- parametervals %>% select(parid,muo) %>% arrange(parid) %>%
  bind_cols(critid = critids[1:1300]) %>%  mutate(muo = -muo)

recoveringmu <- bind_rows(compareboths,compareboths,
                          compareboths) %>% mutate(muo_estimate = c(-gumbeluvsdtest,-gumbelest,uvsdtest),
                                                     # message = c(gumbeluvsdtfit, gumbelfit,uvsdtfit),
                                                  # id = c(ids,ids,ids),
                                                   model = c(rep("Gumbel (free sigma[o])",1300),
                                                             rep("Gumbel EVSDT",1300),
                                                             rep("UVSDT",1300)))

recoveringsigma <- bind_rows(compareboths,compareboths,
                             compareboths) %>% mutate(sigo_estimate = c(gumbeluvsdtsigma,
                                                                        rep(pi/sqrt(6),1300),
                                                                        uvsdtsigma),
                                                      #message = c(gumbeluvsdtfit, gumbelfit,uvsdtfit),
                                                     # id = c(ids,ids,ids),
                                                      model = c(rep("Gumbel (free sigma[o])",1300),
                                                                rep("Gumbel",1300),
                                                                rep("UVSDT",1300)))
#
# sect1 <- intersect(recoveringmu %>% filter(model== "Gumbel EVSDT") %>% filter(message == "relative convergence (4)") %>%
#                                    .$id,
#                                  recoveringmu %>% filter(model== "Gumbel (free sigma[o])") %>% filter(message == "relative convergence (4)") %>%
#                                    .$id)
# relativeconverg_gen <- intersect(sect1,
#                                  recoveringmu %>% filter(model== "UVSDT") %>% filter(message == "relative convergence (4)") %>%
#                                    .$id)
muoestimateall <- ggplot(recoveringmu,
                         aes(y=muo_estimate,
                             x=muo,color=model))+
  # stat_binhex(aes(y=muo_estimate,
  #                           x=muo,fill=model),alpha=0.3) +
  geom_line(aes(y=muo_estimate,
                x=muo,color=model,group=interaction(model,critid)),alpha=0.4)+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.2,size=1) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("Gumbel")),
                              expression(paste("UVSDT"))))+
  # scale_fill_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00"),
  #                    labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
  #                             expression(paste("Gumbel")),
  #                             expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(-d,"'"," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", d,"'")),limits=c(-0.5,15),breaks = c(0:15))+
  guides(color = guide_legend(override.aes = list(size = 1,
                                                  alpha = 1))) +

  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,
    legend.background = element_rect(fill = "transparent"),

    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


ggplot(recoveringmu %>% filter(id %in% relativeconverg_gen),
       aes(y=muo_estimate,
           x=muo,color=model))+
  geom_point(aes(y=muo_estimate,
                 x=muo,color=model),alpha=0.3) +
  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#009E73","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),
                              expression(paste("Gumbel")),
                              expression(paste("UVSDT"))))+
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(-d,"'"," (Gumbel, generating value)")),limits=c(-0.5,5.5),breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", d,"'")),limits=c(-0.5,15),breaks = c(0:15))+
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +
  facet_grid(.~model)+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,

    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

sigoestimateall <- ggplot(recoveringsigma ,
                          aes(y=sigo_estimate,
                              x=muo))+
  geom_line(data = recoveringsigma %>% filter(model != "Gumbel"), aes(y=sigo_estimate,
                                                                       x=muo,color=model,
                                                                      group=interaction(model,critid)),alpha=0.4) +
  geom_point(data = recoveringsigma %>% filter(model != "Gumbel"), aes(y=sigo_estimate,
                                                                       x=muo,color=model),alpha=0.2) +

  scale_color_manual(name = "Fitted model",values = c("#CC79A7","#E69F00"),
                     labels=c(expression(paste("Gumbel (free ", sigma[o],")")),expression(paste("UVSDT"))))+
  geom_line(data =  recoveringsigma %>% filter(model == "Gumbel"),
            aes(y=sigo_estimate, x=muo,linetype=model),color="#009E73",size=1)+
  scale_linetype_manual(name = "Generating model",values="dashed")+
  guides(color = guide_legend(override.aes = list(size = 1,
                                                  alpha = 1)),
         linetype = guide_legend(override.aes = list(size = 1,
                                                     alpha = 1) ) ) +
  # geom_point(data = compareboths, aes(x = muo, y = muo_est_gumbel), alpha=0.1,color="blue") +
  scale_x_continuous(name=expression(paste(-d,"'"," (Gumbel, generating value)")),limits=c(0,5.5),
                     breaks = c(0:5))+
  scale_y_continuous(name=expression(paste("estimated ", sigma[o])),limits=c(0,8),breaks = c(0:8))+
  #annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = c(0.3,0.7),
    legend.text.align = 0,
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),

    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


#  %>% filter(critid %in% unique(critids)[c(1:10)])

uvsdtmuosigo <- ggplot(bestfit %>% filter(model == "GaussianUVSDT"),
                       aes(x=muo,y=sigo, group=critid))+
  # annotate("text", x = 0, y = 4, label = lm_eqn_poly1(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="black",hjust = 0)+
  # annotate("text", x = 0, y = 3.8, label = lm_eqn_poly2(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="blue",hjust = 0)+
  # annotate("text", x = 0, y = 3.6, label = lm_eqn_poly3(bestfit %>% filter(model == "GaussianUVSDT")),
  #          parse = TRUE, color="green",hjust = 0)+
  geom_line(alpha=0.4,color="#E69F00",size=1)+
  geom_point(alpha=0.2,color="#E69F00",size=1) +
  # geom_smooth(method=lm,formula=y~x,color="black",size=1, se=T)+
  # geom_smooth(method=lm,formula=y~x + I(x^2),color="blue",size=1, se=T)+
  # geom_smooth(method=lm,formula=y~x + I(x^2)  + I(x^3),color="green",size=1, se=T)+

  scale_x_continuous(name=expression(paste("d' (estimated by UVSDT)")),limits=c(-0.5,15),breaks = c(0:15))+
  scale_y_continuous(name=expression(paste(sigma[o]," (estimated by UVSDT)")),limits=c(0,8))+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    #legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))


samecrit <- plot_grid(muoestimateall,sigoestimateall,uvsdtmuosigo,rel_widths = c(1,1,1),nrow=1,
                              labels=c("1)","2)","3)"),scale=0.95, label_x = 0, label_y = 0.9)


plot_grid(LargeN,LOWN,samecrit,labels=c("A","B","C"),nrow=3)


#ggsave("Gumbelsimulation_extmvn.png", units="cm", width=30, height=30, dpi=600)

# % Model mimicry ------------------------------------------------------

# LargeN - extmvn
fits <- load_files("FitSimulation_keepcrit_largeN/","gen")
# Gumbel: 81.9, GumbelUVSDT: 16.6, UVSDT: 1.46
# Gumbel: 98.2, UVSDT: 1.75

# LowN - extmvn
# Total (also non-converging fits included for all models)
# Gumbel: 81, GumbelUVSDT: 8.72, UVSDT: 10.2
# Gumbel: 88.8, UVSDT: 11.2

#set.seed(1)
#criticalsets <- paste0(rep("set",10),sample(2:51,10))

#critsets <- paste0("genGumbelEVSDT_",rep(criticalpar,each=10),"_",criticalsets)
DeltaAIC <- readRDS("genGumbel_trials200_bestfit.rds") %>%
  #filter(id %in% critsets) %>%
  arrange(model,id) %>%
  group_by(id) %>%
  mutate(delAIC = AIC - min(AIC))

# LargeN - keepcrits

fits <- bind_rows(load_files("FitsSimulationExt/","gen"),
                  load_files("FitSimulationGumbelUVSDT/","gen"))

# Gumbel: 84.6, GumbelUVSDT: 14.2, UVSDT: 1.25
# Gumbel: 98.5, UVSDT: 1.46



DeltaAIC <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  dplyr::slice(1) %>%
  arrange(model,id) %>%
  group_by(id) %>%
  mutate(delAIC = AIC - min(AIC))


DeltaAIC %>% filter(model == "GaussianUVSDT" & delAIC == 0) %>% .$id

DeltaAIC %>%
  #filter(model %in% c("GumbelEVSDT","GaussianUVSDT")) %>%
  filter(delAIC == 0) %>% # only 1 model can win
  group_by(model) %>%
  summarize(numBest = length(delAIC)) %>%
  complete(model,fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu =wu/sum(wu))

