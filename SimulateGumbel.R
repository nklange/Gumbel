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


#rejectionsampling: for range of generating mu
# set up mu a bit to have a higher chance of high mu, increase sds for all parameters
# keep correlations
# then sample groups across range of 0 - 5 mu

correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
means[1] <- means[1] - 1.5
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1]) * 2
sigmas <- sds %*% t(sds) * correlations


parameters2 <- tmvtnorm::rtmvnorm(n=100000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)

parametersext <- bind_rows(as_tibble(parameters2) %>% filter(V1 < -4.5 & V1 > -5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -4 & V1 > -4.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -3.5 & V1 > -4) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -3 & V1 > -3.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -2.5 & V1 > -3) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -2 & V1 > -2.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -1.5 & V1 > -2) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -1 & V1 > -1.5) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < -0.5 & V1 > -1) %>% dplyr::slice(1:200),
          as_tibble(parameters2) %>% filter(V1 < 0 & V1 > -0.5) %>% dplyr::slice(1:200))

# parameter values from uniform
nmatch <- 1000 #same number of datasets as mvtnorm total

parametersunif <- tibble(muo = runif(nmatch,min=-5,max=.Machine$double.eps),
       c1 = runif(nmatch,min=-3,max=1),
       dc1 = runif(nmatch,min=.01,2),
       dc2 = runif(nmatch,min=.01,2),
       dc3 = runif(nmatch,min=.01,2),
       dc4 = runif(nmatch,min=.01,2))

# Simulate from Gumbel and fit with UVSDT

genmodel <- "GumbelEVSDT"

simulation <- NULL
parametertibble <- NULL

for(i in c(1:dim(parametersext)[[1]])){

  #for(j in c(2:101)){

id <- paste0("gen",genmodel,"_extpar",i,"_set1")
#id <- paste0("gen",genmodel,"_unifpar",i,"_set1")
parid <- paste0("gen",genmodel,"_extpar",i)
pars <- parametersext[i,]#parameters[i,]


sim <- PredictSDT(data = NULL, model = "GumbelEVSDT",par = pars,
           itemspertype = c(5e4,5e4))

# itemspertype choice as high enough where parameter recovery is pretty spot-on.

simtibble <- sim %>% mutate(id = id,
                            parid = parid,
               genmodel = genmodel)

# partibble <- tibble(value = pars,
#                     parameter = names(get_start_par(genmodel))) %>%
#   mutate(parid = parid,
#          genmodel = genmodel)

  partibble <- pars %>%
    set_colnames(names(get_start_par(genmodel))) %>%
  mutate(parid = parid,
         genmodel = genmodel)


simulation <- simulation %>% bind_rows(simtibble)
parametertibble <- parametertibble %>% bind_rows(partibble)

#}

}

#saveRDS(simulation,file="simulate_genGumbelEVSDT_data_trials100k_mvnext.rds")
#saveRDS(parametertibble,file="simulate_genGumbelEVSDT_parametervalues_mvnext.rds")
library("doParallel")
library("foreach")

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


# Large-N simulations ----------------------------------------------------------
# generate from GumbelEVSDT. 50k old trials, 50k new trials per dataset
# 1000 sets of generating parameters in total

################################
# uniformly-dist. gen parameters
################################


fits <- load_files("FitsSimulationunif/","gen")


bestfit <- fits %>%
  mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  dplyr::slice(1)



lm_eqn_poly3 <- function(df){
  m <- lm(sigo ~ muo + I(muo^2) + I(muo^3), df);
  eq <- substitute(sigma[o] == a + b %.% mu[o] + c %.% mu[o]^2 + d %.% mu[o]^3*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        d = format(unname(coef(m)[4]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn_poly2 <- function(df){
  m <- lm(sigo ~ muo + I(muo^2), df);
  eq <- substitute(sigma[o] == a + b %.% mu[o] + c %.% mu[o]^2*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lm_eqn_poly1 <- function(df){
  m <- lm(sigo ~ muo, df);
  eq <- substitute(sigma[o] == a + b %.% mu[o]*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


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


################################################################################
# generated from reasonable Gumbel parameters (based on fits to Spanton & Berry)
# generated from multivariate normal of parameter space
# multivariate sample, no rejection sampling for range
################################################################################

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



################################################################################
# generated from reasonable Gumbel parameters (based on fits to Spanton & Berry)
# generated from multivariate normal of parameter space
# multivariate sample BUTsome changes to values + sds, correlation the same
# rejection sampling for range
# 10 bins, each with 200 data sets, to cover mu_gen of 0 - 5
################################################################################

fits <- load_files("FitsSimulationExt/","gen")

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


parametervals <- readRDS(file="simulate_genGumbelEVSDT_parametervalues_mvnext.rds")


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
                                                                model = c(rep("Gumbel",2000),rep("UVSDT",2000)))

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
  annotate("text",x = 0,y=10,label=paste("2000/2000 data sets"), hjust = 0)+
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
  annotate("text",x = 0,y=10,label=paste(length(relativeconverg_gen),"/2000 data sets - relative convergence (4)"), hjust = 0)+
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

mvnextparameters <- plot_grid(uvsdtmuosigo,muoestimateall,muoestimateconverge,rel_widths = c(1,0.5,0.5),nrow=1)

plot_grid(uniformgenparameters,
          mvnstandardparameters,mvnextparameters,
          nrow=3,labels=c("A","B","C"))

ggsave("DiffGenParameters_UVSDTfittedtoGumbel.png", units="cm", width=30, height=35, dpi=600)
