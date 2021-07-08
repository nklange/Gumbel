# based on Wagenmakers et al 2004 - find optimal decision rule between models
# i.e. to what extent need which models be penalized to find the optimal decision rule
# in other words, how mimic-y are models even under the optimal decision rule

MimicryFits <- NULL
for (gen in c("genGumbelEVSDT","genGaussianUVSDT","genExGaussNormEVSDT")){


  fits <- readRDS(paste0("SimulationFits/",gen,"_mimicry_bestfits.rds")) %>%
  arrange(id,model) %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen)))

data <- readRDS(paste0("SimulationData/simulate_",gen,"_data_mimicry.rds")) %>%
  arrange(id)

LL <- fits %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen))) %>%
  group_by(genmodel,id,model) %>%
  group_nest(keep=T,.key="fit")

observedData <- data %>% group_by(genmodel,id) %>%
  dplyr::slice(rep(1:n(),each=3))

getGsquaredvals <- LL %>% bind_cols(observedData %>% ungroup() %>%  dplyr::select(data)) %>%
  mutate(Gsq = map2(.x = fit,.y = data, .f = ~GetGsquared(.x$objective,.y$Freq))) %>%
  mutate(objective = map(.x = fit,.f = ~mean(.x$objective))) %>%
  mutate(Dev = map(.x = fit,.f = ~mean(2*.x$objective))) %>%
  mutate(AIC = map(.x = fit,.f = ~mean(.x$AIC)))


MimicryFits <- MimicryFits %>% bind_rows(getGsquaredvals)

}

# individual conditions

fits <- readRDS(paste0("SimulationFits/conditiongen_bestfits.rds")) %>%
  arrange(id,model) %>%
  mutate(genmodel = case_when(grepl( "ExGaussNormEVSDT", id, fixed = TRUE) ~ "ExGaussNormEVSDT",
                              grepl( "GumbelEVSDT", id, fixed = TRUE) ~ "GumbelEVSDT",
                              grepl( "GaussianUVSDT", id, fixed = TRUE) ~ "GaussianUVSDT",
                              TRUE~"NA")) %>%
  mutate(condition = case_when(grepl( "cond_A", id, fixed = TRUE) ~ "A",
                            grepl( "cond_B", id, fixed = TRUE) ~ "B",
                            grepl( "cond_C", id, fixed = TRUE) ~ "C",
                            grepl( "cond_D", id, fixed = TRUE) ~ "D",
                            TRUE~"NA"))


MimicryFits <- NULL
for (gen in c("genGumbelEVSDT","genGaussianUVSDT","genExGaussNormEVSDT")){
  for(cond in c("A","B","C","D")){


  data <- readRDS(paste0("SimulationData/simulate_",gen,"_data_mimicry_cond",cond,".rds")) %>%
  arrange(id)

LL <- fits %>% filter(genmodel == substr(gen,start=4,stop=nchar(gen))) %>%
  filter(condition == cond) %>%
  group_by(genmodel,condition,id,model) %>%
  group_nest(keep=T,.key="fit")

observedData <- data %>% group_by(genmodel,condition,id) %>%
  dplyr::slice(rep(1:n(),each=3))

getGsquaredvals <- LL %>% bind_cols(observedData %>% ungroup() %>%  dplyr::select(data)) %>%
  mutate(Gsq = map2(.x = fit,.y = data, .f = ~GetGsquared(.x$objective,.y$Freq))) %>%
  mutate(objective = map(.x = fit,.f = ~mean(.x$objective))) %>%
  mutate(Dev = map(.x = fit,.f = ~mean(2*.x$objective))) %>%
  mutate(AIC = map(.x = fit,.f = ~mean(.x$AIC)))


MimicryFits <- MimicryFits %>% bind_rows(getGsquaredvals)
  }
}



# testdat <- MimicryFits %>% dplyr::select(genmodel,id,model,Gsq,objective,Dev,AIC) %>%
#   unnest(cols=c(Gsq,objective,Dev,AIC))

cond <- "D"
testdat <- MimicryFits %>% dplyr::select(genmodel,condition,id,model,Gsq,objective,Dev,AIC) %>%
  unnest(cols=c(Gsq,objective,Dev,AIC))

#excludeids<-unique(allfits %>% filter(message != "relative convergence (4)") %>% .$id)

#allfits <- allfits %>% filter(!id %in% excludeids)

testdat <- testdat %>% filter(model %in% c("GaussianUVSDT","ExGaussNormEVSDT",
                                                                      "GumbelEVSDT"),
                                                        genmodel %in% c("GaussianUVSDT",
                                                                        "ExGaussNormEVSDT",
                                                                        "GumbelEVSDT")) %>%
  dplyr::select(genmodel,condition,id,model,objective,AIC) %>%
  mutate(objective = 2*objective)



testdat <- testdat %>% ungroup() %>% filter(condition==cond) %>% dplyr::select(-condition)

custompenalty_pair <- function(data,par,modelpair){

  calc<-data %>%
    mutate(penobjective = case_when(model == modelpair[[2]] ~ objective + par,
                                                 model == modelpair[[1]] ~ objective)) %>%
    pivot_longer(cols = c("objective","AIC","penobjective"),names_to = "penalty",values_to="value") %>%
    group_by(genmodel,penalty,id) %>%
    mutate(winning = ifelse(value == min(value),1,0)) %>%
    group_by(genmodel,penalty,model) %>%
    summarize(num = sum(winning)) %>%
    group_by(genmodel,penalty) %>%
    mutate(winprop = num/sum(num))

  return(calc)

}
opt_penalty_maxmeanacc_pair <- function(data,par,modelpair){


  res <- custompenalty_pair(data,par,modelpair)

  meanpw <- sum(res %>% filter(penalty == "penobjective") %>%
                  filter(genmodel==model) %>% .$num)/
    sum(res %>% filter(penalty == "penobjective") %>% .$num)

  return(-meanpw)

}

custompenalty <- function(data,par){

  UVSDT_pen <- par[[1]]
  EG_pen <- par[[2]]
calc<-data %>% mutate(penobjective = case_when(model == "GaussianUVSDT" ~ objective + UVSDT_pen,
                                               model == "ExGaussNormEVSDT" ~ objective + EG_pen,
                                                model == "GumbelEVSDT" ~ objective)) %>%
  pivot_longer(cols = c("objective","AIC","penobjective"),names_to = "penalty",values_to="value") %>%
  group_by(genmodel,penalty,id) %>%
  mutate(winning = ifelse(value == min(value),1,0)) %>%
  group_by(genmodel,penalty,model) %>%
  summarize(num = sum(winning)) %>%
  group_by(genmodel,penalty) %>%
  mutate(winprop = num/sum(num))

return(calc)

}
opt_penalty_maxmeanacc <- function(data,par){

      res <- custompenalty(data,par)


      meanpw <- sum(res %>% filter(penalty == "penobjective") %>%
        filter(genmodel==model) %>% .$num)/
        sum(res %>% filter(penalty == "penobjective") %>% .$num)

  return(-meanpw)

}
# opt_penalty_equalrecovery <- function(data,par){
#
#   res <- custompenalty(data,par)
#
#   diffpw <- abs(diff(res %>% filter(penalty == "penobjective") %>%
#                        filter(genmodel==model) %>% .$winprop))
#
#
#   return(mean(diffpw))
#
# }


#resultsEqual <- NULL
resultsMaxAcc <- NULL

for(i in c(1:20)){

start <- as_tibble(t(c(runif(1,-1,3),runif(1,-1,3)))) %>% set_colnames(c("s_UVSDT_pen","s_ExGN_pen"))

# fits1 <- optim(par = start, opt_penalty_equalrecovery, data = testdat)
#
# res1 <- as_tibble(t(fits1$par)) %>% set_colnames(c("UVSDT_pen","ExGN_pen")) %>%
#   mutate(diffdiagonal = fits1$value,
#          fcteval = fits1$counts[[1]],
#          rep = i) %>%
#   bind_cols(start)
#
# resultsEqual <- resultsEqual %>% bind_rows(res1)

fits2 <- optim(par = start, opt_penalty_maxmeanacc, data = testdat)

res2 <- as_tibble(t(fits2$par)) %>% set_colnames(c("UVSDT_pen","ExGN_pen")) %>%
  mutate(meandiagonal = -fits2$value,
         fcteval = fits2$counts[[1]],
         rep = i) %>%
  bind_cols(start)

resultsMaxAcc <- resultsMaxAcc %>% bind_rows(res2)

}



saveRDS(resultsMaxAcc,file=paste0("SimulationFits/optimizeconfusionmatrix_meanacc_cond",cond,".rds"))
#saveRDS(resultsEqual,file="SimulationFits/optimizeconfusionmatrix_equalrecovery.rds")



# # both parameters vary from -1 to +4
# optgrid <- readRDS("optimizepenalty_gridsearch.rds")
#
# ggplot(optgrid,aes(x=c(1:90601),y=meandiff)) +
#   geom_point()
#
# ggplot(optgrid,aes(x=c(1:90601),y=meandiag)) +
#   geom_point()
#

# Optimize pairwise confusion matrices  ----------------------------------------

#A,B,C in order

models <- c("GumbelEVSDT","GaussianUVSDT")
models <- c("GumbelEVSDT","ExGaussNormEVSDT")
models <- c("ExGaussNormEVSDT","GaussianUVSDT")
models <- c("GumbelEVSDT","GaussianUVSDT","ExGaussNormEVSDT")

makeplotGOF <- function(i){

modelexpression <- list(c("GumbelEVSDT","GaussianUVSDT"),
                        c("GumbelEVSDT","ExGaussNormEVSDT"),
                        c("ExGaussNormEVSDT","GaussianUVSDT"))

modelnames <- list(c("Gumbel","UVSDT"),
                   c("Gumbel","EGNorm"),
                   c("EGNorm","UVSDT"))


prepdata<- testdat %>% filter(genmodel %in% modelexpression[[i]] &
                              model %in% modelexpression[[i]]) %>%
  mutate(model = factor(model,levels=modelexpression[[i]])) %>%
  group_by(genmodel,id) %>% arrange(genmodel,id,model) %>%
  mutate(genmodel= factor(genmodel,levels=modelexpression[[i]]))


dists <- prepdata %>%
  mutate(GOFAB = objective[[1]]-objective)  %>%
  filter(model != modelexpression[[i]][[1]])


#
# FgenA <- ecdf(dists %>% filter(genmodel == modelexpression[[i]][[1]]) %>% .$GOFAB)
# FgenB <- ecdf(dists %>%  filter(genmodel == modelexpression[[i]][[2]]) %>% .$GOFAB)
#
# z <- uniroot(function(z) FgenA(z) + FgenB(z) - 1,
#              interval<-c(min(dists$GOFAB),max(dists$GOFAB)))

# calcempopt <- dists %>% group_by(genmodel) %>%
#   summarize(sayA = length(GOFAB[GOFAB < z$root])/length(GOFAB),
#             sayB = length(GOFAB[GOFAB > z$root])/length(GOFAB)) %>%
#   pivot_longer(names_to="model",values_to="value",cols=c(sayA,sayB)) %>%
#   mutate(crit = paste0("optimal = ",round(z$root,2)),
#          crittype = "threshold")
#
#
# calcemp0 <- dists %>% group_by(genmodel) %>%
#   summarize(sayA = length(GOFAB[GOFAB < 0])/length(GOFAB),
#             sayB = length(GOFAB[GOFAB > 0])/length(GOFAB)) %>%
#   pivot_longer(names_to="model",values_to="value",cols=c(sayA,sayB)) %>%
#   mutate(crit = "nominal = 0",
#          crittype = "threshold")




maxmean <- optimize(f=opt_penalty_maxmeanacc_pair,data = prepdata,interval = c(-100,100),
                    modelpair=modelexpression[[i]])


naxm <- custompenalty_pair(prepdata,maxmean$minimum,modelexpression[[i]])

critmaxdiag <- naxm %>%
  mutate(model = case_when(model == modelexpression[[i]][[1]] ~ "sayA",
                           TRUE~"sayB")) %>% ungroup() %>% dplyr::select(-num) %>%
  rename("value" = winprop) %>%
  mutate(genmodel = as.character(genmodel))



Accuracy <- naxm %>% group_by(penalty) %>%
  filter(genmodel==model) %>%
  summarize(value = mean(winprop))%>%

  mutate(model = "sayA") %>%
  mutate(genmodel = "ZAcc")


Penalty <- naxm %>% group_by(penalty) %>% dplyr::select(penalty,model) %>% distinct() %>%
  mutate(value = case_when(model != "GumbelEVSDT" & penalty == "AIC" ~ ifelse(i == 3,0,2),
                             penalty == "objective" ~ 0,
                             penalty == "penobjective" & model != modelexpression[[i]][[1]] ~ maxmean$minimum,
                             TRUE ~ 0))%>%
  mutate(model = case_when(model == modelexpression[[i]][[1]] ~ "sayA",
                           TRUE~"sayB")) %>%
  mutate(genmodel = "YPen")


tablegraph <- bind_rows(critmaxdiag,Accuracy,Penalty) %>%
mutate(penalty = factor(penalty,levels=c("objective","AIC","penobjective"),
                        labels=c("-2LL","AIC","Max Acc"))) %>%
 # mutate(model = factor(model)) %>%
  mutate(genmodel = factor(genmodel,levels=c("ZAcc","YPen",rev(modelexpression[[i]])),
                           labels = c("Accuracy","Penalty",rev(modelnames[[i]])))) %>%
  mutate(colorunder = ifelse(!genmodel %in% modelnames[[i]],1,0)) %>%
  mutate(colorunder = factor(colorunder))

nameexpression <- c(expression(paste(Delta, GOF[GumUVSDT]," = ", GOF[Gumbel] - GOF[UVSDT])),
 expression(paste(Delta, GOF[GumEGNorm]," = ", GOF[Gumbel] - GOF[EGNorm])),
expression(paste(Delta, GOF[EGNormUVSDT]," = ", GOF[EGNorm] - GOF[UVSDT])))


placeopt <- which(sort(c(seq(-5,15,5),maxmean$minimum)) == maxmean$minimum)
prelabels <- c(-5,0,5,10,15)
xlabels <- append(prelabels,"\nopt",after=placeopt-1)

colorexpression <- list(c("#56B4E9","#E69F00" ),
                        c("#56B4E9", "#009E73"),
                        c("#009E73","#E69F00"))

A <- ggplot(dists,aes(x=GOFAB))+
    geom_histogram(data=dists,aes(x=GOFAB,fill=genmodel,y=..density..),
                  alpha=0.4,binwidth = 0.05,
                 position = "identity")+

  # geom_density_ridges(data=dists,aes(x=GOFAB,fill=genmodel,y = genmodel, height = stat(density)),
  #   stat = "binline", bins = 1000, scale = 1,
  #   draw_baseline = FALSE
  # )+

  scale_color_manual(name="Generating model",values = colorexpression[[i]],
                     labels=modelnames[[i]])+
  scale_fill_manual(name = "Generating model",values = colorexpression[[i]],
                    labels=modelnames[[i]])+

  # geom_density(data=dists,aes(x=GOFAB,color=genmodel),
  #              position = "identity",size=1,bw="ucv",adjust=2,n=1024)+
  geom_vline(size=0.5,color="black",xintercept=0)+
  geom_vline(size=0.5,linetype="dashed",color="black",xintercept=maxmean$minimum)+
  coord_cartesian(xlim=c(-7.5,15),ylim=c(0,1)) +
  annotate("text",label=paste0("Recover: ",modelnames[[i]][[1]]),x=-4,y=0.9)+
  annotate("text",label=paste0("Recover: ",modelnames[[i]][[2]]),x=8,y=0.9)+
  scale_x_continuous(name = nameexpression[i],breaks=sort(c(seq(-5,15,5),maxmean$minimum)),
                     labels=xlabels)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.7,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


B <- ggplot(tablegraph,aes(y = genmodel,x = model,fill=colorunder)) +
  geom_tile(color="white") +
  scale_fill_manual(values=c("#F0F0F0","white"))+
  geom_text(aes(label=ifelse(value != 0, ifelse(genmodel=="Accuracy" | abs(value)<.005,
                                                stringr::str_remove(round(value,3), "^0+"),
                                                stringr::str_remove(round(value,2), "^0+")),"0")),size=4)+
  scale_x_discrete(position = "top",name="Recovered",labels=modelnames[[i]])+
  scale_y_discrete(position = "left",
                   name="Generating")+
  facet_grid(.~penalty)+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(hjust=0.9),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour="black"),
        strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill = "transparent",
                                        colour="transparent"),
        strip.text = element_text(size = 12))

plot_grid(A,B,nrow=2,rel_heights=c(0.6,0.4))

}
library(cowplot)

pwplots <- plot_grid(makeplotGOF(2),
                     makeplotGOF(1),
          makeplotGOF(3),labels=c("A","B","C"),scale=.9,nrow=1,
          label_x = 0, label_y = 1)


# Make tablegraph for 3-model situation
# First: estimate penalties like in two-model case


maxacc <-readRDS(paste0("SimulationFits/optimizeconfusionmatrix_meanacc_cond",cond,".rds"))

optthresh <- maxacc %>% filter(meandiagonal == max(meandiagonal))
penaltybased <- custompenalty(testdat,c(optthresh %>% .$UVSDT_pen,
                                optthresh %>% .$ExGN_pen))

critmaxdiag <- penaltybased %>%
  mutate(model = case_when(model == "GumbelEVSDT" ~ "sayA",
                           model == "ExGaussNormEVSDT" ~ "sayB",
                           TRUE~"sayC")) %>% ungroup() %>% dplyr::select(-num) %>%
  rename("value" = winprop) %>%
  mutate(genmodel = as.character(genmodel))



Accuracy <- penaltybased %>% group_by(penalty) %>%
  filter(genmodel==model) %>%
  summarize(value = mean(winprop))%>%

  mutate(model = "sayA") %>%
  mutate(genmodel = "ZAcc")


Penalty <- penaltybased %>% group_by(penalty) %>% dplyr::select(penalty,model) %>% distinct() %>%
  mutate(value = case_when(model != "GumbelEVSDT" & penalty == "AIC" ~ 2,
                           penalty == "objective" ~ 0,
                           penalty == "penobjective" & model == "GaussianUVSDT" ~ optthresh$UVSDT_pen,
                           penalty == "penobjective" & model == "ExGaussNormEVSDT" ~ optthresh$ExGN_pen,
                           TRUE ~ 0))%>%
  mutate(model = case_when(model == "GumbelEVSDT" ~ "sayA",
                           model == "ExGaussNormEVSDT" ~ "sayB",
                           TRUE~"sayC"))  %>%
  mutate(genmodel = "YPen")



# knn
# riffing on multi-modelPBCM (MMPBCM)  Schultheis & Naidu (2014)
# CogSci proceedings
# https://escholarship.org/content/qt7544w9b0/qt7544w9b0.pdf

# don't use difference GOF but use dimensionality of GOF to find
# use k-nearest-neighbor


makeTibbleCM <- function(cm) {
  tibble(data.frame(num = c(as.numeric(cm[1,c(1:3)]),
                            as.numeric(cm[2,c(1:3)]),
                            as.numeric(cm[3,c(1:3)])),
                    genmodel = rep(row.names(cm),each=3),
                    model = rep(row.names(cm),3)
  ))

}

library(caret)

prep <- testdat %>%
  dplyr::select(-AIC) %>%
  pivot_wider(names_from="model",id_cols=c("genmodel","id"),values_from="objective")


# Splitting data into train
# and test data

datprep <- prep %>% mutate(split = sample(c("TRUE","FALSE"),size=length(dim(prep)[[1]]),replace=T,prob=c(0.7,0.3))) %>%
  mutate(genmodel = factor(genmodel)) %>% ungroup() %>% dplyr::select(-id,-split)

set.seed(300)
#Spliting data as training and test set. Using createDataPartition() function from caret
indxTrain <- createDataPartition(y = datprep$genmodel,p = 0.75,list = FALSE)
training <- datprep[indxTrain,] %>% ungroup()
testing <- datprep[-indxTrain,] %>% ungroup()

set.seed(400)
ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
knnFit <- train(genmodel ~ ., data = training, method = "knn", trControl = ctrl,
                preProcess = c("center","scale"), tuneLength = 20)

#Output of kNN fit
knnFit$results

knnPredict <- predict(knnFit,newdata = testing)
#Get the confusion matrix to see accuracy value and other parameter values
test <- confusionMatrix( testing$genmodel,knnPredict, dnn = c("Generated", "Recovered"))
cm <- test$table


knnacc <- tibble(value = mean(testing$genmodel==knnPredict)) %>%

                   mutate(model = "sayA") %>%
                   mutate(genmodel = "ZAcc") %>%
                    mutate(penalty = "kNN")

knntab <- makeTibbleCM(cm) %>% mutate(penalty = "kNN")

# LDA

library(MASS)

# Estimate preprocessing parameters
preproc.param <- training %>%
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(training)
test.transformed <- preproc.param %>% predict(testing)

model <- lda(genmodel~., data =train.transformed)
predictions <- model %>% predict(test.transformed)

cm <- table(test.transformed$genmodel,predictions$class)
ldatab <- makeTibbleCM(cm) %>% mutate(penalty = "LDA")

ldaacc <- tibble(value = mean(test.transformed$genmodel==predictions$class)) %>%

  mutate(model = "sayA") %>%
  mutate(genmodel = "ZAcc") %>%
  mutate(penalty = "LDA")

# QDA

model <- qda(genmodel~., data =train.transformed)
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$genmodel)

model

cm <- table(test.transformed$genmodel, predictions$class)
qdatab <- makeTibbleCM(cm) %>% mutate(penalty = "QDA")

qdaacc <- tibble(value =mean(predictions$class==test.transformed$genmodel)) %>%

  mutate(model = "sayA") %>%
  mutate(genmodel = "ZAcc") %>%
  mutate(penalty = "QDA")

# MDA

library(mda)
model <- mda(genmodel~., data =train.transformed)
predicted.classes <- model %>% predict(test.transformed)
# Model accuracy
mean(predicted.classes == test.transformed$genmodel)

cm <- table(test.transformed$genmodel, predicted.classes)
mdatab <- makeTibbleCM(cm) %>% mutate(penalty = "MDA")

mdaacc <- tibble(value =mean(predicted.classes == test.transformed$genmodel)) %>%

  mutate(model = "sayA") %>%
  mutate(genmodel = "ZAcc") %>%
  mutate(penalty = "MDA")


mlmax <- bind_rows(knntab,ldatab,qdatab,mdatab) %>% group_by(penalty,genmodel) %>%
  mutate(value = num/sum(num)) %>% dplyr::select(-num) %>%
  mutate(model = case_when(model == "GumbelEVSDT" ~ "sayA",
                           model == "ExGaussNormEVSDT" ~ "sayB",
                           TRUE~"sayC"))

mlaccuracy <- bind_rows(knnacc,ldaacc,qdaacc,mdaacc)



tablegraph <- bind_rows(critmaxdiag,Accuracy,Penalty,mlmax,mlaccuracy) %>%
  mutate(penalty = factor(penalty,levels=c("objective","AIC","penobjective","kNN","LDA","QDA","MDA"),
                          labels=c("-2LL","AIC","Max Acc","k(=17)NN","LDA","QDA","MDA"))) %>%
  # mutate(model = factor(model)) %>%
  mutate(genmodel = factor(genmodel,levels=c("ZAcc","YPen","GaussianUVSDT","ExGaussNormEVSDT","GumbelEVSDT"),
                           labels = c("Accuracy","Penalty","UVSDT","EGNorm","Gumbel"))) %>%
  mutate(colorunder = ifelse(!genmodel %in%c("UVSDT","EGNorm","Gumbel"),1,0)) %>%
  mutate(colorunder = factor(colorunder))






Dgraph <- ggplot(tablegraph,aes(y = genmodel,x = model,fill=colorunder)) +
  geom_tile(color="white") +
  scale_fill_manual(values=c("#E8E8E8","white"))+
  geom_text(aes(label=ifelse(value != 0, ifelse(genmodel=="Accuracy" | abs(value)<.005,
                                                stringr::str_remove(round(value,3), "^0+"),
                                                stringr::str_remove(round(value,2), "^0+")),"0")),size=4)+
  scale_x_discrete(position = "top",name="Recovered",labels=c("Gumbel","EGNorm","UVSDT"))+
  scale_y_discrete(position = "left",
                   name="Generating")+
  facet_grid(.~penalty)+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(hjust=0.75),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour="black"),
        strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill = "transparent",
                                        colour="transparent"),
        strip.text = element_text(size = 12))



#Visualize distributions

testdat <- MimicryFits %>% dplyr::select(genmodel,condition,id,model,Gsq,objective,Dev,AIC) %>%
  filter(condition == cond) %>% ungroup() %>% dplyr::select(-condition) %>%
  unnest(cols=c(Gsq,objective,Dev,AIC))

plotty <- testdat %>% mutate(logobj = log(Gsq)) %>% dplyr::select(-AIC,-objective,-Dev,-Gsq) %>%
  mutate(genmodel = factor(genmodel,levels=c("GumbelEVSDT","ExGaussNormEVSDT","GaussianUVSDT"),
                           labels=c("Gumbel","EGNorm","UVSDT"))) %>%
  pivot_wider(names_from="model",values_from="logobj") %>%
  relocate(GumbelEVSDT,before=ExGaussNormEVSDT) %>%
  rename("Fit: Gumbel"=GumbelEVSDT,
         "Fit: EGNorm" =ExGaussNormEVSDT,
         "Fit: UVSDT" =GaussianUVSDT)



library(GGally)
p1 <- ggpairs(plotty,columns = c(1,2,5),
              mapping = ggplot2::aes(color = genmodel),
              legend = 1,
              upper = list(continuous = "blank"),
              lower = list(continuous = wrap("points", alpha = 0.1),
                           combo = wrap("dot", alpha = 0.4)),
              diag = list(continuous = wrap("densityDiag", alpha = 0.3)),
              xlab = expression(paste("log(",G^2,")")),
              ylab = NULL,
              axisLabels = c("show"))

    for(i in 1:p1$nrow) {
    for(j in 1:p1$ncol){
      p1[i,j] <- p1[i,j] +
        scale_fill_manual(values=c("#56B4E9","#009E73", "#E69F00"), name= "Generating model") +
        scale_color_manual(values=c("#56B4E9","#009E73", "#E69F00"), name= "Generating model")
    }
  }

plotface <- p1 + theme(legend.position = "bottom",
           panel.background = element_rect(fill = "white"),

           panel.border = element_rect(colour = "black", fill = "transparent"),
           strip.background = element_rect(fill = "#F0F0F0",color="black"),
           strip.text = element_text(size=12))


plot_grid(pwplots, ggmatrix_gtable(plotface),Dgraph,labels=c("","D",""),nrow=3,
          rel_heights=c(0.5,0.6,0.2),scale=.9,
          label_x = 0.05, label_y = 1)

ggsave(paste0("MimicryFull_condition",cond,".png"), units="cm", width=45, height=35, dpi=600)
