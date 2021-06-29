# Make Training/Test folds -----------------------------------------------------

source("preprocess_data.R")

trainforCV <- NULL

for(subjid in unique(testphase$id)){
  for(cond in c("A","B","C","D")){


ind <- testphase %>% filter(id == subjid) %>% filter(condition == cond)

folded <- ind %>% arrange(oldnew) %>%
  mutate(Fold = c(sample(factor(rep(c(1:10), length.out=60))),
                  sample(factor(rep(c(1:10), length.out=60)))))

dtrain <- NULL
for (j in c(1:10)) {

  dtrain1 <- folded %>%
    filter(Fold != j) %>%
    mutate(leftout = j) %>%
    mutate(cvid = paste0(id,"_cond",condition,"_LeftOutFold",leftout))

  dtrain <- bind_rows(dtrain,dtrain1)

}

trainforCV <- bind_rows(trainforCV,dtrain)

  }
}


trainforCV10Fold <- trainforCV %>% ungroup() %>%
  group_by(exp,condition,cvid) %>%
  group_nest(keep=T)

#saveRDS(trainforCV10Fold, file = "trainforCV10Fold.rds")

train10Fold <- readRDS("trainforCV10Fold.rds")

trainmainfold <- train10Fold %>% filter(.,grepl('LeftOutFold1',cvid)) %>% filter(.,!grepl('LeftOutFold10',cvid)) %>%
  unnest()
trainfinalfold <- train10Fold %>% filter(.,grepl('LeftOutFold10',cvid)) %>% unnest() %>% filter(Fold == 1)

test10Fold <- bind_rows(trainmainfold,trainfinalfold) %>% select(-exp1,-cvid1,-leftout,-cvid)

testforCV10Fold <- test10Fold %>%
  group_by(exp,id,condition) %>%
  group_nest(keep=T)

#saveRDS(testforCV10Fold, file = "testforCV10Fold.rds")


# Determine deviance for held-out folds ----------------------------------------

# Choose best run of 20 for each model and each training set
# raw data of all 20 fits on 220.222 server/backup harddrive

# load_files <- function(path,pattern) {
#
#   files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
#   # search for test data (dataTEST) file in all subdirectories (recursive) of path
#   tables <- lapply(files, readRDS)
#   do.call(bind_rows, tables)
#   # collate all of them into 1 dataframe
# }
#
# fits <- load_files("FitCV10Fold/","cond")
#
#
# bestfit <- fits %>%
#   mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
#   mutate(AIC = 2 * npar + 2 * objective) %>%
#   group_by(model,cvid,id,condition,leftout) %>%
#   arrange(objective) %>%
#   dplyr::slice(1)
#
# saveRDS(bestfit,file="CV10Fold_trainingfits.rds")

CVTraining <- readRDS("CV10Fold_trainingfits.rds")
CVTraining$CVLL <- NA
# Load data to be fitted

CVTest <- readRDS("testforCV10Fold.rds")


PredictCV <- NULL


for (subj in unique(CVTraining$id)){


  for (cond in c("A","B","C","D")){

    SubjTrain <- CVTraining %>% ungroup() %>%
      filter(id == subj) %>%
      filter(condition == cond)

    SubjTest <- CVTest %>% unnest() %>%
      filter(id == subj) %>%
      filter(condition == cond)

  for (modelN in unique(CVTraining$model)){ # add models to be predicted here


    holdouts <- unique(SubjTrain$leftout)

    for (i in seq_along(holdouts)){


      TrainData <- SubjTrain %>% filter(model ==  modelN) %>%
        filter(leftout == holdouts[[i]])

      if (nrow(TrainData) > 0){
        HoldOutData <- SubjTest %>% filter(Fold == i)

        TrainData$CVLL <- PredictCVHoldout(TrainData,modelN,HoldOutData)

        PredictCV <- bind_rows(PredictCV,TrainData)
      }
    }

  }

}
}

saveRDS(PredictCV,file = "CV10Fold_TestLL.rds")

