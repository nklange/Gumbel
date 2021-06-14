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
saveRDS(trainforCV10Fold, file = "trainforCV10Fold.rds")
