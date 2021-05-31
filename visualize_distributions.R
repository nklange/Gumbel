# Wordsetc






# Look at predicted curves for UVSDT wins --------------------------------------

fullsubj <- NULL
for (subjid in unique(testphase$id)){

  for(cond in c("A","B","C","D")){

    data <- testphase %>%  filter(id == subjid) %>%
      filter(condition == cond)

    dp <- prep_data(data)

    ROCs <- NULL
    partot <- NULL



    for (modelname in models){

      par <- bestfit %>%
        filter(model == modelname) %>%
        filter(id == subjid) %>%
        filter(condition == cond) %>%
        ungroup() %>%
        select(names(get_start_par(modelname)))

      makePredictions <- PredictSDT(data,modelname,as.numeric(par))

      predicted <- tibble(
        exp = dp$exp,
        id = dp$id,
        condition = cond,
        model = modelname,
        pred = "predicted",
        confidence = rep(dp$confidence,2),
        oldnew = c(rep("New",length(dp$confidence)),rep("Old",length(dp$confidence))),
        freq = makePredictions$predicted
      )

      ROCs <- ROCs %>% bind_rows(predicted)

    }



    observed <-  tibble(
      exp = dp$exp,
      id = dp$id,
      condition = cond,
      model = "Data",
      pred = "observed",
      confidence = rep(dp$confidence,2),
      oldnew = c(rep("New",length(dp$confidence)),rep("Old",length(dp$confidence))),
      freq =c(dp$datalist$New,dp$datalist$Old)
    )

    raw <- bind_rows(ROCs,observed) %>%
      mutate(confidence = factor(confidence,levels=c(6:1))) %>%
      group_by(exp,id,condition,model,pred,oldnew) %>%
      arrange(exp,id,condition,model,pred,oldnew,confidence) %>%
      mutate(roc = (cumsum(freq)/sum(freq))) %>%
      mutate(zroc =  qnorm(roc)) %>%
      mutate(zroc = ifelse(!is.finite(zroc),NA,zroc)) %>%
      mutate(condid = paste0(id,"_",condition))

    fullsubj <- fullsubj %>% bind_rows(raw)

  }

}

#saveRDS(fullsubj,file = "SB2021_individualpredictions.rds")

allfits <- bestfit %>%
  filter(model %in% models) %>%
  group_by(exp,id,condition) %>%
  mutate(delAIC = AIC - min(AIC)) %>%
  filter(condition == cond) %>%
  filter(exp == "SB2021_e1")

winning <- allfits %>%
  filter(delAIC == 0)

UVSDTWin <- winning %>% filter(model == "GaussianUVSDT" & delAIC == 0) %>%
  relocate(sigo, .after = muo) %>% .$id

GumbelWin <- winning %>% filter(model == "GumbelEVSDT" && delAIC == 0) %>% .$id




UVSDTpreds <- fullsubj %>%
  filter(id %in% UVSDTWin) %>%
  filer(condition %in% cond) %>%
  ungroup() %>% select(-exp,-id) %>%  unnest(c(ROCinfo)) %>%
  select(exp,id,condition,model,confidence,oldnew,zroc,roc) %>%
  pivot_wider(names_from = oldnew, values_from = c(roc,zroc))
Datonly <- UVSDTpreds %>% filter(model == "Data")

ggplot(Datonly,aes(x = zroc_New,y = zroc_Old)) +
  # geom_point(data = UVSDTpreds %>% filter(model != "Data"),
  #            aes(x = zroc_New,y = zroc_Old,color=model),size=0.5)+
  geom_smooth(data = UVSDTpreds %>% filter(model != "Data"),
              aes(x = zroc_New,y = zroc_Old,color=model),size=1,
              method="lm",formula = y ~ poly(x, 2), se = FALSE)+
  geom_point(data = Datonly,aes(x = zroc_New,y = zroc_Old), color= "black",size=2)+

  facet_wrap(.~id)



