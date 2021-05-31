library(ordinal)

prep_data <- function(data){

  preprepdat <- data %>% group_by(oldnew,response) %>%
    summarize(Freq = length(response)) %>%
    mutate(oldnew = factor(oldnew,levels = c("New","Old")),
           response = factor(response, levels = c(1,2,3,4,5,6)))

  fullresp <- preprepdat %>% expand(response)

  prepdat <- preprepdat %>% right_join(fullresp) %>%
    mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
    arrange(oldnew,response)

  dl <- split(prepdat$Freq, f = prepdat$oldnew)
  trialsNew <- sum(prepdat %>% filter(oldnew == "New") %>% .$Freq)
  trialsOld <- sum(prepdat %>% filter(oldnew == "Old") %>% .$Freq)

  id <- unique(data$id)
  confidence <- sort(unique(prepdat$response))
  if (is.null(unique(data$cvid))) {
    cvid <- id
  }
  else {
    cvid <- unique(data$cvid)
  }
  if (is.null(unique(data$exp))) {
    exp <- NA
  }
  else {
    exp <- unique(data$exp)
  }
  if (is.null(unique(data$leftout))) {
    leftout <- 0
  }
  else {
    leftout <- unique(data$leftout)
  }

  return(list(datalist = dl, confidence = confidence, trialsNew = trialsNew, trialsOld = trialsOld, id = id,
              cvid = cvid, exp = exp, leftout = leftout))

}
get_start_par <- function(model){


  if (model %in% c("GaussianUVSDT","LaplaceUVSDT")){
    pstartunif <- tibble(
      muo   = runif(1, min=0, max=3),
      sigo = runif(1, min=1, max=4),
      c1   = runif(1, min=-2, max=0),
      dc1  = runif(1, min=0, max=1),
      dc2  = runif(1, min=0, max=1),
      dc3  = runif(1, min=0, max=1),
      dc4  = runif(1, min=0, max=1)
    )

    pstart <- bind_rows(pstartunif)

  } else if (model %in% c("GaussianEVSDT","GumbelEVSDT","GumbelFlipEVSDT","LaplaceEVSDT")){
    pstartunif <- tibble(
      muo    = runif(1, min=0, max=3),
      c1   = runif(1, min=-2, max=0),
      dc1  = runif(1, min=0, max=1),
      dc2  = runif(1, min=0, max=1),
      dc3  = runif(1, min=0, max=1),
      dc4  = runif(1, min=0, max=1)
    )

    pstart <- bind_rows(pstartunif)

  } else if (model %in% c("GumbelEVSDT","GumbelFlipEVSDT")){
    pstartunif <- tibble(
      muo    = runif(1, min=-3, max=0),
      c1   = runif(1, min=-2, max=0),
      dc1  = runif(1, min=0, max=1),
      dc2  = runif(1, min=0, max=1),
      dc3  = runif(1, min=0, max=1),
      dc4  = runif(1, min=0, max=1)
    )

    pstart <- bind_rows(pstartunif)

  }

 return(pstart)
}
get_par_limits <- function(model){

  if (model %in% c("GaussianUVSDT","LaplaceUVSDT")){
      # "d"    "sigo" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf, .Machine$double.eps,-Inf, .Machine$double.eps,
               .Machine$double.eps, .Machine$double.eps, .Machine$double.eps)
    upper <- Inf

  } else if (model %in% c("GaussianEVSDT","GumbelEVSDT","GumbelFlipEVSDT","LaplaceEVSDT")){
      # "d"  "c1"   "dc1"  "dc2"  "dc3"  "dc4"

      lower <- c(-Inf,-Inf, .Machine$double.eps,
                 .Machine$double.eps, .Machine$double.eps, .Machine$double.eps)
      upper <- Inf

  }

  list(lower,upper)
}

optfunction <- function(model, data_list, par){

  if(model %in% c("GaussianUVSDT")){
    out <- gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GumbelEVSDT")){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GumbelFlipEVSDT")){
    out <- gumbelFlip_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GaussianEVSDT")){
    out <- gaussian_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  }

  return(-out)
}


gaussian_uvsdt_opt <- function(data_list,par,predictorLL){

DataN <- data_list$New
DataO <- data_list$Old

d     <- par[1]
sigo  <- par[2]
c     <- vector()
c[1]  <- par[3]
c[2]  <- c[1] + par[4]
c[3]  <- c[2] + par[5]
c[4]  <- c[3] + par[6]
c[5]  <- c[4] + par[7]

# Constraints
sign  <- 1
I <- c(-Inf,c,Inf) # put criteria into larger array

# Likelihood of every trial
# New items
NlikJ <- vector()
pNlikJ <- vector()
for (i in 1:length(DataN)){
  pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  NlikJ[i] <- DataN[i] * log( pNlikJ[i] )
}

# Old items
OlikJ <- vector()
pOlikJ <- vector()
for (i in 1:length(DataO)){

  pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  OlikJ[i] <- DataO[i] * log( pOlikJ[i] )
}

if(predictorLL == "LL"){
  return(sum(c(NlikJ,OlikJ)))
} else {
  return(c(pNlikJ,pOlikJ))
}


}
gaussian_evsdt_opt <- function(data_list,par,predictorLL){

  DataN <- data_list$New
  DataO <- data_list$Old

  d     <- par[1]
  sigo  <- 1
  c     <- vector()
  c[1]  <- par[2]
  c[2]  <- c[1] + par[3]
  c[3]  <- c[2] + par[4]
  c[4]  <- c[3] + par[5]
  c[5]  <- c[4] + par[6]

  # Constraints
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  NlikJ <- vector()
  pNlikJ <- vector()
  for (i in 1:length(DataN)){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
    NlikJ[i] <- DataN[i] * log(pNlikJ[i] )
  }

  # Old items
  pOlikJ <- vector()
  OlikJ <- vector()
  for (i in 1:length(DataO)){
    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
    OlikJ[i] <- DataO[i] * log(pOlikJ[i])
  }

  if(predictorLL == "LL"){
    return(sum(c(NlikJ,OlikJ)))
  } else {
    return(c(pNlikJ,pOlikJ))
  }

}

gumbel_evsdt_opt <- function(data_list,par,predictorLL){

  DataN <- data_list$New
  DataO <- data_list$Old

  d     <- par[1]
  sigo  <- 1
  c     <- vector()
  c[1]  <- par[2]
  c[2]  <- c[1] + par[3]
  c[3]  <- c[2] + par[4]
  c[4]  <- c[3] + par[5]
  c[5]  <- c[4] + par[6]

  # Constraints
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  NlikJ <- vector()
  pNlikJ <- vector()
  for (i in 1:length(DataN)){
    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=FALSE)-ordinal::pgumbel(I[i],location=0,scale=sign,max=FALSE)
    NlikJ[i] <- DataN[i] * log(pNlikJ[i])

  }
ordinal::pgumbel

  # Old items
  OlikJ <- vector()
  pOlikJ <- vector()
  for (i in 1:length(DataO)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=FALSE)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=FALSE)
    OlikJ[i] <- DataO[i] * log( pOlikJ[i])
  }

  if(predictorLL == "LL"){
    return(sum(c(NlikJ,OlikJ)))
  } else {
    return(c(pNlikJ,pOlikJ))
  }
}
gumbelFlip_evsdt_opt <- function(data_list,par,predictorLL){

  DataN <- data_list$New
  DataO <- data_list$Old

  d     <- par[1]
  sigo  <- 1
  c     <- vector()
  c[1]  <- par[2]
  c[2]  <- c[1] + par[3]
  c[3]  <- c[2] + par[4]
  c[4]  <- c[3] + par[5]
  c[5]  <- c[4] + par[6]

  # Constraints
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  NlikJ <- vector()
  pNlikJ <- vector()
  for (i in 1:length(DataN)){
    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=TRUE)-ordinal::pgumbel(I[i],location=0,scale=sign,max=TRUE)
    NlikJ[i] <- DataN[i] * log(pNlikJ[i])

  }


  # Old items
  OlikJ <- vector()
  pOlikJ <- vector()
  for (i in 1:length(DataO)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=FALSE)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=FALSE)
    OlikJ[i] <- DataO[i] * log( pOlikJ[i])
  }

  if(predictorLL == "LL"){
    return(sum(c(NlikJ,OlikJ)))
  } else {
    return(c(pNlikJ,pOlikJ))
  }
}


fit_nlminb <- function(data, model, rep, startpar){


  dp <- prep_data(data)
  out_list <- vector("list",rep)

  if (is.null(startpar)) {
    startpar <- NULL
    for (reps in c(1:rep)) {
      startpar[[reps]] <- get_start_par(model)
    }
    startpar <- data.frame(matrix(unlist(startpar), nrow = length(startpar),
                                  byrow = T))
    colnames(startpar) <- names(get_start_par(model))
  }

  for (i in seq_len(rep)){


    start <- startpar %>% dplyr::slice(i) %>% unlist()
    tic <- Sys.time()
    tmp <- tryCatch(nlminb(start,objective = optfunction,
                           data_list = dp$datalist,
                           lower = get_par_limits(model)[[1]],
                           upper = get_par_limits(model)[[2]],
                           model = model,
                           control = list(eval.max = 300, iter.max = 300, trace = 1,
                                          rel.tol = 1e-05, x.tol = 1.5e-08)), error = function(e) NA)
    if (is.list(tmp)) {
      out_list[[i]] <- bind_cols(tidyr::spread(tibble::enframe(tmp$par),
                                               name, value) %>% setNames(sort(names(get_start_par(model)))),
                                 tibble::tibble(npar = length(names(get_start_par(model)))),
                                 tibble::as_tibble(tmp[c(2, 3, 4, 6)]), tibble::tibble(time = Sys.time() -
                                                                                         tic), tibble::tibble(model = model), tibble::tibble(id = dp$id),
                                 tibble::tibble(leftout = dp$leftout),
                                 tibble::tibble(rep = i),
                                 tibble::tibble(exp = dp$exp), tibble::tibble(cvid = dp$cvid),
                                 tidyr::spread(tibble::enframe(start, "start"),
                                               start, value) %>% setNames(paste0("s_", sort(names(get_start_par(model))))))
    }
  }
  dplyr::bind_rows(out_list)
}

FitSDT <- function (data, model, rep = rep, startpar = NULL) {
  res <- fit_nlminb(data = data, rep = rep, model = model, startpar = startpar)
  return(res)
}




PredictSDT <- function(data, model, par){

  dp <- prep_data(data)

  predfreq <- predict_frequencies(data_list = dp$datalist, model = modelname, par = par) *
    c(rep(dp$trialsNew,length(dp$confidence)),rep(dp$trialsOld,length(dp$confidence)))

  # Multiply each probability with the total number of new and old items
  # complicated here only to future-proof it for data sets with different numbers of targets/lures

  obsfreq<- c(dp$datalist$New,dp$datalist$Old)

  return(list(
    predicted = predfreq,
    observed = obsfreq
  ))
}

predict_frequencies <- function(data_list,model,par){

  out <- vector()
  if(model %in% c("GaussianUVSDT")){
    out <- gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("LaplaceUVSDT")){
    out <- laplace_uvsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelEVSDT")){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GaussianEVSDT")){
    out <- gaussian_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  }

  return(out)

}

#
# GetGsquared <- function(model,LL,observedData){
#
#   temp <- observedData * log(observedData/60)
#   temp[observedData == 0] <- 0
#
#   f <- tibble( GSq = 2*(LL - -sum(temp) )) # Calculate G^2 (better approximator to Chi-sq dist than Pearson Chi-Sq)
#   f$df     <- length(observedData)-length(names(get_start_par(model)))-1            # df model (L&F,p.63) n_observations - n_freeparas - 1
#   f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq
#
#   return(f)
# }
#
# predict_roc <- function(model,par){
#
#
#
#   if (model %in% c("GaussianUVSDT")){
#
#     par <- get_start_par("GaussianUVSDT")
#
#
#
#     results             <- tibble(zfa = stats::qnorm(c(seq(0.001, 0.99999, 0.01), 0.99999), 0, 1)
#     )
#   # DPSD # Yonelinas (1999), p. 1420
#
#   results$zfa[which(results$zfa > 1)] <- NA
#   results$zhit        <- 1/par[[2]] * results$zfa + par[[1]]/par[[2]]
#   results$zhit[results$zhit > 1] <- NA
#
#   } if (model %in% c("GumbelEVSDT")){
#
#
#     par <- get_start_par("GumbelEVSDT")
#
#
#
#     results             <- tibble(zfa = ordinal::qgumbel(c(seq(0.001, 0.99999, 0.01), 0.99999), location = 0, scale = 1,max=F)
#     )
#
#     scale <- 1
#     meanold <- -par[[1]] - scale*-digamma(1)
#     meannew <- 0 - scale*-digamma(1)
#     sig <- sqrt(((pi^2)/6) * scale^2)
#     results$zhit        <- results$zfa - (meanold - meannew)
#
#   }
#
#   # DPSD # Yonelinas (1999), p. 1420
#   criterion           <- stats::qnorm(results$fa/(1 - rl), - mul, sdl)
#   criterion[which(results$fa/(1 - rl) > 1)] <- 1000
#   results$hit        <- results$fa + rt + (1 - rt) * stats::pnorm(criterion, - mut, sdt) - (1 - rl) * stats::pnorm(criterion, - mul, sdl)
#   results$hit[results$hit > 1] <- 1
# }
#
# dats <- raw %>% .$ROCinfo %>% .[[1]] %>% filter(model == "Data")
# obs <- tibble(zfa = dats %>% filter(oldnew == "New") %>% .$zroc,
#               zhit = dats %>% filter(oldnew == "Old") %>% .$zroc)
#
# ggplot(data = results,aes(x = zfa, y = zhit))+
#   geom_line(data = results,aes(x = zfa, y = zhit),color="red") +
#   geom_line(data = resultsGauss,aes(x = zfa, y = zhit),color="blue") +
#   coord_cartesian(xlim = c(-3,2),ylim=c(-3,2))+
#   geom_point(data = obs, aes(x =zfa,
#              y = zhit), shape = 21, color="black")
#
