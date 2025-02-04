library(ordinal) #for Gumbel with small extremes (max = F)
#library(brms) #for ExGaussian
#source("implementgumbel.R")

prep_data <- function(data,freqdat){

  if(freqdat){

    prepdat <- data %>%
      mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
      arrange(oldnew,response)

  } else {

    preprepdat <- data %>% group_by(oldnew,response) %>%
      summarize(Freq = length(response)) %>%
      ungroup() %>%
      mutate(oldnew = factor(oldnew,levels = c("New","Old")),
             response = factor(response, levels = c(1,2,3,4,5,6)))

    fullresp <- preprepdat %>% expand(oldnew,response)

    prepdat <- preprepdat %>% right_join(fullresp) %>%
      mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
      #mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      arrange(oldnew,response)

  }

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


  if (model %in% c("GaussianUVSDT")){
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

  } else if (model %in% c("GaussianEVSDT","GumbelLargeEVSDT","GumbelLargeNormEVSDT")){
    pstartunif <- tibble(
      muo    = runif(1, min=0, max=3),
      c1   = runif(1, min=-2, max=0),
      dc1  = runif(1, min=0, max=1),
      dc2  = runif(1, min=0, max=1),
      dc3  = runif(1, min=0, max=1),
      dc4  = runif(1, min=0, max=1)
    )

    pstart <- bind_rows(pstartunif)

  } else if (model %in% c("GumbelEVSDT","GumbelFlipEVSDT","GumbelNormEVSDT")){
    pstartunif <- tibble(
      muo    = runif(1, min=-3, max=0),
      c1   = runif(1, min=-2, max=0),
      dc1  = runif(1, min=0, max=1),
      dc2  = runif(1, min=0, max=1),
      dc3  = runif(1, min=0, max=1),
      dc4  = runif(1, min=0, max=1)
    )

    pstart <- bind_rows(pstartunif)

  }   else if (model %in% c("ExGaussNormEVSDT")){
    pstartunif <- tibble(
      muo   = runif(1, min=0, max=3),
      betao = runif(1, min=0.1, max=3),
      c1   = runif(1, min=-2, max=0),
      dc1  = runif(1, min=0, max=1),
      dc2  = runif(1, min=0, max=1),
      dc3  = runif(1, min=0, max=1),
      dc4  = runif(1, min=0, max=1)
    )

    pstart <- bind_rows(pstartunif)

  }  else if (model %in% c("GumbelUVSDT")){
    pstartunif <- tibble(
      muo   = runif(1, min=0, max=3),
      betao = runif(1, min=1, max=3),
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

  if (model %in% c("GaussianUVSDT")){
    # "d"    "sigo" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf, .Machine$double.eps,-Inf, .Machine$double.eps,
               .Machine$double.eps, .Machine$double.eps, .Machine$double.eps)
    upper <- Inf

  } else if (model %in% c("GaussianEVSDT","GumbelEVSDT","GumbelFlipEVSDT",
                          "GumbelNormEVSDT","GumbelLargeEVSDT","GumbelLargeNormEVSDT")){
    # "d"  "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf,-Inf, .Machine$double.eps,
               .Machine$double.eps, .Machine$double.eps, .Machine$double.eps)
    upper <- Inf

  } else if (model %in% c("ExGaussNormEVSDT","GumbelUVSDT")){
    # "d" "betao" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf,.Machine$double.eps,-Inf, .Machine$double.eps,
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
  } else if(model %in% c("GumbelUVSDT")){
    out <- gumbel_uvsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GumbelLargeEVSDT")){
    out <- gumbelLarge_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GaussianEVSDT")){
    out <- gaussian_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GumbelFlipEVSDT")){
    out <- gumbelFlip_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GumbelNormEVSDT")){
    out <- gumbelNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("GumbelLargeNormEVSDT")){
    out <- gumbelLargeNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("ExGaussNormEVSDT")){
    out <- exGaussNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  }

  return(-out)
}




gaussian_uvsdt_opt <- function(data_list,par,predictorLL){




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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  }

  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){

    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }


}
gaussian_evsdt_opt <- function(data_list,par,predictorLL){


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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  }

  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }

}
gumbel_uvsdt_opt <- function(data_list,par,predictorLL){

  d     <- par[1]
  betao  <- par[2]
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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=FALSE)-ordinal::pgumbel(I[i],location=0,scale=sign,max=FALSE)


  }


  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=betao,max=FALSE)-ordinal::pgumbel(I[i],location=d,scale=betao,max=FALSE)

  }

  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }
}
gumbel_evsdt_opt <- function(data_list,par,predictorLL){

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

  # likelihood by response category

  pNlikJ <- vector()
  pOlikJ <- vector()

  for (i in 1:(length(c)+1)){

    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=F)-ordinal::pgumbel(I[i],location=0,scale=sign,max=F)
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=F)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=F)

  }

  if(predictorLL == "LL"){

    Data <-  c(data_list$New,data_list$Old)
    llk <- Data * log(c(pNlikJ,pOlikJ))
    llk[Data == 0] <- 0

    llk <- sum(llk)
    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pNlikJ,pOlikJ))
  }
}
gumbelLarge_evsdt_opt <- function(data_list,par,predictorLL){


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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=T)-ordinal::pgumbel(I[i],location=0,scale=sign,max=T)

  }


  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=T)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=T)
  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }
} # large-extremes Gumbel
gumbelNorm_evsdt_opt <- function(data_list,par,predictorLL){

  d     <- par[1]
  betao  <- sqrt(6)/pi # since SD of Gumbel = pi/sqrt(6) * beta
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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],0,sign)-pnorm(I[i],0,sign)

  }

  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=betao,max=FALSE)-ordinal::pgumbel(I[i],location=d,scale=betao,max=FALSE)
  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }
} # new: pnorm, old: pgumbel small-extremes
gumbelLargeNorm_evsdt_opt <- function(data_list,par,predictorLL){

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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],0,sign)-pnorm(I[i],0,sign)

  }

  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=T)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=T)
  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }
} # new: pnrom, old: pgumbel large-extremes
gumbelFlip_evsdt_opt <- function(data_list,par,predictorLL){

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
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=TRUE)-ordinal::pgumbel(I[i],location=0,scale=sign,max=TRUE)

  }


  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=FALSE)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=FALSE)
  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)

    return(sum(c(NlikJ,OlikJ)))

  } else {
    return(c(pNlikJ,pOlikJ))
  }
} #new: pgumbel Large, old: pgumbel Small

exGaussNorm_evsdt_opt <- function(data_list,par,predictorLL){

  d     <- par[1]
  sigo <- 1
 # lambda <- 1/par[2]
 # sigo  <- sqrt(sigma +  1/lambda^2)

  betao <- par[2]
  c     <- vector()
  c[1]  <- par[3]
  c[2]  <- c[1] + par[4]
  c[3]  <- c[2] + par[5]
  c[4]  <- c[3] + par[6]
  c[5]  <- c[4] + par[7]

  # Constraints for equal-variance!
  sign  <- sigo
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
     pNlikJ[i] <- pnorm(I[i+1],0,sign)-pnorm(I[i],0,sign)

     # beta has to be > 0 (maybe try with super small beta -> super big lambda -> as gaussian
     # as possible)
    # evaluateboundaryi <- brms::pexgaussian(I[i],mu=0,sigma=sign,beta=0)
    # evaluateboundaryiplus <- brms::pexgaussian(I[i+1],mu=0,sigma=sign,beta=0)
    #
    # pNlikJ[i] <- ifelse(is.na(evaluateboundaryiplus),0,evaluateboundaryiplus) -
    #   ifelse(is.na(evaluateboundaryi),0,evaluateboundaryi)
  }




  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){


    evaluateboundaryi <- brms::pexgaussian(I[i],mu=d,sigma=sigo,beta=betao)
    evaluateboundaryiplus <- brms::pexgaussian(I[i+1],mu=d,sigma=sigo,beta=betao)

    pOlikJ[i] <- ifelse(is.na(evaluateboundaryiplus),0,evaluateboundaryiplus) -
      ifelse(is.na(evaluateboundaryi),0,evaluateboundaryi)

    #pOlikJ[i] <- evaluateboundaryiplus - evaluateboundaryi

    # this is a hack for quick test. pexGaussian(-Inf,mu,sigma,beta) evaluates to NaN while
    # pnorm(-Inf,mead,sd) evaluates to 0, so does pgumbel(-Inf,location,scale,max=F)
    # pexGaussian(Inf,mu,sigma,beta) evaluates to 1


  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)


    llk <- sum(c(NlikJ,OlikJ))

    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pNlikJ,pOlikJ))
  }
} #new: pnorm, old: exgaussian


fit_nlminb <- function(data, model, rep, startpar, freqdat){


  dp <- prep_data(data,freqdat)
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

FitSDT <- function (data, model, rep = rep, startpar = NULL, freqdat = F) {

  res <- fit_nlminb(data = data, rep = rep, model = model, startpar = startpar, freqdat = freqdat)

  return(res)
}

PredictCVHoldout <- function(TrainData, model, HoldOutData){

  pars <- TrainData %>% filter(model == model) %>%
    select(names(get_start_par(model))) %>% unlist()

  dp <- prep_data(HoldOutData,freqdat=F)

  out <- optfunction(model = model, data_list = dp$datalist,par = pars)
  return(out)

}


PredictSDT <- function(data = NULL, model, par, itemspertype = NULL){



  if (is.null(data)){


    # using rmultinom
    # alternative: using sample(x,n,prob,replace=T) but unsampled categories are not listed in results

    probs <- predict_frequencies(data_list = NULL, model = model, par = par)

    predfreq <- c(rmultinom(c(1:6),itemspertype[[1]],prob = probs[1:6]),
                  rmultinom(c(1:6),itemspertype[[2]],prob = probs[7:12]))

    simulated <- tibble(Freq = predfreq,
                        oldnew = rep(c("New","Old"),each = 6),
                        response = rep(c(1:6),2))

    return(simulated)

  } else {

    dp <- prep_data(data,freqdat=F)
    obsfreq<- c(dp$datalist$New,dp$datalist$Old)

    predfreq <- predict_frequencies(data_list = dp$datalist, model = model, par = par) *
      c(rep(dp$trialsNew,length(dp$confidence)),rep(dp$trialsOld,length(dp$confidence)))

    # Multiply each probability with the total number of new and old items
    # complicated here only to future-proof it for data sets with different numbers of targets/lures


    return(list(
      predicted = predfreq,
      observed = obsfreq
    ))

  }





}

predict_frequencies <- function(data_list,model,par){

  out <- vector()
  if(model %in% c("GaussianUVSDT")){
    out <- gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelEVSDT")){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelUVSDT")){
    out <- gumbel_uvsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelLargeEVSDT")){
    out <- gumbelLarge_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GaussianEVSDT")){
    out <- gaussian_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  }  else if(model %in% c("GumbelFlipEVSDT")){
    out <- gumbelFlip_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelNormEVSDT")){
    out <- gumbelNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelLargeNormEVSDT")){
    out <- gumbelLargeNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("ExGaussNormEVSDT")){
    out <- exGaussNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  }

  return(out)

}


GetGsquared <- function(model,LL,observedData){

  temp <- c(observedData[1:6] * log(observedData[1:6]/sum(observedData[1:6])),
            observedData[7:12] * log(observedData[7:12]/sum(observedData[7:12])))
  temp[observedData == 0] <- 0

  f <- tibble( GSq = 2*(LL - -sum(temp) )) # Calculate G^2
  f$df     <- length(observedData)-length(names(get_start_par(model)))-1
  f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq

  return(f)
}
#
#
# # sample_auc_withcritssimulation<-function(par,model,nitems){
# #
# #   if(model == "GaussianUVSDT"){
# #     strength_old <- rnorm(nitems,mean = par[[1]],sd = par[[2]])
# #     strength_new <- rnorm(nitems,mean = 0, sd = 1)
# #     crits <- c(-Inf,cumsum(as.numeric(par[c(3:7)])),Inf)
# #
# #   } else if(model == "GumbelEVSDT"){
# #
# #     # switch the sign on muo to get equivalent results to the pgumbel fitting
# #     strength_old <- ordinal::rgumbel(nitems,location = -par[[1]],scale = 1,max=F)
# #     strength_new <- ordinal::rgumbel(nitems,location = 0,scale = 1,max=F)
# #     crits <- c(-Inf,cumsum(as.numeric(par[c(2:6)])),Inf)
# #   }
# #
# #
# #
# #   rating_old <- vector()
# #   rating_new <- vector()
# #
# #   for(i in c(1:nitems)){
# #
# #
# #     rating_old[[i]] <- match(strength_old[[i]],sort(c(strength_old[[i]],crits)))-1
# #     rating_new[[i]] <- match(strength_new[[i]],sort(c(strength_new[[i]],crits)))-1
# #
# #   }
# #
# #   dp <- prep_data(tibble(oldnew = rep(c("Old","New"),each=nitems),
# #          response = c(rating_old,rating_new),
# #          id = "simulate"),freqdat=F)
# #
# #   rocs <-  tibble(
# #     id = dp$id,
# #     confidence = rep(dp$confidence,2),
# #     oldnew = c(rep("New",length(dp$confidence)),rep("Old",length(dp$confidence))),
# #     freq =c(dp$datalist$New,dp$datalist$Old)
# #   ) %>%
# #     mutate(confidence = factor(confidence,levels=c(6:1))) %>%
# #     group_by(id,oldnew) %>%
# #     arrange(id,oldnew,confidence) %>%
# #     mutate(roc = (cumsum(freq)/sum(freq)))
# #
# #   # taken from pROC package auc function
# #
# #   se <- rocs %>% filter(oldnew=="Old") %>% .$roc #roc$sensitivities #hit rate
# #   sp <- rocs %>% filter(oldnew=="New") %>% .$roc #roc$specificities #fa rate
# #
# #
# #   diffs.x <- sp[-1] - sp[-length(sp)]
# #   means.vert <- (se[-1] + se[-length(se)])/2
# #   auc <- sum(means.vert * diffs.x)
# #   auc
# # }
# #
# #
# # genpar <- bestfit %>% filter(model %in% c("GaussianUVSDT","GumbelEVSDT")) %>% arrange(id,condition,model) %>%
# #   select(id,model,muo,sigo,c1,dc1,dc2,dc3,dc4) %>% ungroup()
# #
# # gaussianpar <- genpar %>% filter(model=="GaussianUVSDT") %>% select(names(get_start_par("GaussianUVSDT")))
# # gumbelpar <- genpar %>% filter(model=="GumbelEVSDT") %>% select(names(get_start_par("GumbelEVSDT")))
# #
# #
# # multiples <- NULL
# #
# # for(j in c(1:dim(gaussianpar)[[1]])){
# #
# #   uvsdt_auc <- vector()
# #   gumbel_auc <- vector()
# #   for(i in c(1:100)){
# #
# #     uvsdt_auc[[i]] <- sample_auc_withcritssimulation(as.numeric(gaussianpar[j,]),"GaussianUVSDT",100)
# #     gumbel_auc[[i]] <- sample_auc_withcritssimulation(as.numeric(gumbelpar[j,]),"GumbelEVSDT",100)
# #
# #   }
# #
# #   results <- tibble(run = j,
# #                     model = c("GaussianUVSDT","GumbelEVSDT"),
# #                     mAUC = c(mean(uvsdt_auc),mean(gumbel_auc)),
# #                     sdAUC = c(sd(uvsdt_auc),sd(gumbel_auc)))
# #
# #
# #   multiples <- multiples %>% bind_rows(results)
# #
# # }
# #
# #
# #
# #
# sample_auc_critfree <-function(par,model){
#
#
#   if(model == "GaussianUVSDT"){
#
#     results             <- tibble(zfa = qnorm(seq(0.00001, 0.99999, 0.00001)))
#     results$zhit        <- 1/par[[2]] * results$zfa + par[[1]]/par[[2]]
#
#     res2 <- results %>% mutate(fa = pnorm(zfa), hit = pnorm(zhit))
#
#   } else if(model == "GumbelEVSDT"){
#
#     results             <- tibble(zfa = ordinal::qgumbel(seq(0.00001, 0.99999, 0.00001),
#                                                          location = 0, scale = 1,max=F))
#     results$zhit        <- results$zfa + par[[1]]
#
#     res2 <- results %>% mutate(fa = pnorm(qnorm(pgumbel(zfa))), hit = pnorm(qnorm(pgumbel(zhit)))) #?
#
#   } else if(model == "GaussianEVSDT"){
#
#       results             <- tibble(zfa = qnorm(seq(0.00001, 0.99999, 0.00001)))
#       results$zhit        <- results$zfa + par[[1]]
#
#       res2 <- results %>% mutate(fa = pnorm(zfa), hit = pnorm(zhit))
#
#   }
#
#
#   #taken from pROC package
#
#   se <- res2$hit#rocs %>% filter(oldnew=="Old") %>% .$roc #roc$sensitivities #hit rate
#   sp <- res2$fa#rocs %>% filter(oldnew=="New") %>% .$roc #roc$specificities #fa rate
#
#
#   diffs.x <- sp[-1] - sp[-length(sp)]
#   means.vert <- (se[-1] + se[-length(se)])/2
#   auc <- sum(means.vert * diffs.x)
#   auc
#
# }
# #
# #
# #
# SDcomp <- NULL
# for (j in c(1:500)){
#
# par<- tibble(muo = truncnorm::rtruncnorm(100,mean=2,sd=1,a=0,b=Inf),
#              sigo = truncnorm::rtruncnorm(100,mean=2,sd=1,a=1,b=3))
#
# auc<-vector()
# auc2 <- vector()
# auc3 <- vector()
#
# for (i in 1:100){
# auc[[i]] <- sample_auc_critfree(par[i,],"GaussianUVSDT")
# auc2[[i]] <- sample_auc_critfree(par[i,],"GumbelEVSDT")
# auc3[[i]] <- sample_auc_critfree(par[i,],"GaussianEVSDT")
# }
#
#
# collect <- tibble(run = j,
#                   model = c("GaussianUVSDT","GumbelEVSDT","GaussianEVSDT"),
#                   sds = c(sd(auc), sd(auc2), sd(auc3)),
#                   means = c(mean(auc), mean(auc2),mean(auc3)))
#
# SDcomp <- SDcomp %>% bind_rows(collect)
#
# }
#
#
# ggplot(SDcomp,aes(x = sdAUC,color=model))+
#   geom_density()
#
#
# # translating gumbel on gumbel coordinates (log-log) to gumbel on
# # normal coordinates
#
# # res3 <- res2 %>% mutate(zfa = qnorm(zfa), zhit =qnorm(zhit))
# #
# # ggplot(results,aes(zfa,zhit))+
# #   geom_point() +
# #   geom_abline(slope=1,intercept=0)+
# #   scale_y_continuous(name="pgumbel(hit = qgumbel(fa) + d)")+
# #   scale_x_continuous(name="qgumbel(fa)")+
# #   coord_cartesian(xlim=c(-4,4),ylim=c(-4,4))
# #
# #
# # ggplot(res2,aes(zfa,zhit))+
# #   geom_point() +
# #   geom_abline(slope=1,intercept=0)+
# #   scale_y_continuous(name="pgumbel(hit = qgumbel(fa) + d)")+
# #   scale_x_continuous(name="pgumbel(qgumbel(fa))")+
# #   coord_cartesian(xlim=c(0,1),ylim=c(0,1))
# #
# #
# # ggplot(res3,aes(zfa,zhit))+
# #   geom_line() +
# #   geom_abline(slope=1,intercept=0)+
# #   scale_y_continuous(name="qnorm(pgumbel(hit = qgumbel(fa) + d))")+
# #   scale_x_continuous(name="qnorm(pgumbel(qgumbel(fa)))")+
# #   coord_cartesian(xlim=c(-4,4),ylim=c(-4,4))
# #
# #
# # res2
# #
# # se <- res2$zhit# rocs %>% filter(oldnew=="Old") %>% .$roc #roc$sensitivities #hit rate
# # sp <- res2$zfa#rocs %>% filter(oldnew=="New") %>% .$roc #roc$specificities #fa rate
# #
# #
# # diffs.x <- sp[-1] - sp[-length(sp)]
# # means.vert <- (se[-1] + se[-length(se)])/2
# # auc <- sum(means.vert * diffs.x)
# # auc
