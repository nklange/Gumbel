# checking wikipedia parameterizations for pdf/cdf of gumbel and exgauss against
# packages

# Note: mu in gumbel is location parameter, mean is mu + scale*-digamma(1) (mascheroni constant)

# mu in brms::pexgauss is
# mu_input = location (mu of gaussian component) + scale (inverse rate, exp component)
# i.e., uses the mean of the distribution for input/output to the function

# checking ordinal package implementation of gumbel -----------------------------

# x <- 0.5 # for small extremes -x
# mu <- 1
# beta <- 1

dgumbel_manual <- function(x, mu, beta) {
  z <- (x - mu) / beta
  1 / beta * exp(-(z + exp(-z))) #pdf
}

pgumbel_manual <- function(x, mu, beta) {
  z <- (x - mu) / beta
  exp(-exp(-z)) #cdf

}


# setting max = F is equivalent to setting x = -x with max = T
# ordinal::dgumbel(x,location = mu, scale = beta,max=F)
# ordinal::pgumbel(x,location = mu, scale = beta,max=T)
#
#
# ordinal::pgumbel(-1,location = mu, scale = beta,max=T) -
#   ordinal::pgumbel(-1.7, location = mu, scale = beta,max=T)
#
# pgumbel_manual(-1,mu,beta) - pgumbel_manual(-1.7,mu,beta)

# exgaussian --------------------------------------------------------------------

# x <- 0.5
#
# eta <- 0.5
# sigma <- 1
# lambda <- 0.75
# beta <- 1/lambda
# mu <- eta + beta

dexgaus_manual <- function(x,mu,sigma,lambda){

  erfc_term <- (mu + lambda * sigma^2 - x)/(sqrt(2) * sigma)

  # below is identical to pracma::erfc(erfc_term)
  #integrand <- function(x){exp(-x^2)}
  #erfc <- 2/sqrt(pi) * integrate(integrand,lower=erfcx_term,upper=Inf)[[1]]

  lambda/2 * exp(lambda/2 * (2 * mu + lambda * sigma^2 - 2*x)) * pracma::erfc(erfc_term)

}


dexgaus_manual_brmspara <- function(x,eta,sigma,beta){

  erfc_term <- (eta + sigma^2/beta - x)/(sqrt(2)*sigma)


  1/(2*beta) * exp(1/(2*beta)*(2 * eta + sigma^2/beta - 2*x)) * pracma::erfc(erfc_term)


}


# brms uses a different parameterization...vignette("brms_families")
# mu = eta + lambda, to reflect the mean of the distribution
# beta = 1/lambda (scale = inverse rate)

# the below all evaluate to the same
#
# dexgaus_manual(x,mu = eta,sigma,lambda) # standard parameterization
# dexgaus_manual_brmspara(x,eta,sigma,beta)
# brms::dexgaussian(x,mu,sigma,beta)


# x = 0.5
# eta <- 0.5
# sigma <- 1
# lambda <- 0.75
# beta <- 1/lambda
# mu <- eta + beta



pexgaus_manual <- function(x,mu,sigma,lambda){

  u <- lambda*(x-mu)
  v <- lambda * sigma

  pnorm(u,0,v) - exp(((-u + v^2)/2)+log(pnorm(u,v^2,v)))

}

# if x goes towards -Inf (i.e., gets towards -1e5-ish onwards), brms:pex evaluates to NaN
# in the manual implementation only seems to happen at -Inf exactly, otherwise
# evaluates to 0

# pexgaus_manual(x,mu = eta,sigma,lambda) # standard parameterization
# brms::pexgaussian(x,mu,sigma,beta)


# WHat is going on with Gumbel? What does the best-fitting model look-like?

# Gumbel-EVSDT SB2021_e1_1, condition A
# bestfit :-1.22 -1.80  1.14 0.549 0.362 0.772
source("preprocess_data.R")

data <- testphase %>% filter(id == "SB2021_e1_1") %>% filter(condition=="A")
dp <- prep_data(data,freqdat=F)


genpar <- bestfit %>% filter(model %in% c("GaussianUVSDT","GumbelEVSDT")) %>% arrange(id,condition,model) %>%
  select(id,model,muo,sigo,c1,dc1,dc2,dc3,dc4) %>% ungroup()

gumbelpar <- genpar %>% filter(model=="GumbelEVSDT") %>% select(names(get_start_par("GumbelEVSDT")))
par <- gumbelpar[1,]
par[1] <- par[1] * -1
newpar <- c(as.numeric(par[1]),rev(cumsum(as.numeric(par[2:6])) * -1))

gumbel_evsdt_opt(dp$datalist,gumbelpar[1,],"predict")
gumbel_evsdt_opt(dp$datalist,gumbelpar[1,],"LL")
gumbelLarge_evsdt_opt(dp$datalist,gumbelpar[1,],"LL")


I <- c(-Inf,cumsum(as.numeric(par[2:6])),Inf) # put criteria into larger array

# Likelihood of every trial
# New items
pNlikJ <- vector()
for (i in 1:6){

  pNlikJ[i] <- ordinal::pgumbel(-I[i+1],location=0,scale=1,max=T,lower.tail=F)-
    ordinal::pgumbel(-I[i],location=0,scale=1,max=T,lower.tail=F)

}
ordinal::pgumbel(-1.8,location=0,scale=1,max=T)
pNlikJ <- vector()
for (i in 1:6){
  pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=1,max=F)-
    ordinal::pgumbel(I[i],location=0,scale=1,max=F)

}


# what happens when this gumbel model is fit:


pgumbel_manualMin <- function(q, location, scale, lower.tail=T,max=F) {

  q <- -q # for small-extremes
  z <- (q - location) / scale
  p <- exp(-exp(-z)) #cdf
  if (!lower.tail) p else 1-p # for small extremes

}

pgumbel_manualMax <- function(q, location, scale, lower.tail=T,max=T) {

  q <- q # for large-extremes
  z <- (q - location) / scale
  p <- exp(-exp(-z)) #cdf
  if (!lower.tail) 1-p else p # for small extremes


}



par <- gumbelpar[1,]

strength_old <- ordinal::rgumbel(1e5,location = -par[[1]],scale = 1,max=F)
strength_new <- ordinal::rgumbel(1e5,location = 0,scale = 1,max=F)
crits <- c(-Inf,cumsum(as.numeric(par[c(2:6)])),Inf)


rating_old <- vector()
rating_new <- vector()

for(i in c(1:1e5)){


  #if (model == "GaussianUVSDT"){
  rating_old[[i]] <-  match(strength_old[[i]],sort(c(strength_old[[i]],crits)))
  rating_new[[i]] <-  match(strength_new[[i]],sort(c(strength_new[[i]],crits)))

  # } else if (model == "GumbelEVSDT"){
  #
  #   rating_old[[i]] <- 8 - match(strength_old[[i]],sort(c(strength_old[[i]],crits)))
  #   rating_new[[i]] <- 8 - match(strength_new[[i]],sort(c(strength_new[[i]],crits)))
  #
  # }
}


  dats <- tibble(type = rep(c("old","new"),each=1e5),
                 value = c(rating_old,rating_new))

dats %>% group_by(type,value) %>% summarize(freq=length(value)) %>%
  group_by(type) %>% summarize(prop = freq/sum(freq))


x <- seq(-5,5,.01)
dars2 <- tibble(x = x,
       y = ordinal::pgumbel(x,0,1,max=F))

ordinal::pgumbel(-1.8,0,1,max=F)
ordinal::pgumbel(-Inf,0,1,max=F)

ggplot(dars2, aes(x = x,y=y))+
  geom_point() +
  geom_vline(xintercept = crits[[2]]) +
  geom_vline(xintercept = crits[[3]]) +
  geom_vline(xintercept = crits[[4]]) +
  geom_vline(xintercept = crits[[5]]) +
  geom_vline(xintercept = crits[[6]])


dats <- tibble(type = rep(c("old","new"),each=1e5),
               value = c(strength_old,strength_new))

ggplot(dats,aes(x = value,color=type))+
  geom_density() +
  geom_vline(xintercept = crits[[2]]) +
  geom_vline(xintercept = crits[[3]]) +
  geom_vline(xintercept = crits[[4]]) +
  geom_vline(xintercept = crits[[5]]) +
  geom_vline(xintercept = crits[[6]])



par <- as.numeric(gumbelpar[1,])
