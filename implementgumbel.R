# checking wikipedia parameterizations for pdf/cdf of gumbel and exgauss against
# packages

# Note: mu in gumbel is location parameter, mean is mu + scale*-digamma(1) (mascheroni constant)

# mu in brms::pexgauss is
# mu_input = location (mu of gaussian component) + scale (inverse rate, exp component)
# i.e., uses the mean of the distribution for input/output to the function

# checking ordinal package implementation of gumbel -----------------------------

x <- 0.5 # for small extremes -x
mu <- 1
beta <- 1

dgumbel_manual <- function(x, mu, beta) {
  z <- (x - mu) / beta
  1 / beta * exp(-(z + exp(-z))) #pdf
}

pgumbel_manual <- function(x, mu, beta) {
  z <- (x - mu) / beta
  exp(-exp(-z)) #cdf

}

# setting max = F is equivalent to setting x = -x with max = T
ordinal::dgumbel(x,location = mu, scale = beta,max=F)
ordinal::pgumbel(x,location = mu, scale = beta,max=T)


ordinal::pgumbel(-1,location = mu, scale = beta,max=T) -
  ordinal::pgumbel(-1.7, location = mu, scale = beta,max=T)

pgumbel_manual(-1,mu,beta) - pgumbel_manual(-1.7,mu,beta)

# exgaussian --------------------------------------------------------------------

x <- 0.5

eta <- 0.5
sigma <- 1
lambda <- 0.75
beta <- 1/lambda
mu <- eta + beta

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

dexgaus_manual(x,mu = eta,sigma,lambda) # standard parameterization
dexgaus_manual_brmspara(x,eta,sigma,beta)
brms::dexgaussian(x,mu,sigma,beta)


x = 0.5
eta <- 0.5
sigma <- 1
lambda <- 0.75
beta <- 1/lambda
mu <- eta + beta



pexgaus_manual <- function(x,mu,sigma,lambda){

  u <- lambda*(x-mu)
  v <- lambda * sigma

  pnorm(u,0,v) - exp(((-u + v^2)/2)+log(pnorm(u,v^2,v)))

}

# if x goes towards -Inf (i.e., gets towards -1e5-ish onwards), brms:pex evaluates to NaN
# in the manual implementation only seems to happen at -Inf exactly, otherwise
# evaluates to 0

pexgaus_manual(x,mu = eta,sigma,lambda) # standard parameterization
brms::pexgaussian(x,mu,sigma,beta)


