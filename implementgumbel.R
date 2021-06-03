# checking ordinal package implementation of gumbel

x <- 0.5 # for small extremes -x
mu <- 1
beta <- 1

z <- (x - mu)/beta
dgumbel_manual <- 1/beta * exp(-(z + exp(-z)))
pgumbel_manual <- exp(-exp(-z))

# setting max = F is equivalent to setting x = -x with max = T
ordinal::dgumbel(x,location = mu, scale = beta,max=F)
ordinal::pgumbel(x,location = mu, scale = beta,max=T)
