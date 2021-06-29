test <- gumbelpar[1,] %>% .$c1 +
location <-gumbelpar[1,] %>% .$muo

scale <- 1

dgumbelR <- function(x, location = 0, scale = 1, log = FALSE)
  ### dgumbel in R
{
  q <- (x - location)/scale
  log.d <- -exp(-q) - q - log(scale)
  if (!log) exp(log.d) else log.d
}

dgumbel2R <- function(x, location = 0, scale = 1, log = FALSE)
{
  q <- (-x - location)/scale
  log.d <- -exp(-q) - q - log(scale)
  if (!log) exp(log.d) else log.d
}



(-critsprop[[5]] - location)/scale

-location - critsprop[[5]]
