match_lnorm_meanlog <- function(mu, sigma){
  log(mu)-0.5*log((sigma^2 + mu^2)/mu^2)
}

match_lnorm_sdlog <- function(mu, sigma){
  sqrt(log((sigma^2 + mu^2)/mu^2))
}

match_beta_shape1 <- function(mu, sigma){
 (mu^2 - mu^3 - mu*sigma^2)/sigma^2
}

match_beta_shape2 <- function(mu, sigma){
  (mu - 2*mu^2 + mu^3 - sigma^2 + mu*sigma^2)/sigma^2
}

