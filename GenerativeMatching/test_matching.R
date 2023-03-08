# Some refactoring

rm(list = ls())
library('R2jags')

## How many iterations?
nIterations <- 1000

## Collect all priors
priors <- list(mean_alpha = 0.0,
               sd_alpha   = 1.0,
               mean_beta  = 0.0,
               sd_beta    = 1.0,
               shape_tau  = 2.0,
               rate_tau   = 0.5)

## Collect information about the experiment
design <- list(n_sessions = 10)

# Get parameters from prior
get_parameters_from_prior <- function(priors) {
  alpha0 <- truncnorm::rtruncnorm(n    =  1, 
                                  a    = -2.0,
                                  b    =  2.0,
                                  mean = priors$mean_alpha, 
                                  sd   = priors$sd_alpha)
  beta0  <- truncnorm::rtruncnorm(n    =  1, 
                                  a    = -2.0,
                                  b    =  2.0,
                                  mean = priors$mean_beta,
                                  sd   = priors$sd_beta )
  tau0   <- rgamma(n     = 1, 
                   shape = priors$shape_tau,
                   rate  = priors$rate_tau)

  parameters <- list(alpha0 = alpha0, 
                     beta0  = beta0 ,
                     tau0   = tau0  )  

  return(parameters)
}

# Get data from parameters
get_data_from_parameters <- function(parameters, design) {
  # Sample data 'Br, Bl' from the likelihood given prior sample 'th0'
  Wr <- rpois(n = design$n_sessions, lambda = 50)
  Wl <- rpois(n = design$n_sessions, lambda = 80)
  
  rMean <-  parameters$alpha0 / 2 + parameters$beta0 * log(Wr) / 2
  lMean <- -parameters$alpha0 / 2 - parameters$beta0 * log(Wl) / 2
  stdev <- 1 / sqrt(parameters$tau0)
  
  lambda_Br0 <- rlnorm(n = design$n_sessions, meanlog = rMean, sdlog = stdev)
  lambda_Bl0 <- rlnorm(n = design$n_sessions, meanlog = lMean, sdlog = stdev)
  
  Br <- rpois(n = design$n_sessions, lambda = lambda_Br0)
  Bl <- rpois(n = design$n_sessions, lambda = lambda_Bl0)
  
  data <- list(Wr = Wr, Wl = Wl, Br = Br, Bl = Bl)
  
  return(data)
}


get_mcmc_from_data <- function(data, priors, design) {
  # Compute a (to-be-tested) posterior sample [th1,th2,...,thL]
  # given the sample from the likelihood
  observed <- list(
    Br               = data$Br,
    Bl               = data$Bl,
    Wr               = data$Wr,
    Wl               = data$Wl,
    mean_alpha_prior = priors$mean_alpha,
    sd_alpha_prior   = priors$sd_alpha,
    mean_beta_prior  = priors$mean_beta,
    sd_beta_prior    = priors$sd_beta,
    shape_tau_prior  = priors$shape_tau,
    rate_tau_prior   = priors$rate_tau,
    n_obs            = design$n_sessions
  )
  unobserved <- c('alpha', 'beta', 'tau', 'lambda_Br', 'lambda_Bl')
  write(
    'model{
         alpha ~ dnorm(mean_alpha_prior, pow(sd_alpha_prior, -2))T(-2,2)
         beta  ~ dnorm(mean_beta_prior , pow(sd_beta_prior , -2))T(-2,2)
         tau   ~ dgamma(shape_tau_prior, rate_tau_prior)T(0.01,)
         for(i in 1:n_obs){
             lambda_Br[i] ~ dlnorm( alpha/2 + beta * log(Wr[i])/2, tau)
             lambda_Bl[i] ~ dlnorm(-alpha/2 - beta * log(Wl[i])/2, tau)
             Br[i] ~ dpois(lambda_Br[i])
             Bl[i] ~ dpois(lambda_Bl[i])
         }
     }',
    'generative_matching.bug'
  )
  bayes <- jags(
    data = observed,
    parameters.to.save = unobserved,
    model.file = 'generative_matching.bug'
  )
  
  return(bayes)
}

get_quantiles_from_mcmc <- function(mcmc, parameters) {

  alpha <- mcmc$BUGSoutput$sims.list$alpha
  beta  <- mcmc$BUGSoutput$sims.list$beta
  tau   <- mcmc$BUGSoutput$sims.list$tau

  q0_alpha <- mean(parameters$alpha0 > alpha)
  q0_beta  <- mean(parameters$beta0  > beta )
  q0_tau   <- mean(parameters$tau0   > tau  )

  quantiles <- list(
    q0_alpha = q0_alpha,
    q0_beta  = q0_beta,
    q0_tau   = q0_tau
  )
  
  return(quantiles)
}


# Do the thing

## Compute quantiles
quantiles_alpha <- NA
quantiles_beta  <- NA
quantiles_tau   <- NA

for (i in 1:nIterations) {
  parameters <- get_parameters_from_prior(priors = priors)
  data       <- get_data_from_parameters (parameters = parameters, design = design)
  mcmc       <- get_mcmc_from_data       (data = data, priors = priors, design = design)
  quantiles  <- get_quantiles_from_mcmc  (mcmc = mcmc, parameters = parameters)

  quantiles_alpha[i] <- quantiles$q0_alpha
  quantiles_beta[i]  <- quantiles$q0_beta
  quantiles_tau[i]   <- quantiles$q0_tau
  print(i)
}

## Save cache
timestamp <- format(Sys.time() , "%Y%m%d_%H%M%S")
save(file = paste0('test_matching_', timestamp, '.RData'), 'quantiles_alpha', 'quantiles_beta', 'quantiles_tau')

## Plot histograms
pdf(file = paste0('test_matching_', timestamp, '.pdf'), width = 6, height = 4)
layout(1:3)
hist(quantiles_alpha, breaks = 50)
hist(quantiles_beta , breaks = 50)
hist(quantiles_tau  , breaks = 50)
dev.off()
