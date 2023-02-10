rm(list=ls())
library('R2jags')

posterior_quantile <- function(mean_alpha_prior,sd_alpha_prior,
			       mean_beta_prior,sd_beta_prior,
			       shape_tau_prior,rate_tau_prior,
			       n_sessions=10){
	# Sample 'alpha0', 'beta0', and 'tau0' from their priors
	alpha0 <- rnorm(n=1,mean=mean_alpha_prior,sd=sd_alpha_prior)
	beta0 <- rlnorm(n=1,mean=mean_beta_prior,sd=sd_beta_prior)
	tau0 <- rgamma(n=1,shape=shape_tau_prior,rate=rate_tau_prior)
	# Sample 'y' from the likelihood given prior sample 'th0'
	Br <- NA
	Bl <- NA
	Wr <- rpois(n=n_sessions,lambda=10) # Fixed data
	Wl <- rpois(n=n_sessions,lambda=20) # Fixed data
	for(i in 1:n_sessions){
		lambda_Br[i] <- dlnorm(alpha0/2+beta0*log(Wr[i])/2,tau0)	 
		lambda_Bl[i] <- dlnorm(-alpha0/2-beta0*log(Wl[i])/2,tau0)	 
		Br[i] <- rpois(n=1,lambda=lambda_Br[i])
		Bl[i] <- rpois(n=1,lambda=lambda_Bl[i])
	}
	# Compute a (to-be-tested) posterior sample [th1,th2,...,thL]
	# given the sample from the likelihood
	n_obs <- n_sessions
	observed <- list('Br','Bl',
			 'mean_alpha_prior','sd_alpha_prior',
			 'mean_beta_prior','sd_beta_prior',
			 'sd_prior',
			 'logmean_prior','logsd_prior')
	unobserved <- c('mean','sd')
	write('model{	      
		mean~dnorm(mean_prior,1/sd_prior^2)
		sd~dnorm(logmean_prior,1/logsd_prior^2)
		y~dnorm(mean,1/sd^2)
		}','normal.bug')
		bayes <- jags(data=observed,
			      parameters.to.save=unobserved,
			      model.file='normal.bug')
		nds<- bayes$BUGSoutput$sims.list
	mean <- nds$mean # Posterior sample
	sd <- nds$sd # Posterior sample
	L <- dim(mean)[1] # Size of the posterior sample
	# Compute quantile (1/L)*sum(I_{th0>th_m})
	q0_mean <- sum(mean0>mean)/L
	q0_sd <- sum(sd0>sd)/L
	return(list(q0_mean=q0_mean,q0_sd=q0_sd))
}

quantiles_mean <- NA
quantiles_sd <- NA
for(i in 1:1000){
	iteration <- posterior_quantile(mean_prior=10,sd_prior=1,
				logmean_prior=5,logsd_prior=3)
	quantiles_mean[i] <- iteration$q0_mean 
	quantiles_sd[i] <- iteration$q0_sd 
}

layout(1:2)
hist(quantiles_mean)
hist(quantiles_sd)
