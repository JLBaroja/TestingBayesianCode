rm(list=ls())
library('R2jags')

posterior_quantile <- function(mean_prior,sd_prior,
			       logmean_prior,logsd_prior){
	# Sample 'mean0' and 'sd0' from their priors
	mean0 <- rnorm(n=1,mean=mean_prior,sd=sd_prior)
	sd0 <- rlnorm(n=1,meanlog=logmean_prior,sdlog=logsd_prior)
	# Sample 'y' from the likelihood given prior sample 'th0'
	y <- rnorm(n=1,mean=mean0,sd=sd0)
	# Compute a (to-be-tested) posterior sample [th1,th2,...,thL]
	# given the sample from the likelihood
	n_obs <- length(y)
	observed <- list('y','mean_prior','sd_prior',
			 'logmean_prior','logsd_prior')
	unobserved <- c('mean','sd')
	write('model{	      
		mean~dnorm(mean_prior,1/sd_prior^2)
		sd~dlnorm(logmean_prior,1/logsd_prior^2)
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
