rm(list=ls())
library('R2jags')

posterior_quantile <- function(shape1_prior,shape2_prior){
	# Sample 'th0' from the prior
	th0 <- rbeta(n=1,shape1=shape1_prior,shape2=shape2_prior)
	# Sample 'y' from the likelihood given prior sample 'th0'
	y <- rbinom(n=1,size=20,prob=th0)
	# Compute a (to-be-tested) posterior sample [th1,th2,...,thL]
	# given the sample from the likelihood
	n_obs <- length(y)
	observed <- list('y','shape1_prior','shape2_prior')
	unobserved <- c('theta')
	write('model{	      
		theta~dbeta(shape1_prior,shape2_prior)
		y~dbin(theta,20)
		}','beta_binom.bug')
		bayes <- jags(data=observed,
			      parameters.to.save=unobserved,
			      model.file='beta_binom.bug')
		nds<- bayes$BUGSoutput$sims.list
	theta <- nds$theta # Posterior sample
	L <- dim(theta)[1] # Size of the posterior sample
	# Compute quantile (1/L)*sum(I_{th0>th_m})
	q0 <- sum(th0>theta)/L
}

quantiles <- NA
for(i in 1:1000){
	quantiles[i] <- posterior_quantile(5,3) 
}

hist(quantiles)
