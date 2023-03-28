rm(list=ls())
library('R2jags')

posterior_quantile <- function(mean_alpha_prior,sd_alpha_prior,
			       mean_beta_prior,sd_beta_prior,
			       shape_tau_prior,rate_tau_prior,
			       n_sessions=10){
	# Sample 'th0'=['alpha0', 'beta0', 'tau0'] from their priors
	alpha0 <- rnorm(n=1,mean=mean_alpha_prior,sd=sd_alpha_prior)
	beta0 <- rnorm(n=1,mean=mean_beta_prior,sd=sd_beta_prior)
	tau0 <- rgamma(n=1,shape=shape_tau_prior,rate=rate_tau_prior)
	sd0 <- 1/sqrt(tau0)
	#tau0 <- rlnorm(n=1,meanlog=shape_tau_prior,sdlog=1/sqrt(rate_tau_prior))
	# Sample data 'Br, Bl' from the likelihood given prior sample 'th0'
	Br <- NA
	Bl <- NA
	Wr <- rpois(n=n_sessions,lambda=50)+1 # Fixed data, avoid 0 reward counts
	Wl <- rpois(n=n_sessions,lambda=80)+1 # Fixed data, avoid 0 reward counts
	lambda_Br0 <<- NA
	lambda_Bl0 <<- NA
	for(i in 1:n_sessions){
		lambda_Br0[i] <<- rlnorm(n=1,alpha0/2+beta0*log(Wr[i])/2,1/sqrt(tau0))	 
		lambda_Bl0[i] <<- rlnorm(n=1,-alpha0/2-beta0*log(Wl[i])/2,1/sqrt(tau0))	 
		Br[i] <- rpois(n=1,lambda=lambda_Br0[i])+1
		Bl[i] <- rpois(n=1,lambda=lambda_Bl0[i])+1
	}
		# Compute a (to-be-tested) posterior sample [th1,th2,...,thL]
		# given the sample from the likelihood
		n_obs <- n_sessions
		observed <- list('Br','Bl',
				 'Wr','Wl',
				 'mean_alpha_prior','sd_alpha_prior',
				 'mean_beta_prior','sd_beta_prior',
				 'shape_tau_prior','rate_tau_prior',
				 'n_obs')
		unobserved <- c('alpha','beta','tau','sd','lambda_Br','lambda_Bl')
		write('model{	      
			alpha ~ dnorm(mean_alpha_prior,1/(sd_alpha_prior)^2)
			beta ~ dnorm(mean_beta_prior,1/(sd_beta_prior)^2)
			tau ~ dgamma(shape_tau_prior,rate_tau_prior)
			sd <- 1/pow(tau,1/2) 
			#tau ~ dlnorm(shape_tau_prior,rate_tau_prior)
			for(i in 1:n_obs){
				lambda_Br[i] ~ dlnorm(alpha/2+beta*log(Wr[i])/2,tau)
				lambda_Bl[i] ~ dlnorm(-alpha/2-beta*log(Wl[i])/2,tau)	 
				Br[i] ~ dpois(lambda_Br[i])
				Bl[i] ~ dpois(lambda_Bl[i])
			}
			}','generative_matching.bug')
			bayes <- jags(data=observed,
				      parameters.to.save=unobserved,
				      model.file='generative_matching.bug',
					inits=list(list(alpha=0,beta=0,tau=1),
					list(alpha=0,beta=0,tau=1),
					list(alpha=0,beta=0,tau=1)),
				      n.chains=3)
			nds<- bayes$BUGSoutput$sims.list
		alpha <- nds$alpha # Posterior sample
		beta <- nds$beta # Posterior sample
		tau <- nds$tau # Posterior sample
		sd <- nds$sd # Posterior sample
		L <- dim(alpha)[1] # Size of the posterior sample
		# Compute quantile (1/L)*sum(I_{th0>th_m})
		q0_alpha <- sum(alpha0>alpha)/L
		q0_beta <- sum(beta0>beta)/L
		q0_tau <- sum(tau0>tau)/L
		q0_sd <- sum(sd0>sd)/L
		return(list(q0_alpha=q0_alpha,q0_beta=q0_beta,q0_tau=q0_tau,q0_sd=q0_sd))
	}

	quantiles_alpha <- NA
	quantiles_beta <- NA
	quantiles_tau <- NA
	quantiles_sd <- NA
	for(i in 1:100){
		iteration <- posterior_quantile(mean_alpha_prior=0,sd_alpha_prior=1/sqrt(1),
						mean_beta_prior=0,sd_beta_prior=1/sqrt(1),
						shape_tau_prior=1,rate_tau_prior=1,
						#shape_tau_prior=5,rate_tau_prior=1,
						n_sessions=5)
		quantiles_alpha[i] <- iteration$q0_alpha 
		quantiles_beta[i] <- iteration$q0_beta 
		quantiles_tau[i] <- iteration$q0_tau 
		quantiles_sd[i] <- iteration$q0_sd 
		print(i)
	}

	layout(1:3)
	hist(quantiles_alpha,breaks=50)
	hist(quantiles_beta,breaks=50)
	hist(quantiles_tau,breaks=50)
