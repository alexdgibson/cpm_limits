# cumulative_meta_analysis_discrimination.R
# run a cumulative Bayesian meta-analysis
# January 2025

# to run, the data should be named 'df'

# load libraries
library(dplyr)
library(nimble)
library(metamisc) # for valmeta function
library(ggplot2)

# test with the EuroSCORE data set
# data(EuroSCORE)
# df <- EuroSCORE

# load in the data

df <- read.csv(file = )

### part 1: complete a meta-analysis in metamisc ###
# Meta-analysis of the c-statistic (random effects)
fit_discr <- valmeta(cstat = ,
                     cstat.se = ,
                     cstat.cilb = ,
                     cstat.ciub = ,
                     cstat.cilv = 0.95,
                     N = ,
                     O = ,
                     slab = ,
                     data = df)

# check summary of the metamisc meta-analysis
fit_discr

# plot the standard meta-analysis through metamisc
plot(fit_discr)

# check variance against sample size
plot(log2(df$n), fit_discr$data$theta.se)
median(df$n) # median study size



### part 2: create a cumulative meta-analysis for the discrimination ###
# cumulative completed with NIMBLE in a bayesian framework

# set the NimbleCode to run the meta-analysis
## code; copied from metamisc, bugsmodels.r
code_discr <- nimbleCode({
  ## Likelihood
  for (i in 1:N){ # loop through studies
    theta[i] ~ dnorm(alpha[i], wsprec[i])
    alpha[i] ~ dnorm(mu.tobs, bsprec) # study-specific estimate centred on overall mean
    wsprec[i] <- 1/(theta.var[i]) # study-specific precision using observed variance
  }
  # priors
  bsTau ~ T(dt(0, (0.5^2), df = 3), 0, 10) # truncated, see paper for prior set-up
  mu.tobs ~ dnorm(0, 1/1000) # overall mean (logit)
  pred.tobs ~ dnorm(mu.tobs, bsprec) # prediction interval (logit)
  # back-transform to [0,1] and scalar:
  logit(mu.obs) <- mu.tobs
  logit(pred.obs) <- pred.tobs
  bsprec <- 1/(bsTau*bsTau) # transform SD to precision
})

# set the constants
constants_discr <- list(N = nrow(df))

# set the data from the metamisc valmeta function
data_discr <- list(theta.var = fit_discr$data$theta.se^2, # use results from valmeta | this is the variance, it is calculated from different data when not directly available
                   theta = fit_discr$data$theta) # this is the logit c-stat

## initial values
inits_discr <- list(bsTau = 1,
                    mu.tobs = 0)

# parameters to store
parms_discr = c('mu.obs', 'pred.obs', 'bsTau')


# specify the chain details from the protocol
niter = 100000 # one hundred thousand
thin = 5 # five
n.chains = 2 # two
seeds = c(1234,5678) # one seed per chain | specified in the protocol
n.burnin = 10000 # ten thousand


# create a list of data lists in cumulative to make into multiple models
# create an empty list to store data
data_list <- list()

# create n data lists in cumulative to complete meta-analysis each time
for(i in 1:nrow(df)){
  
  a <- as.numeric(data_discr[[1]][1:i])
  b <- as.numeric(data_discr[[2]][1:i])
  name <- paste('item', i, sep='')
  tmp <- list(theta.var = a, theta = b)
  data_list[[name]] <- tmp
}

# create n models in cumulative to then be run through NimbleMCMC
# create an empty list to store data
cum_discr <- list()

# create n lists of models to complete meta-analysis each time
for(i in 1:length(data_list)){
  
  model_discr <- nimbleModel(code_discr, 
                             data = data_list[[i]], 
                             inits = inits_discr, 
                             constants = list(N = i))
  name <- paste('item', i, sep='')
  tmp <- list(model_1 = model_discr)
  cum_discr[[name]] <- tmp
}

# run the meta-analysis for each model in NimbleMCMC
# create an empty list to store outputs
cum_meta_discr <- list()

# run each cumulative model through MCMC
for(i in 1:length(cum_discr)){
  
  mcmc_out_discr <- nimbleMCMC(model = cum_discr[[i]][[1]],
                               inits = inits_discr,
                               monitors = parms_discr,
                               niter = niter, # times 2 for burn-in 
                               thin = thin,
                               nchains = n.chains, 
                               nburnin = n.burnin,
                               summary = TRUE, 
                               setSeed = seeds,
                               WAIC = FALSE)
  
  name <- paste('meta', i, sep='')
  tmp <- list(model = mcmc_out_discr)
  cum_meta_discr[[name]] <- tmp
}


# check that the meta-analysis was completed for each of the cumulative models
for (i in 1:length(cum_meta_discr)){
  print(cum_meta_discr[[i]]$model$summary$all.chains)
}

# take the output of all chains from the model
# create an empty list to store the outputs
cum_discr_plot <- list()

# take the outputs of each model
for (i in 1:length(cum_meta_discr)){
  
  mu <- cum_meta_discr[[i]]$model$summary$all.chains
  mu_df <- as.data.frame(mu)
  name <- paste('study', i, sep='')
  tmp <- list(summary = mu_df)
  cum_discr_plot[[name]] <- tmp
}


# from all the data take mu and confidence intervals to then be plotted
# create an empty list to store data
cum_discr_forest_plot <- data.frame()

# take mu and confidence intervals
for (i in 1:length(cum_discr_plot)){
  
  df <- as.data.frame(cum_discr_plot[[i]][[1]][2,1:5])
  cum_discr_forest_plot <- rbind(cum_discr_forest_plot, df)
}

# plot a forest plot for the discrimination
# for each of the cumulative models plot the mean and confidence interval
cum_discr_forest_plot %>% 
  mutate(name = 23:1) %>% 
  ggplot()+
  geom_errorbar(aes(x = Mean, y = name,
                    xmin = `95%CI_low`,
                    xmax = `95%CI_upp`))+
  geom_point(aes(x = Mean, y = name))+
  coord_cartesian(xlim = c(0.5, 1))+
  theme_classic()




