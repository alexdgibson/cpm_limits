# cumulative_meta_analysis_calibration.R
# run a cumulative Bayesian meta-analysis
# January 2025
# Alexander D Gibson

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

### part 2: complete a meta-analysis in metamisc ###
# Meta-analysis of the calibration (random effects)
# change the variable names to those in the data set
fit_cal <- valmeta(measure = "OE",
                   N = n,
                   O = n.events,
                   E = e.events,
                   slab = Study,
                   data = df)

# check summary of the metamisc meta-analysis
fit_cal

# plot the standard meta-analysis through metamisc
plot(fit_cal)

# check variance against sample size
plot(log2(df$n), fit_cal$data$theta.se)
median(df$n) # median study size



### part 2: create a cumulative meta-analysis for the calibration ###
# cumulative completed with NIMBLE in a bayesian framework

## code; copied from metamisc, bugsmodels.r
code_cal <- nimbleCode({
  for (i in 1:N){
    theta[i] ~ dnorm(alpha[i], wsprec[i])
    alpha[i] ~ dnorm(mu.tobs, bsprec)
    wsprec[i] <- 1/(theta.var[i])
  }
  # priors
  bsTau ~ T(dt(0, (1.5^2), df = 3), 0, 10) # truncated, see paper for prior set-up
  mu.tobs ~ dnorm(0, 1/1000) # overall mean
  pred.tobs ~ dnorm(mu.tobs, bsprec) # prediction interval
  # back-transform from natutal log
  log(mu.oe) <- mu.tobs
  log(pred.oe) <- pred.tobs
  bsprec <- 1/(bsTau*bsTau) # transform SD to precision
})

## set the constatnts
constants_cal <- list(N = nrow(df))

# set the data from the metamisc valmeta function
# it will calculate the O:E from other params when not directly available - through the valmeta function
data_cal <- list(theta.var = fit_cal$data$theta.se^2,
                 theta = fit_cal$data$theta)

## initial values
inits_cal <- list(bsTau = 1,
                  mu.tobs = 0)

# parameters to store
parms_cal = c('mu.oe', 'pred.oe', 'bsTau')


# specify the chain details from the protocol
niter = 100000 # one hundred thousand
thin = 5 # five
n.chains = 2 # two
seeds = c(1234,5678) # one seed per chain | specified in the protocol
n.burnin = 10000 # ten thousand



# create an empty list to store data
data_list_cal <- list()

# create n data lists in cumulative to complete meta-analysis each time
for(i in 1:nrow(df)){
  
  a <- as.numeric(data_cal[[1]][1:i])
  b <- as.numeric(data_cal[[2]][1:i])
  name <- paste('item', i, sep='')
  tmp <- list(theta.var = a, theta = b)
  data_list_cal[[name]] <- tmp
}

# create n models in cumulative to then be run through NimbleMCMC
# create an empty list to store data
cum_cal <- list()

# create n lists of models to complete meta-analysis each time
for(i in 1:length(data_list_cal)){
  
  model_cal <- nimbleModel(code_cal, 
                           data = data_list_cal[[i]], 
                           inits = inits_cal, 
                           constants = list(N = i))
  name <- paste('item', i, sep='')
  tmp <- list(model_1 = model_cal)
  cum_cal[[name]] <- tmp
}


# run the meta-analysis for each model in NimbleMCMC
# create an empty list to store outputs
cum_meta_cal <- list()

# run each cumulative model through MCMC
for(i in 1:length(cum_cal)){
  
  mcmc_out_cal <- nimbleMCMC(model = cum_cal[[i]][[1]],
                             inits = inits_cal,
                             monitors = parms_cal,
                             niter = niter, # times 2 for burn-in 
                             thin = thin,
                             nchains = n.chains, 
                             nburnin = n.burnin,
                             summary = TRUE, 
                             setSeed = seeds,
                             WAIC = FALSE)
  
  name <- paste('meta', i, sep='')
  tmp <- list(model = mcmc_out_cal)
  cum_meta_cal[[name]] <- tmp
}


# check that the meta-analysis was completed for each of the cumulative models
for (i in 1:length(cum_meta_cal)){
  print(cum_meta_cal[[i]]$model$summary$all.chains)
}

# take the output of all chains from the model
# create an empty list to store the outputs
cum_cal_plot <- list()

# take the outputs of each model
for (i in 1:length(cum_meta_cal)){
  
  mu <- cum_meta_cal[[i]]$model$summary$all.chains
  mu_df <- as.data.frame(mu)
  name <- paste('study', i, sep='')
  tmp <- list(summary = mu_df)
  cum_cal_plot[[name]] <- tmp
}


# from all the data take mu and confidence intervals to then be plotted
# create an empty list to store data
cum_cal_forest_plot <- data.frame()

# take mu and confidence intervals
for (i in 1:length(cum_cal_plot)){
  
  df <- as.data.frame(cum_cal_plot[[i]][[1]][2,1:5])
  cum_cal_forest_plot <- rbind(cum_cal_forest_plot, df)
}

# plot a forest plot for the discrimination
# for each of the cumulative models plot the mean and confidence interval
cum_cal_forest_plot %>% 
  mutate(name = 23:1) %>% 
  ggplot()+
  geom_errorbar(aes(x = Mean, y = name,
                    xmin = `95%CI_low`,
                    xmax = `95%CI_upp`))+
  geom_point(aes(x = Mean, y = name))+
  coord_cartesian(xlim = c(0, 10))+
  theme_classic()
