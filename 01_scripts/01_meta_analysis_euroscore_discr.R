# 01_meta_analysis_euroscore_discr.R
# run a Bayesian meta-analysis of AUC studies
# August 2024

# load libraries
library(dplyr)
library(nimble)
library(metamisc)
library(ggplot2)
library(boot)
library(renv)
library(viridis)

# renv things to do and check before commits
# renv::install('package')
# renv::status()
# renv::snapshot()
# renv::restore()

#### part 1: run the standard meta-analysis as a comparison ####
# bring in euroscore data set from metamisc
data(EuroSCORE)

# Meta-analysis of the c-statistic (random effects) through metamisc valmeta
fit_discr <- valmeta(cstat = c.index, cstat.se = se.c.index, cstat.cilb = c.index.95CIl,
                     cstat.ciub = c.index.95CIu, cstat.cilv = 0.95, N = n, O = n.events,
                     slab = Study, data = EuroSCORE)

# create a forest plot of the meta-analysis
plot(fit_discr)


# get summary
fit_discr

# check variance against sample size
plot(log2(EuroSCORE$n), fit_discr$data$theta.se^2)
median(EuroSCORE$n) # get the median study size

# save the output
saveRDS(fit_discr, file = "02_data/01_euro_valmeta_discr.RDA")



#### part 2a: run in nimble for the discrimination ####
## code; copied from metamisc, bugsmodels.r
code_discr <- nimbleCode({
 ## Likelihood
 for (i in 1:N){ # loop through studies
  theta[i] ~ dnorm(alpha[i], wsprec[i])
  alpha[i] ~ dnorm(mu.tobs, bsprec) # study-specific estimate centered on overall mean
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

## data
constants_discr <- list(N = nrow(EuroSCORE))

data_discr <- list(theta.var = fit_discr$data$theta.se^2, # use results from valmeta | square the se to get the var
             theta = fit_discr$data$theta) # this is the logit c-stat

## initial values
inits_discr <- list(bsTau = 1,
              mu.tobs = 0)

# parameters to store
parms_discr = c('mu.obs', 'pred.obs', 'bsTau')

# models
model_discr <- nimbleModel(code_discr, 
                           data = data_discr, 
                           inits = inits_discr, 
                           constants = constants_discr)

# chain details
n.chains = 2
thin = 5
MCMC = 10000
seeds = c(1234,5678) # one per chain

# MCMC samples
mcmc_out_discr <- nimbleMCMC(model = model_discr,
                             inits = inits_discr,
                             monitors = parms_discr,
                             niter = MCMC*2*thin, # times 2 for burn-in 
                             thin = thin,
                             nchains = n.chains, 
                             nburnin = MCMC,
                             summary = TRUE, 
                             setSeed = seeds,
                             WAIC = FALSE)

# view the summary of all chains
nimble_meta_discr <- mcmc_out_discr$summary$all.chains

# compare with valmeta result
fit_discr

# view the summary of a single chain
hist(mcmc_out_discr$samples$chain1[,1])
hist(mcmc_out_discr$samples$chain1[,2])

hist(mcmc_out_discr$samples$chain2[,1])
hist(mcmc_out_discr$samples$chain2[,2])

# save the output as a RDA file
saveRDS(nimble_meta_discr, file = "02_data/02_euro_nimble_meta_discr.RDA")





#### part 2b: create a cumulative meta-analysis for the discrimination ####
# create a list of data lists in cumulative to make into multiple models
# create an empty list to store data
data_list <- list()

# create n data lists in cumulative to complete meta-analysis each time
for(i in 1:nrow(EuroSCORE)){
  
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
                               niter = MCMC*2*thin, # times 2 for burn-in 
                               thin = thin,
                               nchains = n.chains, 
                               nburnin = MCMC,
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

# store the prediction interval of the meta-analysis
# 23 as it is the 23rd meta-analysis
discr_pred_nit <- cum_discr_plot[[23]][[1]][3,1:5]


# join the estimates and prediction
cummeta_euro_discr <- rbind(cum_discr_forest_plot, discr_pred_nit)

# save the output
saveRDS(cummeta_euro_discr, file = "02_data/03_euro_nimble_cummeta_discr.RDA")


# plot a forest plot for the discrimination
# for each of the cumulative models plot the mean and confidence interval
# ensure remove the first study
cum_discr_forest_plot %>% 
  mutate(name = 23:1) %>% 
  filter(name < 23) %>% # removes the first study so it is k-1 studies
  ggplot()+
  geom_point(aes(x = Mean, y = name), shape = 15)+
  geom_linerange(aes(y = name,
                     xmin = `95%CI_low`,
                     xmax = `95%CI_upp`))+
  geom_linerange(data = discr_pred_nit,aes(y = 0,
                                           xmin = `95%CI_low`,
                                           xmax = `95%CI_upp`))+
  coord_cartesian(xlim = c(0.25, 1))+
  scale_y_continuous(
    breaks = 0:22,  # Include 0 (prediction interval) and studies 1-22
    labels = c("Prediction Interval", 23:2)
  ) +
  theme_classic()+
  labs(x = "AUC Summary Estimate",
       y = "Number of Studies (EuroSCORE)")+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")


ggsave(filename = "03_figures/cumulative_meta_discr_dash_zero.jpg",
       width = 6,
       height = 4)



#### Part 2c: create a recursive cumulative meta-analysis plot ####

# load in data from part 2b
readRDS(file = "02_data/03_euro_nimble_cummeta_discr.RDA")

# Take the data from the meta-analysis and calculate the recursive data
# Plot into a graph
cummeta_euro_discr %>% head(23) %>% # renove the prediction interval
  mutate(recursive = (Mean-lag(Mean)),
         step = 1:23) %>% 
  ggplot(aes(x = step, y = recursive))+
  geom_line(linewidth = 1)+
  theme_classic()+
  labs(x = "Number of Studies (EuroSCORE)",
       y = "Change in AUC Summary Estimate")+
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")

# save the graph
ggsave(filename = "03_figures/recursive_cumulative_meta_discr_dash_zero.jpg",
       width = 6,
       height = 4)










#### Part 5a: simulate a new study for the discrimination ####
# bootstrap 1000 resamples of the c-statistic and it's variance
# set a seed for reproducible results
set.seed(42)

# store the data of logit c-stat in a new data frame to bootstrap
discr_boot <- fit_discr$data$theta

# create function to get the mean of the samples
mean_function <- function(data, indices) {
  sample_data <- data[indices]  # Resample the data using the indices
  return(mean(sample_data))    # Compute the mean
}

# boot strap the samples 1000 times
discr_boot_out <- boot(data = discr_boot,
                          statistic = mean_function,
                          R = 1000)
# check the output
discr_boot_out

# view the output as a plot
plot(discr_boot_out)

# store the mean value of the c-stat
boot_study_discr <- discr_boot_out$t0


# save the output
saveRDS(object = boot_study_discr, file = "02_data/04_euro_bootstrap_discr.RDA")


### complete now for the variance
# store the data of variance in a new data frame to bootstrap
discr_var_boot <- fit_discr$data$theta.se^2

# create function to get the mean of the samples
mean_function <- function(data, indices) {
  sample_data <- data[indices]  # Resample the data using the indices
  return(mean(sample_data))    # Compute the mean
}

# boot strap the samples 1000 times
discr_var_boot_out <- boot(data = discr_var_boot,
                          statistic = mean_function,
                          R = 1000)
# check the output
discr_var_boot_out

# view the output as a plot
plot(discr_var_boot_out)

# store the mean value of the variance
boot_study_var_discr <- discr_var_boot_out$t0

# save the output
saveRDS(object = boot_study_var_discr, file = "02_data/05_euro_bootstrap_var_discr.RDA")





#### part 6: completing a new meta-analysis with the bootstrapped simulated study for discrimination ####
# update the model and the data
# re-use the inits, monitors, niter, nchains, nburnin from previous

# add the simulated study to the data for the c-stat and varaince
data_discr_boot <- within(data_discr, {
  theta <- append(theta, boot_study_discr)
  theta.var <- append(theta.var, boot_study_var_discr)
})

## data
constants_discr_boot <- list(N = length(data_discr_boot[[1]])) # changed to add the new simulated study

## initial values
inits_discr_boot <- list(bsTau = 1,
                    mu.tobs = 0)

# parameters to store
parms_discr_boot = c('mu.obs', 'pred.obs', 'bsTau')

# models
model_discr_boot <- nimbleModel(code_discr, 
                           data = data_discr_boot, 
                           inits = inits_discr_boot, 
                           constants = constants_discr_boot)

# chain details
n.chains = 2
thin = 5
MCMC = 10000
seeds = c(1234,5678) # one per chain

# MCMC samples
mcmc_out_discr_boot <- nimbleMCMC(model = model_discr_boot,
                                  inits = inits_discr_boot,
                                  monitors = parms_discr_boot,
                                  niter = MCMC*2*thin, # times 2 for burn-in 
                                  thin = thin,
                                  nchains = n.chains, 
                                  nburnin = MCMC,
                                  summary = TRUE, 
                                  setSeed = seeds,
                                  WAIC = FALSE)

# view the summary of all chains
mcmc_out_discr_boot$summary$all.chains

# compare with the model without the simulated study
mcmc_out_discr$summary$all.chains


# save the summary
meta_bootstrap_discr <- mcmc_out_discr_boot$summary$all.chains

# save the file
saveRDS(object = meta_bootstrap_discr, file = "02_data/06_euro_meta_bootstrap_discr.RDA")



#### part 7a: simulate new data where the discrimination is 0.01 better ####
### while the sample size is at 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600 and 51200

# take the AUC and increment 0.01 improvements till 1.0
# create an empty list
discr_imp_sim <- list()

# loop auc improvements by 0.01 till 1.0
for (i in 1:((1 - round(mcmc_out_discr$summary$all.chains[2], digits = 2)) * 100)){
  
   tmp <- round(mcmc_out_discr$summary$all.chains[2], digits = 2) + (i * 0.01)
   
   discr_imp_sim[i] <- tmp
}

# simulate the sample sizes
discr_imp_sim_samp <- list(100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200)

# join the c-stat an the standard error together
# also where the AUC value is 1.0 change it to 0.99999 as 1.0 goes to infinity
discr_combined <- data.frame(expand.grid(auc = discr_imp_sim, sample = discr_imp_sim_samp)) %>%
  mutate(auc = ifelse(auc == 1.0, 0.99, auc))


# simulate standard errors from the different samples sizes
# formula from Newcombe 2005 to estimate standard error of the c-statistic

# average observed events
euro_observed <- mean(EuroSCORE$Po)
euro_non_observed <- 1-euro_observed

sim_data <- discr_combined %>% mutate(
  cstat = as.numeric(auc),
  theta = logit(cstat),
  sample = as.numeric(sample),
  n = sample*euro_observed, # using the average observed events from all the studies in the meta-analysis
  m = sample*euro_non_observed, # using the average non-observed events from all the studies in the meta-analysis
  star = (0.5*sample)-1,
  logit.se.c = sqrt(((1+star)*((1-cstat)/(2-cstat))+((star*cstat)/(1+cstat)))/(n*m*cstat*(1-cstat)))
) %>% select(
  cstat,
  sample,
  m,
  n,
  logit.se.c,
  theta
)

# convert the data frame to a list 
sim_data_list <- as.list(sim_data %>% select(theta, logit.se.c))

# check the sample size to variance
plot(log2(sim_data$sample), sim_data$logit.se.c)
plot(log2(sim_data$sample), sim_data$logit.se.c^2)


# Create an empty list to store all the study data
all_study_data <- list()

# Loop over the studies
for (i in 1:length(sim_data_list[[1]])) {
  
  # Create a list for the current study's data
  new_study_data <- list(
    theta.var = c(fit_discr$data$theta.se^2, sim_data_list$logit.se.c[i]^2), # variance
    theta = c(fit_discr$data$theta, sim_data_list$theta[i]) # logit c-stat
  )
  
  # Append the current study's list to the all_study_data list
  all_study_data[[i]] <- new_study_data
}

# data
sim_constants_discr <- list(N = nrow(EuroSCORE)+1) # add one for the new simulated study

# initial values
inits_discr <- list(bsTau = 1,
                    mu.tobs = 0)

# parameters to store
parms_discr = c('mu.obs', 'pred.obs', 'bsTau')

# chain details
n.chains = 2
thin = 5
MCMC = 10000
seeds = c(1234,5678) # one per chain


# create n models to then be run through NimbleMCMC
# create an empty list to store data
sim_study_model <- list()

# create n lists of models to complete meta-analysis each time with the new simulate study
for(i in 1:length(sim_data_list[[1]])){
  
  sim_model_discr <- nimbleModel(code_discr, 
                           data = all_study_data[[i]], 
                           inits = inits_discr, 
                           constants = sim_constants_discr)
  name <- paste('item', i, sep='')
  tmp <- list(model_1 = sim_model_discr)
  sim_study_model[[name]] <- tmp
}



# run the meta-analysis for each model in NimbleMCMC
# create an empty list to store outputs
sim_meta_discr <- list()

# run each cumulative model through MCMC
for(i in 1:length(sim_study_model)){
  
  mcmc_out_discr <- nimbleMCMC(model = sim_study_model[[i]][[1]],
                               inits = inits_discr,
                               monitors = parms_discr,
                               niter = MCMC*2*thin, # times 2 for burn-in 
                               thin = thin,
                               nchains = n.chains, 
                               nburnin = MCMC,
                               summary = TRUE, 
                               setSeed = seeds,
                               WAIC = FALSE)
  
  name <- paste('meta', i, sep='')
  tmp <- list(model = mcmc_out_discr)
  sim_meta_discr[[name]] <- tmp
}


# check that the meta-analysis was completed for each new simulated study
for (i in 1:length(sim_meta_discr)){
  print(sim_meta_discr[[i]]$model$summary$all.chains)
}

# take the output of all chains from the model
# create an empty list to store the outputs
sim_discr_plot <- list()

# take the outputs of each model with the sample size and increment improvement
for (i in 1:length(sim_meta_discr)){

  mu <- sim_meta_discr[[i]]$model$summary$all.chains
  mu_df <- as.data.frame(mu)
  name <- paste('study', i, sep='')
  tmp <- list(summary = mu_df)
  sim_discr_plot[[name]] <- tmp
}


# from all the data take mu and confidence intervals to then be plotted
# create an empty list to store data
sim_discr_square_plot <- data.frame()

# take mu and confidence intervals
for (i in 1:length(sim_discr_plot)){

  df <- as.data.frame(sim_discr_plot[[i]][[1]][2,1:5])
  sim_discr_square_plot <- rbind(sim_discr_square_plot, df)
}


# join this data with the sample size and incremental improvement data frame
sim_discr_square_plot_final <- cbind(sim_discr_square_plot, sim_data)

# save the data
saveRDS(sim_discr_square_plot_final, file = "02_data/07_eruo_sim_meta_discr.RDA")

# plot square plot of the sample size, improvement and overall 
# for each of the cumulative models plot the mean and confidence interval

# a heat mat tile plot
sim_discr_square_plot_final %>%
  ggplot()+
    geom_tile(aes(x = as.factor(sample), y = cstat, fill = Mean))+
  scale_x_discrete()+
  scale_fill_viridis(option="viridis")+
  labs(x = "Sample Size",
       y = "Simulated Study C-statistic",
       fill = "Meta-Analysis Summary Estiamte")+
  theme_classic()

ggsave(filename = "03_figures/sim_new_study_discr_euroscore_tile.jpg",
       width = 8,
       height = 4)


# a line plot
sim_discr_square_plot_final %>% 
  ggplot()+
  geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  labs(x = "Simulated Study AUC Value",
       y = "Meta-Analysis Summary Estiamte (EuroSCORE)",
       colour = "Sample Size")+
  scale_color_viridis(discrete=TRUE)+
  scale_x_continuous(breaks = c(0.8, 0.85, 0.9, 0.95, 0.99))+
  scale_y_continuous(breaks = c(0.79, 0.795, 0.80, 0.805, 0.81), limits = c(0.788, 0.81))


ggsave(filename = "03_figures/sim_new_study_discr_euroscore_line.jpg",
       width = 6,
       height = 4)











#### Part 4 create a funnel plot for the discrimination and calibration ####
# use the metamisc package function fat() to do so

# create the objects needed for the fat function
b_discr <- EuroSCORE$c.index
b_se_discr <- EuroSCORE$se.c.index
n_total_discr <- EuroSCORE$n
d_total_discr <- EuroSCORE$n.events

# complete the meta-regression using fat for the discrimination
funnel_discr <- EuroSCORE %>% fat(b = b_discr,
                                  b.se = b_se_discr,
                                  n.total = n_total_discr,
                                  d.total = n_total_discr,
                                  method = "D-FIV")

# plot the funnel plot for the discrimination
plot(x = funnel_discr,
     ref = 0.5,
     xlab = "C-Statistic")



# create the objects needed for the fat function
b_cal <- EuroSCORE$Po/EuroSCORE$Pe
b_se_cal <- EuroSCORE$SD.Pe/sqrt(EuroSCORE$n)
n_total_cal <- EuroSCORE$n
d_total_cal <- EuroSCORE$n.events

# complete the meta-regression using fat for the calibration
funnel_cal <- EuroSCORE %>% fat(b = b_cal,
                                b.se = b_se_cal,
                                n.total = n_total_cal,
                                d.total = d_total_cal,
                                method = "M-FIV")

# plot the funnel plot for the calibration
plot(x = funnel_cal,
     ref = 1,
     xlab = "Calibration (O:E)")





