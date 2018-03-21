# R code to analyse size-structured demographic vital rates as functions
#   of traits. This code uses Gaussian Processes (GP) to model vital rates
#   as continuous funtions of size. The model uses a sparse GP approximation.
#
# Requires:
#   -pre-compiled .RData data files for each growth form
#   -helper functions (helpers.R)
#   -reticulate R package
#   -python3 installed locally with GPflow, numpy, feather and pandas libraries

# optional: clear workspace
# rm(list = ls())

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
library(reticulate)
library(parallel)
library(viridis)
library(ggplot2)

# required: set python version
use_python("/usr/local/bin/python3", required = TRUE)

# set up helper functions
source("./code/helpers.R")

# create vector of possible growth forms to loop over
gform_opt <- c("all", "tree", "shrub", "herb", "palm")

# define GP model settings
num_iter <- 5000         # total number of iterations
num_burn <- 2500         # number of burn-in iterations to discard
num_cv <- "loo"          # blocks by species, "loo" holds out one species per fold
num_iter_cv <- 5000      # number of iterations for cross validation
num_burn_cv <- 2500      # number of burn-in iterations to discard in cross validation
max_it_optim_set <- 1000 # number of iterations to optimise starting values

# number of inducing points for sparse GP model
num_induce <- rbind(c(70, 20, 10),
                    c(70, 20, 10),
                    c(30, 20, 10),
                    c(30, 100, 70),
                    c(30, 20, 10))
rownames(num_induce) <- gform_opt
colnames(num_induce) <- c("growth", "fecundity", "survival")

# set traits to plot for each growth form
gform_traits <- list(c("all"),
                     c("sla", "wood_den", "seed_mass", "ht"),
                     c("sla", "wood_den", "seed_mass", "ht"),
                     c("sla", "seed_mass", "ht"),
                     c("seed_mass", "ht"))

# set tolerances for each vital rate (growth, fecundity, survival)
eps_set <- c(0.002, 0.005, 0.005)

# loop through all growth forms and fit GP models
for (i in seq_along(gform_opt)) {
  
  # load growth data for growth form i
  data_set <- load_data(growth_form = gform_opt[i],
                        param_type = "growth",
                        n_induce = num_induce[i, 1])
  
  # remove 1s and 0s
  data_set$Y <- ifelse(data_set$Y >= 1, 0.9999, data_set$Y)
  data_set$Y <- ifelse(data_set$Y <= 0, 0.0001, data_set$Y)
  data_set$Y <- ifelse(is.na(data_set$Y), mean(data_set$Y, na.rm = TRUE), data_set$Y)
  
  # inverse logistic transform growth data to convert (0, 1) data to (-Inf, Inf)
  data_set$Y <- qlogis(data_set$Y)
  
  # fit model using fit_gp_mod() helper function
  mod_grow_tmp <- fit_gp_mod(data_set,
                             n_iter = num_iter,
                             n_burn = num_burn,
                             eps = eps_set[1],
                             lmin = 1, lmax = 50,
                             max_it_optim = 1000)
  
# load fecundity data for growth form i
  data_set <- load_data(growth_form = gform_opt[i],
                        param_type = "fec",
                        n_induce = num_induce[i, 2])
  
  # log-transform fecundity data
  data_set$Y <- log((data_set$Y / max(data_set$Y)) + 0.001)
  
  # fit model
  mod_fec_tmp <- fit_gp_mod(data_set,
                            n_iter = num_iter,
                            n_burn = num_burn,
                            eps = eps_set[2],
                            lmin = 1, lmax = 50,
                            max_it_optim = 1000)
  
  # load survival data for growth form i
  data_set <- load_data(growth_form = gform_opt[i],
                        param_type = "surv",
                        n_induce = num_induce[i, 3])
  
  # remove 1s and 0s
  data_set$Y <- ifelse(data_set$Y >= 1, 0.9999, data_set$Y)
  data_set$Y <- ifelse(data_set$Y <= 0, 0.0001, data_set$Y)
  data_set$Y <- ifelse(is.na(data_set$Y), mean(data_set$Y, na.rm = TRUE), data_set$Y)
  
  # inverse logistic transform survival data to convert (0, 1) data to (-Inf, Inf)
  data_set$Y <- qlogis(data_set$Y)
  
  # fit model
  mod_surv_tmp <- fit_gp_mod(data_set,
                             n_iter = num_iter,
                             n_burn = num_burn,
                             eps = eps_set[3],
                             lmin = 1, lmax = 50,
                             max_it_optim = 1000)

  # create plot outputs for fitted models
  for (j in seq_along(gform_traits[[i]])) {
    
    # prepare outputs to plot growth
    out_grow <- prepare_plot_outputs(mod_grow_tmp, trait_name = gform_traits[[i]][j])
    
    # prepare outputs to plot fecundity
    out_fec <- prepare_plot_outputs(mod_fec_tmp, trait_name = gform_traits[[i]][j])
    
    # prepare outputs to plot survival
    out_surv <- prepare_plot_outputs(mod_surv_tmp, trait_name = gform_traits[[i]][j])
    
    # save outputs to outputs directory
    save(out_grow, file = paste0("./outputs/pre-plot/growth_", gform_opt[i], "-", gform_traits[[i]][j], ".RData"))
    save(out_fec, file = paste0("./outputs/pre-plot/fec_", gform_opt[i], "-", gform_traits[[i]][j], ".RData"))
    save(out_surv, file = paste0("./outputs/pre-plot/surv_", gform_opt[i], "-", gform_traits[[i]][j], ".RData"))
    
  }
  
}
