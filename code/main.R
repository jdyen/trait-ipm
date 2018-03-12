# R code to analyse size-structured demographic vital rates as functions
#   of traits. This code uses Gaussian Processes (GP) to model vital rates
#   as continuous funtions of size. The model uses a sparse GP approximation.
#
# Requires:
#   -pre-compiled .feather data files for each growth form
#   -helper functions (helpers.R)
#   -reticulate R package
#   -python3 installed locally with GPflow, numpy, feather and pandas libraries

# optional: clear workspace
# rm(list = ls())

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
library(reticulate)

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
  
  # inverse logistic transform growth data to convert (0, 1) data to (-Inf, Inf)
  data_set$Y <- qlogis(data_set$Y)
  
  # fit model using fit_gp_mod() helper function
  mod_grow_tmp <- fit_gp_mod(data_set,
                             n_iter = num_iter,
                             n_burn = num_burn,
                             eps = eps_set[1],
                             lmin = 1, lmax = 50,
                             max_it_optim = 1000)
  
  # cross validated fitted growth model
  mod_grow_tmp_cv <- gp_cv(mod_grow_tmp,
                           n_cv = num_cv,
                           n_iter = num_iter_cv,
                           n_burn = num_burn_cv,
                           eps = eps_set[1],
                           lmin = 1, lmax = 50,
                           max_it_optim = 1000,
                           par_run = FALSE)
  
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
  
  # cross validated fitted fecundity model
  mod_fec_tmp_cv <- gp_cv(mod_fec_tmp,
                          n_cv = num_cv,
                          n_iter = num_iter_cv,
                          n_burn = num_burn_cv,
                          eps = eps_set[2],
                          lmin = 1, lmax = 50,
                          max_it_optim = 1000,
                          par_run = FALSE)
  
  # load survival data for growth form i
  data_set <- load_data(growth_form = gform_opt[i],
                        param_type = "surv",
                        n_induce = num_induce[i, 3])
  
  # remove 1s and 0s
  data_set$Y <- ifelse(data_set$Y >= 1, 0.9999, data_set$Y)
  data_set$Y <- ifelse(data_set$Y <= 0, 0.0001, data_set$Y)
  
  # inverse logistic transform survival data to convert (0, 1) data to (-Inf, Inf)
  data_set$Y <- qlogis(data_set$Y)
  
  # fit model
  mod_surv_tmp <- fit_gp_mod(data_set,
                             n_iter = num_iter,
                             n_burn = num_burn,
                             eps = eps_set[3],
                             lmin = 1, lmax = 50,
                             max_it_optim = 1000)

  # cross validated fitted survival model
  mod_surv_tmp_cv <- gp_cv(mod_surv_tmp,
                           n_cv = num_cv,
                           n_iter = num_iter_cv,
                           n_burn = num_burn_cv,
                           eps = eps_set[3],
                           lmin = 1, lmax = 50,
                           max_it_optim = 1000,
                           par_run = FALSE)
  
  # collate full model and cross-validated model in a single list
  mod_grow <- list(mod = mod_grow_tmp, mod_cv = mod_grow_tmp_cv)
  mod_fec <- list(mod = mod_fec_tmp, mod_cv = mod_fec_tmp_cv)
  mod_surv <- list(mod = mod_surv_tmp, mod_cv = mod_surv_tmp_cv)
  
  # summary fitted models using gp_mod_summary() helper function
  out_grow <- gp_mod_summary(mod_grow)
  out_fec <- gp_mod_summary(mod_fec)
  out_surv <- gp_mod_summary(mod_surv)
  
  # save outputs to outputs directory
  save(out_grow, file = paste0("./outputs/fitted_growth_", gform_opt[i], ".R"))
  save(out_fec, file = paste0("./outputs/fitted_fec_", gform_opt[i], ".R"))
  save(out_surv, file = paste0("./outputs/fitted_surv_", gform_opt[i], ".R"))
  
}
