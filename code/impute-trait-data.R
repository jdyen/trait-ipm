# R code to fit a Stan model that imputes missing trait data. This model
#   projects the entire species matrix into a lower dimensional space
#   and uses correlations among species to impute missing values.
# See Schrodt et al. (2015, Global Ecol. Biogeogr. 24: 1510-1521) for further details.

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
library(rstan)
library(BIEN)

# load trait data (WARNING: could take > 20 min)
# source("./code/load-trait-data-bien.R")

# set response variable for imputation
y <- all.traits

# pull out species with no observed trait data
rows_to_rm <- which(apply(y, 1, function(x) all(is.na(x))))

# create a data frame matching species to genera to families
sp_to_fam <- data.frame(sp = factor(sp.list),
                        gn = factor(gn.sp),
                        fm = fm.sp)

# remove species with no trait data
sp_to_fam <- sp_to_fam[-rows_to_rm, ]

# convert factors to integers to use at indices in Stan model
sp_to_fam$sp_int <- as.integer(as.factor(as.integer(sp_to_fam$sp)))
sp_to_fam$gn_int <- as.integer(as.factor(as.integer(sp_to_fam$gn)))
sp_to_fam$fm_int <- as.integer(as.factor(as.integer(sp_to_fam$fm)))

# convert response variable to a vector
y <- c(y[-rows_to_rm, ])
missing_ids <- which(is.na(y))
y <- y[-missing_ids]

# calculate indices for Stan model
nsp <- nrow(sp_to_fam)
n <- length(y)
ntrait <- ncol(all.traits)

# set number of latent variables
nlatent <- 6

# convert species, genus and family identities to vectors matching the elements of y
sp <- rep(sp_to_fam$sp_int, times = ntrait)
sp <- sp[-missing_ids]
gen <- rep(NA, times = length(unique(sp_to_fam$sp_int)))
gen[unique(sp_to_fam$sp_int)] <- sp_to_fam$gn_int[match(unique(sp_to_fam$sp_int), sp_to_fam$sp_int)]
gen[which(is.na(gen))] <- max(gen, na.rm = TRUE) + 1
trait <- rep(1:ntrait, each = nsp)
trait <- trait[-missing_ids]
fam <- rep(NA, times = length(unique(sp_to_fam$gn_int)))
fam[unique(sp_to_fam$gn_int)] <- sp_to_fam$fm_int[match(unique(sp_to_fam$gn_int), sp_to_fam$gn_int)]

# collate Stan data set
data_set <- list(n = n,
                 y = y,
                 nsp = nsp,
                 ntrait = ntrait,
                 nlatent = nlatent,
                 ngen = length(unique(gen)),
                 nfam = length(unique(fam)),
                 sp = sp,
                 trait = trait,
                 gen = gen,
                 fam = fam)

# compile Stan model
mod_compiled <- stan_model("./code/impute-trait-data.stan")

# set MCMC settings
n_chains <- 3
n_iter <- 2000

# fit Stan model
mod <- sampling(object = mod_compiled,
                data = data_set,
                chains = n_chains,
                iter = n_iter,
                init = "0",
                control = list(adapt_delta = 0.96, max_treedepth = 20),
                cores = n_chains)

# optional: save fitted model
# save(mod, file = "./outputs/imputed-trait-data.R")

# summarise fitted model and store mean values in a matrix
out <- summary(mod, pars = "trait_mat")$summary
out_mat <- matrix(out[, "mean"], ncol = ncol(all.traits), byrow = TRUE)
