# R code to 
#
# requires:
#   - load-data-helpers.R
#   - COMPADRE database
#   - trait data in 

# optional: set working directory
# setwd("PATH/TO/DIR")

# load libraries
library(feather)

# source function to load all data
source("./code/load-data-helpers.R")

# load COMPADRE database (check version and rename)
compadre_data <- get(load("./data/COMPADRE_v.X.X.X.RData"))

# load trait data
trait_data <- get(load("./data/imputed_trait_data.R"))

# load growth data for all growth forms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = NULL)
write_feather(out, "./data/all_growth_data.feather")

# load fecundity data for all growth forms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = NULL)
write_feather(out, "./data/all_fec_data.feather")

# load survival data for all growth forms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = NULL)
write_feather(out, "./data/all_surv_data.feather")

# load growth data for trees, shrubs, herbs and palms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Tree")
write_feather(out, "./data/tree_growth_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Shrub")
write_feather(out, "./data/shrub_growth_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Herbaceous perennial")
write_feather(out, "./data/herb_growth_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Palm")
write_feather(out, "./data/palm_growth_data.feather")

# load fecundity data for trees, shrubs, herbs and palms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Tree")
write_feather(out, "./data/tree_fec_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Shrub")
write_feather(out, "./data/shrub_fec_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Herbaceous perennial")
write_feather(out, "./data/herb_fec_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Palm")
write_feather(out, "./data/palm_fec_data.feather")

# load survival data for trees, shrubs, herbs and palms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Tree")
write_feather(out, "./data/tree_surv_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Shrub")
write_feather(out, "./data/shrub_surv_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Herbaceous perennial")
write_feather(out, "./data/herb_surv_data.feather")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Palm")
write_feather(out, "./data/palm_surv_data.feather")
