# R code to 
#
# requires:
#   - load-data-helpers.R
#   - COMPADRE database
#   - trait data in 

# optional: set working directory
# setwd("PATH/TO/DIR")

# source function to load all data
source("./code/load-data-helpers.R")

# load COMPADRE database (check version and rename)
compadre_data <- get(load("./data/COMPADRE_v.4.0.1.RData"))

# load trait data
# trait_data <- read.csv("./data/bien_trait_data_imputed_mice.csv")
trait_data <- read.csv("./data/bien_trait_data_imputed.csv", row.names = 1)

# load growth data for all growth forms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = NULL)
save(out, file = "./data/all_growth_data.RData")

# load fecundity data for all growth forms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = NULL)
save(out, file = "./data/all_fec_data.RData")

# load survival data for all growth forms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = NULL)
save(out, file = "./data/all_surv_data.RData")

# load growth data for trees, shrubs, herbs and palms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Tree")
save(out, file = "./data/tree_growth_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Shrub")
save(out, file = "./data/shrub_growth_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Herbaceous perennial")
save(out, file = "./data/herb_growth_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "growth",
                          gform = "Palm")
save(out, file = "./data/palm_growth_data.RData")

# load fecundity data for trees, shrubs, herbs and palms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Tree")
save(out, file = "./data/tree_fec_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Shrub")
save(out, file = "./data/shrub_fec_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Herbaceous perennial")
save(out, file = "./data/herb_fec_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "fec",
                          gform = "Palm")
save(out, file = "./data/palm_fec_data.RData")

# load survival data for trees, shrubs, herbs and palms
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Tree")
save(out, file = "./data/tree_surv_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Shrub")
save(out, file = "./data/shrub_surv_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Herbaceous perennial")
save(out, file = "./data/herb_surv_data.RData")
out <- compadre_data_load(compadre_data = compadre_data,
                          trait_data = trait_data,
                          vital = "surv",
                          gform = "Palm")
save(out, file = "./data/palm_surv_data.RData")
