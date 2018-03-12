# R code to impute missing trait data using mice. This approach
#   creates multiple imputations for a multivariate data set with 
#   missing values. 

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
library(mice)

# load trait data (WARNING: could take > 20 min)
# source("./code/load-trait-data-bien.R")
all.traits <- read.csv("./data/bien_trait_data_not_imputed.csv",
                       row.names = 1)
sp.names <- read.csv("./data/bien_species_names.csv", row.names = 1)
gn.sp <- read.csv("./data/bien_genus_to_species.csv", row.names = 1)
fm.sp <- read.csv("./data/bien_family_to_species.csv", row.names = 1)

# set response variable for imputation
y <- all.traits

# pull out species with no observed trait data
rows_to_rm <- which(apply(y, 1, function(x) all(is.na(x))))
y <- y[-rows_to_rm, ]

# impute missing data using mice
y_imp <- mice(y)

# create complete output matrix
out_mat <- complete(y_imp)

# optional: save imputed traits
# write.csv(out_mat, file = "./data/bien_trait_data_imputed_mice.csv")
