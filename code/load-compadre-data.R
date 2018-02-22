
setwd("~/Dropbox/research/spline-surface-models/")

library(feather)

source("./code/load-compadre-data.R")


gform.set = NULL
out <- compadre_data_load(SURV = FALSE, FEC = FALSE, gform.sub = gform.set)$out
write_feather(out, "./data/all_growth_data.feather")

gform.set = NULL
out <- compadre_data_load(SURV = FALSE, FEC = TRUE, gform.sub = gform.set)$out
sp_sub <- unique(out$species[which(out$y > 1000)])
out <- out[-which(!is.na(match(out$species, sp_sub))), ]
write_feather(out, "./data/all_fec_data.feather")

gform.set = NULL
out <- compadre_data_load(SURV = TRUE, FEC = FALSE, gform.sub = gform.set)$out
write_feather(out, "./data/all_surv_data.feather")

gform.set <- "Tree"
out <- compadre_data_load(SURV = FALSE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/tree_growth_data.feather")

gform.set <- "Shrub"
out <- compadre_data_load(SURV = FALSE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/shrub_growth_data.feather")

gform.set <- "Herbaceous perennial"
out <- compadre_data_load(SURV = FALSE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/herb_growth_data.feather")

gform.set <- "Palm"
out <- compadre_data_load(SURV = FALSE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/palm_growth_data.feather")

gform.set <- "Tree"
out <- compadre_data_load(SURV = FALSE, FEC = TRUE, gform.sub = gform.set)$out
sp_sub <- unique(out$species[which(out$y > 1000)])
out <- out[-which(!is.na(match(out$species, sp_sub))), ]
write_feather(out, "./data/tree_fec_data.feather")

gform.set <- "Shrub"
out <- compadre_data_load(SURV = FALSE, FEC = TRUE, gform.sub = gform.set)$out
sp_sub <- unique(out$species[which(out$y > 1000)])
out <- out[-which(!is.na(match(out$species, sp_sub))), ]
write_feather(out, "./data/shrub_fec_data.feather")

gform.set <- "Herbaceous perennial"
out <- compadre_data_load(SURV = FALSE, FEC = TRUE, gform.sub = gform.set)$out
sp_sub <- unique(out$species[which(out$y > 1000)])
out <- out[-which(!is.na(match(out$species, sp_sub))), ]
write_feather(out, "./data/herb_fec_data.feather")

gform.set <- "Palm"
out <- compadre_data_load(SURV = FALSE, FEC = TRUE, gform.sub = gform.set)$out
write_feather(out, "./data/palm_fec_data.feather")

gform.set <- "Tree"
out <- compadre_data_load(SURV = TRUE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/tree_surv_data.feather")

gform.set <- "Shrub"
out <- compadre_data_load(SURV = TRUE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/shrub_surv_data.feather")

gform.set <- "Herbaceous perennial"
out <- compadre_data_load(SURV = TRUE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/herb_surv_data.feather")

gform.set <- "Palm"
out <- compadre_data_load(SURV = TRUE, FEC = FALSE, gform.sub = gform.set)$out

write_feather(out, "./data/palm_surv_data.feather")

