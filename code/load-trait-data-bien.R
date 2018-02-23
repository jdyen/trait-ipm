# R code to download trait data to match the COMPADRE database

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
library(BIEN)

# load COMPADRE data
load("./data/COMPADRE_v.4.0.1.RData")

# filter by species with eight or more size classes
cp.sub <- which((compadre$metadata$MatrixCriteriaSize != "No") &
                  (compadre$metadata$MatrixDimension > 7) &
                  !is.na(sapply(compadre$mat, function(x) x[[2]][1, 1])))

# compile list of species, genera, and family names
sp.names <- unique(compadre$metadata$SpeciesAccepted[cp.sub])
gn.names <- unique(compadre$metadata$Genus[cp.sub])
fm.names <- unique(compadre$metadata$Family[cp.sub])

# store information for citing/acknowledging data sources
sp.citations <- BIEN_metadata_citation(BIEN_occurrence_species(sp.names))
gn.citations <- BIEN_metadata_citation(BIEN_occurrence_genus(gn.names))
fm.citations <- BIEN_metadata_citation(BIEN_occurrence_family(as.character(fm.names)))

# store BEIN database version
bien_db <- BIEN_metadata_database_version()

# download BIEN trait data
bien.data.sp <- BIEN_trait_species(sp.names, all.taxonomy = TRUE)
bien.data.gn <- BIEN_trait_genus(gn.names, all.taxonomy = TRUE)
bien.data.fm <- BIEN_trait_family(as.character(fm.names), all.taxonomy = TRUE)

# prepare full list of available traits
trait.list <- unique(bien.data.sp$trait_name)

# create a matrix with of species by traits
sp.list <- sp.names
sp.traits <- matrix(NA,
                    nrow = length(sp.list),
                    ncol = length(trait.list))
for (i in seq(along = sp.list)) {
  sp_match <- which(bien.data.sp$scrubbed_species_binomial == sp.list[i])
  if (length(sp_match) > 0) {
    data_tmp <- bien.data.sp[sp_match, ]
    trait.vals <- tapply(suppressWarnings(as.numeric(data_tmp$trait_value)),
                         data_tmp$trait_name,
                         mean,
                         na.rm = TRUE)
    cols.to.fill <- match(names(trait.vals), trait.list)
    sp.traits[i, cols.to.fill] <- trait.vals
  }
}
rownames(sp.traits) <- sp.list
colnames(sp.traits) <- trait.list

# create a matrix with of genus by traits
gn.list <- gn.names
gn.traits <- matrix(NA,
                    nrow = length(gn.list),
                    ncol = length(trait.list))
for (i in seq(along = gn.list)) {
  gn_match <- which(bien.data.gn$scrubbed_genus == gn.list[i])
  if (length(gn_match) > 0) {
    data_tmp <- bien.data.gn[gn_match, ]
    trait.vals <- tapply(suppressWarnings(as.numeric(data_tmp$trait_value)),
                         data_tmp$trait_name,
                         mean,
                         na.rm = TRUE)
    cols.to.fill <- match(names(trait.vals), trait.list)
    gn.traits[i, cols.to.fill] <- trait.vals
  }
}
rownames(gn.traits) <- gn.list
colnames(gn.traits) <- trait.list

# create a matrix with of family by traits
fm.list <- fm.names
fm.traits <- matrix(NA,
                    nrow = length(fm.list),
                    ncol = length(trait.list))
for (i in seq(along = fm.list)) {
  fm_match <- which(bien.data.fm$scrubbed_family == fm.list[i])
  if (length(fm_match) > 0) {
    data_tmp <- bien.data.fm[fm_match, ]
    trait.vals <- tapply(suppressWarnings(as.numeric(data_tmp$trait_value)),
                         data_tmp$trait_name,
                         mean,
                         na.rm = TRUE)
    cols.to.fill <- match(names(trait.vals), trait.list)
    fm.traits[i, cols.to.fill] <- trait.vals
  }
}
rownames(fm.traits) <- fm.list
colnames(fm.traits) <- trait.list

# subset to the four traits we"re interested in and align the species, genus and family data
trait_sub <- c("seed mass", "whole plant leaf area per whole plant leaf dry mass", "whole plant height", "stem wood density")
all.traits <- sp.traits[, trait_sub]
sp.traits <- sp.traits[, trait_sub]
gn.traits <- gn.traits[, trait_sub]
fm.traits <- fm.traits[, trait_sub]
gn.sp <- compadre$metadata$Genus[match(sp.list, compadre$metadata$SpeciesAccepted)]
fm.sp <- compadre$metadata$Family[match(sp.list, compadre$metadata$SpeciesAccepted)]

# pull out genus and family data to match the rows in the full species matrix
gn.new <- gn.traits[match(gn.sp, rownames(gn.traits)), ]
fm.new <- fm.traits[match(fm.sp, rownames(fm.traits)), ]

# identify rows with all missing data in each trait matrix
na_sum <- apply(all.traits, 1, function(x) sum(is.na(x)))
na_sum_gn <- apply(gn.new, 1, function(x) sum(is.na(x)))
na_sum_fm <- apply(fm.new, 1, function(x) sum(is.na(x)))
sp_missing <- which(na_sum == 4)
gn_missing <- which(na_sum_gn == 4)
fm_missing <- which(na_sum_fm == 4)

# work out which rows are missing species data but not genus data
gn_sub <- which(!is.na(match(rownames(all.traits),
                             names(sp_missing[which(is.na(match(sp_missing,
                                                                gn_missing)))]))))
gn_replace <- gn.new[gn_sub, ]

# work out which rows are missing species data and genus data but not family data
fm_sub <- which(!is.na(match(rownames(all.traits),
                             names(sp_missing[which(is.na(match(sp_missing,
                                                                fm_missing)))]))))
fm_sub <- fm_sub[which(is.na(match(fm_sub, gn_sub)))]
fm_replace <- fm.new[fm_sub, ]

# prepare full lists of species data with genus and family data included
all.traits <- sp.traits
all.traits[gn_sub, ] <- gn_replace
rownames(all.traits)[gn_sub] <- paste0(rownames(all.traits)[gn_sub], "_gen_sub")
all.traits[fm_sub, ] <- fm_replace
rownames(all.traits)[fm_sub] <- paste0(rownames(all.traits)[fm_sub], "_fam_sub")

# return reduced COMPADRE data set
compadre_mat <- compadre$mat[cp.sub]
compadre_meta <- compadre$metadata[cp.sub, ]

# tidy workspace
rm(bien_db, bien.data.sp, bien.data.gn, bien.data.fm,
   cols.to.fill, compadre, compadre_data, cp.sub,
   data_tmp, fm_match, fm_missing, fm_replace, fm_sub,
   fm.list, fm.new, gn_match, gn_missing, gn_replace, gn_sub,
   gn.list, gn.new, i, na_sum, na_sum_fm, na_sum_gn,
   sp_match, sp_missing, sp.list, trait_sub, trait.list,
   trait.vals)
