# R script to create supplementary plots for analysis

# load fitted models (summarised for plotting)

# create a list of all pre-plotting outputs
file_list <- dir("./outputs/pre-plot")

# load all outputs into a single vector
out <- vector("list", length = length(file_list))
for (i in seq_along(out))
  out[[i]] <- get(load(paste0("./outputs/pre-plot", file_list[i])))

# first plot main files
# Fig. 1: growth; all growth forms

# Fig. 2: trees; fec; SLA

# Fig. 3 herbs; growth; max_ht

# Fig. 4: herbs; survival; SLA

# Fig. 5: palms; growth; seed_mass

# work out names and then plot in order

# PLOT AS ONE FILE WITH CAPTIONS (use rmarkdown or latex?)

# Fig. S1: fecundity; all growth forms
# Fig. S2: survival; all growth forms

# S3: Tree growth_sd_mass
# S4: Tree growth_SLA
# S5: Tree growth_max_ht
# S6: Tree growth_woodden

# S7: Tree fec_sd_mass
# S8: Tree fec_max_ht
# S9: Tree fec_woodden

# S10: Tree surv_sd_mass
# S11: Tree surv_SLA
# S12: Tree surv_max_ht
# S13: Tree surv_woodden

# S14: Shrub growth_sd_mass
# S15: Shrub growth_SLA
# S16: Shrub growth_max_ht
# S17: Shrub growth_woodden

# S18: Shrub fec_sd_mass
# S19: Shrub fec_SLA
# S20: Shrub fec_max_ht
# S21: Shrub fec_woodden

# S22: Shrub surv_sd_mass
# S23: Shrub surv_SLA
# S24: Shrub surv_max_ht
# S25: Shrub surv_woodden

# S26: Herb growth_sd_mass
# S27: Herb growth_SLA

# S28: Herb fec_sd_mass
# S29: Herb fec_SLA
# S30 Herb fec_max_ht

# S31 Herb surv_sd_mass
# S32: Herb surv_max_ht

# S33: Palm growth_max_ht

# S34: Palm fec_sd_mass
# S35: Palm fec_max_ht

# S36: Palm surv_sd_mass
# S37: Palm surv_max_ht
