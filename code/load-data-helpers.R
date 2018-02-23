# helper function to tidy compadre data and trait data
compadre_data_load <- function(compadre_data = NULL,
                               trait_data = NULL,
                               vital = "growth",
                               gform = NULL) {
  
  # check all data available
  if (is.null(compadre_data))
    stop("COMPADRE data must be provided", call. = FALSE)
  if (is.null(trait_data))
    stop("trait data must be provided", call. = FALSE)
  
  # check vital rate
  if (!(vital %in% c("growth", "fec", "surv")))
    stop("vital must be one of growth, fec or surv", call. = FALSE)
  
  # subset compadre data to size-structured matrices with > 7 classes
  cp.sub <- which((compadre_data$metadata$MatrixCriteriaSize != "No") &
                    (compadre_data$metadata$MatrixDimension > 7) &
                    !is.na(sapply(compadre_data$mat, function(x) x[[2]][1, 1])))
  
  # pull out species names
  sp.names <- unique(compadre_data$metadata$SpeciesAccepted[cp.sub])

  # check trait data (replace all this with trait_data)
  all.traits <- trait_data

  # pull out data on growth form for each species
  gform.vec <- compadre$metadata$GrowthType[cp.sub][match(sp.names,
                                                          compadre$metadata$SpeciesAccepted[cp.sub])]
  
  # pull out population matrices for each vital rate
  cp.mats <- compadre$mat[cp.sub]
  cp.surv <- lapply(cp.mats, function(x) apply(x[[2]], 2, sum))
  cp.grow <- lapply(cp.mats, function(x) x[[2]])
  cp.fec <- lapply(cp.mats, function(x) apply(x[[3]] + x[[4]], 2, sum, na.rm = TRUE))
  
  # pull out metadata for subsetted species
  comp.meta <- compadre$metadata[cp.sub, ]

  # tidy survival data if needed (remove matrices with all zeros and truncate values > 1)  
  if (vital == "surv") {
    cp_rm <- which(sapply(cp.surv, function(x) max(x) == 0))
    cp.surv <- cp.surv[-cp_rm]
    mats.out <- cp.mats[-cp_rm]
    n.bins <- sapply(cp.surv, length)
    y <- matrix(NA, nrow=length(cp.surv), ncol=max(n.bins))
    for (i in seq(along=cp.surv)) {
      y[i, seq_len(n.bins[i])] <- ifelse(cp.surv[[i]] > 1, 1, cp.surv[[i]])
    }
    n.obs <- length(cp.surv)
  }
  
  # tidy fecundity data (remove matrices with all zeros and round to nearest whole number)
  if (vital == "fec") {
    cp_rm <- which(sapply(cp.fec, function(x) max(x) == 0))
    cp.fec <- cp.fec[-cp_rm]
    mats.out <- cp.mats[-cp_rm]
    n.bins <- sapply(cp.fec, length)
    y <- matrix(NA, nrow = length(cp.fec), ncol = max(n.bins))
    for (i in seq_along(cp.fec)) {
      y[i, seq_len(n.bins[i])] <- round(cp.fec[[i]])
    }
    n.obs <- length(cp.fec)
  }
  
  # tidy growth data (remove matrices with all zeros)
  if (vital == "growth") {
    cp_rm <- which(sapply(cp.fec, function(x) max(x) == 0))
    cp.grow <- cp.grow[-cp_rm]
    mats.out <- cp.mats[-cp_rm]
    n.bins <- sapply(cp.grow, ncol)
    y <- array(NA, dim = c(max(sapply(cp.grow, ncol)),
                           max(sapply(cp.grow, ncol)),
                           length(cp.grow)))
    for (i in seq_len(dim(y)[3])) {
      y[1:n.bins[i], 1:n.bins[i], i] <- cp.grow[[i]]
    }
  }
  
  # subset growth form, species, study and metadata objects to included species
  growth <- as.integer(compadre$metadata$GrowthType[cp.sub])
  growth <- growth[-cp_rm]
  growth <- ifelse(is.na(growth), max(growth, na.rm = TRUE) + 1, growth)
  growth <- as.integer(as.factor(growth))
  species <- as.integer(as.factor(compadre$metadata$SpeciesAccepted[cp.sub]))
  species <- species[-cp_rm]
  species <- as.integer(as.factor(species))
  study <- as.integer(as.factor(compadre$metadata$DOI.ISBN[cp.sub]))
  study <- study[-cp_rm]
  study <- as.integer(as.factor(study))
  comp.meta <- comp.meta[-cp_rm, ]
  
  # create trait matrix
  trait <- matrix(growth, nrow = 1)
  trait.store <- NULL
  for (i in seq_len(nrow(trait))) {
    trait.store <- rbind(trait.store, t(model.matrix(~factor(trait[i, ])))[-1, ])
  }
  x <- t(trait.store)
  colnames(x) <- NULL
  groups <- matrix(species, ncol = 1)
  
  # pull out a subset of traits for the included species
  sp.list <- compadre$metadata$SpeciesAccepted[cp.sub]
  sp.list <- sp.list[-cp_rm]
  all.traits <- all.traits[match(sp.list, rownames(all.traits)), ]

  # pull out rows that match the chosen growth form
  gform.list <- compadre$metadata$GrowthType[cp.sub][-cp_rm]  
  if (is.null(gform)) {
    gform.sub <- seq_len(nrow(all.traits))
  } else {
    gform.sub <- which(gform.list == gform)
  }
  
  # subset to the selected gform
  if (vital == "growth") {
    y.temp <- y[, , gform.sub]
    mats.out <- mats.out[gform.sub]
  } else {
    y.temp <- y[gform.sub, ]
    mats.out <- mats.out[gform.sub]
  }
  gform.list <- gform.list[gform.sub]
  x.temp <- x[gform.sub, ]
  groups.temp <- as.matrix(groups[gform.sub, ], ncol = 1)
  all.traits <- all.traits[gform.sub, ]
  comp.meta <- comp.meta[gform.sub, ]

  # remove any rows with missing trait data  
  final.sub <- which(apply(all.traits, 1, function(x) sum(is.na(x))) == 0)
  if (vital == "growth") {
    y <- y.temp[, , final.sub]
    mats.out <- mats.out[final.sub]
  } else {
    y <- y.temp[final.sub, ]
    mats.out <- mats.out[final.sub]
  }
  gform.list <- gform.list[final.sub]
  x <- traits.temp[final.sub, ]
  x <- sweep(x, 2, apply(x, 2, mean), "-")
  x <- sweep(x, 2, apply(x, 2, sd), "/")
  groups <- cbind(as.integer(as.factor(groups.temp[final.sub, ])),
                  as.integer(as.factor(growth[final.sub])))
  comp.meta <- comp.meta[final.sub, ]
  
  # create final output data frame
  out <- NULL
  if (vital == "growth") {
    nj <- apply(y[1, , ], 2, function(x) sum(!is.na(x)))
    for (i in seq_len(dim(y)[3])) {
      tmp <- data.frame(y = c(y[seq_len(nj[i]), seq_len(nj[i]), i]),
                        bins1 = rep(round(seq(0, 1, length = nj[i]), 2), times = nj[i]),
                        bins2 = rep(round(seq(0, 1, length = nj[i]), 2), each = nj[i]),
                        seed_mass = rep(x[i, 1], times = (nj[i] * nj[i])),
                        sla = rep(x[i, 2], times = (nj[i] * nj[i])),
                        ht = rep(x[i, 3], times = (nj[i] * nj[i])),
                        wood_den = rep(x[i, 4], times = (nj[i] * nj[i])),
                        gform = as.factor(rep(groups[i, 2], times = (nj[i] * nj[i]))),
                        species = as.factor(rep(groups[i, 1], times = (nj[i] * nj[i]))),
                        mat_id = as.factor(rep(i, times = (nj[i] * nj[i]))),
                        gform_name = rep(gform.list[i], times = (nj[i] * nj[i])))
      out <- rbind(out, tmp)
    }
  } else {
    nj <- apply(y, 1, function(x) sum(!is.na(x)))
    for (i in seq_len(nrow(y))) {
      tmp <- data.frame(y = c(y[i, seq_len(nj[i])]),
                        bins = round(seq(0, 1, length = nj[i]), 2),
                        seed_mass = rep(x[i, 1], times = nj[i]),
                        sla = rep(x[i, 2], times = nj[i]),
                        ht = rep(x[i, 3], times = nj[i]),
                        wood_den = rep(x[i, 4], times = nj[i]),
                        gform = as.factor(rep(groups[i, 2], times = nj[i])),
                        species = as.factor(rep(groups[i, 1], times = nj[i])),
                        mat_id = as.factor(rep(i, times = nj[i])),
                        gform_name = rep(gform.list[i], times = nj[i]))
      out <- rbind(out, tmp)
    }
  }
  
  # return outputs
  out <- list(out = out,
              mat = mats.out,
              traits = traits.temp[final.sub, ],
              comp.meta = comp.meta)
  out
  
}
