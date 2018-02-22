# function to load data sets from COMPADRE data base with full trait data
compadre_data_load <- function(param = "growth", gform = "tree",
                               trait_data, compadre_mat,
                               compadre_meta, sp.list) {
  
  comp.meta <- compadre_meta
  gform.vec <- compadre_meta$GrowthType[match(sp.list, compadre_meta$SpeciesAccepted)]
  
  tree.trait <- trait_data[which(gform.vec == "Tree"), ]
  herb.trait <- trait_data[which(gform.vec == "Herbaceous perennial"), ]
  palm.trait <- trait_data[which(gform.vec == "Palm"), ]
  shrub.trait <- trait_data[which(gform.vec == "Shrub"), ]
  
  cp.mats <- compadre_mat
  cp.surv <- lapply(cp.mats, function(x) apply(x[[2]], 2, sum))
  cp.grow <- lapply(cp.mats, function(x) x[[2]])
  cp.fec <- lapply(cp.mats, function(x) apply(x[[3]] + x[[4]], 2, sum,
                                              na.rm = TRUE))
  
  if (param == 'surv') {
    cp.sub2 <- which(sapply(cp.surv, function(x) max(x) == 0))
    cp.surv <- cp.surv[-cp.sub2]
    mats.out <- cp.mats[-cp.sub2]
    n.bins <- sapply(cp.surv, length)
    y <- matrix(NA, nrow=length(cp.surv), ncol=max(n.bins))
    for (i in seq(along=cp.surv)) {
      y[i, 1:n.bins[i]] <- ifelse(cp.surv[[i]] > 1, 1, cp.surv[[i]])
    }
    n.obs <- length(cp.surv)
  } else {
    if (param == 'fec') {
      cp.sub2 <- which(sapply(cp.fec, function(x) max(x) == 0))
      cp.fec <- cp.fec[-cp.sub2]
      mats.out <- cp.mats[-cp.sub2]
      n.bins <- sapply(cp.fec, length)
      y <- matrix(NA, nrow=length(cp.fec), ncol=max(n.bins))
      for (i in seq(along=cp.fec)) {
        y[i, 1:n.bins[i]] <- round(cp.fec[[i]])
      }
      n.obs <- length(cp.fec)
    } else {
      if (param != "growth") {
        warning(paste0(param, " is not a known demographic parameter; param = 'growth used by default"),
                call. = FALSE)
      }
      cp.sub2 <- which(sapply(cp.fec, function(x) max(x) == 0))
      cp.grow <- cp.grow[-cp.sub2]
      mats.out <- cp.mats[-cp.sub2]
      n.bins <- sapply(cp.grow, ncol)
      y <- array(NA, dim = c(max(sapply(cp.grow, ncol)),
                             max(sapply(cp.grow, ncol)),
                             length(cp.grow)))
      for (i in 1:dim(y)[3]) {
        y[1:n.bins[i], 1:n.bins[i], i] <- cp.grow[[i]]
      }
    }
  }
  
  GFORM <- as.integer(as.factor(compadre_meta$GrowthType)[-cp.sub2])
  GFORM <- ifelse(is.na(GFORM), max(GFORM, na.rm = TRUE) + 1, GFORM)
  SPECIES <- as.integer(as.factor(compadre_meta$SpeciesAccepted[-cp.sub2]))
  STUDY <- as.integer(as.factor(compadre_meta$DOI.ISBN[-cp.sub2]))
  comp.meta <- comp.meta[-cp.sub2, ]
  
  trait <- matrix(GFORM, nrow = 1)
  trait.store <- NULL
  for (i in 1:nrow(trait)) {
    trait.store <- rbind(trait.store, t(model.matrix( ~ factor(trait[i, ])))[-1, ])
  }
  x <- t(trait.store)
  colnames(x) <- NULL
  groups <- matrix(SPECIES, ncol = 1)
  
  sp.list <- sp.list[-cp.sub2]
  traits3 <- all.traits[match(sp.list, rownames(all.traits)), ]
  
  gform.list <- compadre_meta$GrowthType[-cp.sub2]  
  if (gform == "all") {
    tree.sub <- 1:nrow(traits3)
  } else {
    if (is.na(match(gform, c("tree", "shrub", "herb", "palm")))) {
      warning(paste0("'", gform, "'", " is not one of 'tree, 'shrub', 'herb' or 'palm'; gform = 'tree' by default"),
              call. = FALSE)
      gform <- "tree"
    }
    gform_set <- switch(gform,
                        "tree" = "Tree",
                        "shrub" = "Shrub",
                        "herb" = "Herbaceous perennial",
                        "palm" = "Palm")
    tree.sub <- which(gform.list == gform)
  }
  
  if (param == "growth") {
    y.temp <- y[, , tree.sub]
    mats.out <- mats.out[tree.sub]
  } else {
    y.temp <- y[tree.sub, ]
    mats.out <- mats.out[tree.sub]
  }
  gform.list <- gform.list[tree.sub]
  x.temp <- x[tree.sub, ]
  groups.temp <- as.matrix(groups[tree.sub, ], ncol = 1)
  traits.temp <- traits3[tree.sub, ]
  comp.meta <- comp.meta[tree.sub, ]
  
  final.sub <- which(apply(traits.temp, 1, function(x) sum(is.na(x))) == 0)
  
  if (param == 'growth') {
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
                  as.integer(as.factor(GFORM[final.sub])))
  comp.meta <- comp.meta[final.sub, ]
  
  out <- NULL
  if (param == "growth") {
    nj <- apply(y[1, , ], 2, function(x) sum(!is.na(x)))
    for (i in 1:dim(y)[3]) {
      tmp <- data.frame(y = c(y[1:nj[i], 1:nj[i], i]),
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
    for (i in 1:nrow(y)) {
      tmp <- data.frame(y = c(y[i, 1:nj[i]]),
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
  
  return(list(out = out,
              mat = mats.out,
              traits = traits.temp[final.sub, ],
              comp.meta = comp.meta))
  
}

