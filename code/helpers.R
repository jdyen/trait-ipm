# import Python modules
GPflow <- import("GPflow")
np <- import("numpy")
feather <- import("feather")
pd <- import("pandas")

# define operators for GP kernels (overload R "+" and "*")
`+.GPflow.kernels.Kern` <- function (x, y) {
  GPflow$kernels$Add(list(x, y))
}
`*.GPflow.kernels.Kern` <- function (x, y) {
  GPflow$kernels$Prod(list(x, y))
}

# load data from pre-compiled data sets
load_data <- function(growth_form = "tree", param_type = "growth", n_induce = 10) {
  
  filename <- paste0("./data/", growth_form, "_", param_type, "_data.feather")
  df <- feather(filename)
  if (param_type == "growth") {
    if (growth_form == "tree" | growth_form == "shrub") {
      xx <- data.frame("bins1" = df$bins1,
                       "bins2" = df$bins2,
                       "seed_mass" = df$seed_mass,
                       "sla" = df$sla,
                       "ht" = df$ht,
                       "wood_den" = df$wood_den)
      z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    } else {
      if (growth_form == "herb") {
        xx <- data.frame("bins1" = df$bins1,
                         "bins2" = df$bins2,
                         "seed_mass" = df$seed_mass,
                         "sla" = df$sla,
                         "ht" = df$ht)
        z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      } else {
        if (growth_form == "palm") {
          xx <- data.frame("bins1" = df$bins1,
                           "bins2" = df$bins2,
                           "seed_mass" = df$seed_mass,
                           "ht" = df$ht)
          z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        } else {
          xx <- data.frame("bins1" = df$bins1,
                           "bins2" = df$bins2)
          z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        }
      }
    }
  } else {
    if (growth_form == "tree" | growth_form == "shrub") {
      xx <- data.frame("bins" = df$bins,
                       "seed_mass" = df$seed_mass,
                       "sla" = df$sla,
                       "ht" = df$ht,
                       "wood_den" = df$wood_den)
      z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    } else {
      if (growth_form == "herb") {
        xx <- data.frame("bins" = df$bins,
                         "seed_mass" = df$seed_mass,
                         "sla" = df$sla,
                         "ht" = df$ht)
        z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      } else {
        if (growth_form == "palm") {
          xx <- data.frame("bins" = df$bins,
                           "seed_mass" = df$seed_mass,
                           "ht" = df$ht)
          z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        } else {
          xx <- data.frame("bins" = df$bins)
          z_range <- apply(xx, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        }
      }
    }
  }
  xxsp <- as.numeric(as.character(df$species))
  xxmat <- as.numeric(as.character(df$mat_id))
  yy <- matrix(df$y, ncol = 1)
  if (growth_form == "all") {
    xxgform <- as.numeric(as.character(df$gform))
    xx <- data.frame(xx, "gform" = xxgform, "species" = xxsp, "matid" = xxmat)
  } else {
    xx <- data.frame(xx, "species" = xxsp, "matid" = xxmat)
  }
  if (growth_form == "all") {
    z_tmp <- c(c(max(xxgform, na.rm = TRUE) - min(xxgform, na.rm = TRUE)),
               c(max(xxsp, na.rm = TRUE) - min(xxsp, na.rm = TRUE)),
               c(max(xxmat, na.rm = TRUE) - min(xxmat, na.rm = TRUE)))
  } else {
    z_tmp <- c(c(max(xxsp, na.rm = TRUE) - min(xxsp, na.rm = TRUE)),
               c(max(xxmat, na.rm = TRUE) - min(xxmat, na.rm = TRUE)))
  }
  z_range <- c(z_range, z_tmp)
  zz <- sweep(matrix(runif((n_induce * ncol(xx))), ncol = ncol(xx)),
              2, z_range, "*")
  xx <- as.matrix(xx)
  
  return(list(Y = yy, X = xx, Z = zz, gform = growth_form, param = param_type))
  
}

# define IPM kernels as GPs
kernel <- function(n_dim = 6, gform = "tree") {
  
  kern <- GPflow$kernels$RBF(as.integer(n_dim), ARD = TRUE, active_dims = c(0:(n_dim - 1))) +
    GPflow$kernels$White(1, active_dims = matrix(n_dim)) +
    GPflow$kernels$White(1, active_dims = matrix(c(n_dim + 1)))
  
  return(kern = kern)
  
}

# fit a sparse GP model from a pre-loaded data set
fit_gp_mod <- function(data_set, n_iter = 1000, n_burn = 1000,
                       eps = 0.005, lmin = 1, lmax = 25,
                       max_it_optim = 1000) {
  
  yy <- data_set$Y
  xx <- data_set$X
  zz <- data_set$Z
  
  # set likelihood
  if (data_set$param == "fec") {
    gp_likl <- GPflow$likelihoods$Gaussian()
  } else {
    gp_likl <- GPflow$likelihoods$Gaussian()
  }
  
  # define mod
  mod <- GPflow$sgpmc$SGPMC(X = xx,
                            Y = yy,
                            kern = kernel(n_dim = (ncol(xx) - 2),
                                          gform = data_set$gform),
                            likelihood = gp_likl,
                            Z = zz)

  # set priors
  mod$kern$rbf$lengthscales$prior <- GPflow$priors$Gamma(2.0, 0.1)
  mod$kern$rbf$lengthscales$fixed <- FALSE
  mod$kern$rbf$variance$fixed <- TRUE
  mod$kern$white_1$variance$prior <- GPflow$priors$Gamma(0.1, 0.1)
  mod$kern$white_2$variance$prior <- GPflow$priors$Gamma(0.1, 0.1)
  mod$kern$white_1$variance$fixed <- FALSE
  mod$kern$white_2$variance$fixed <- FALSE

  # optimise for initial conditions
  mod$optimize(maxiter = max_it_optim)

  # full MCMC
  samples <- mod$sample(as.integer(n_iter),
                        epsilon = eps,
                        Lmin = as.integer(lmin),
                        Lmax = as.integer(lmax),
                        burn = as.integer(n_burn),
                        verbose = TRUE)

  # return fitted model
  return(list(mod = mod,
              samples = samples,
              param = data_set$param,
              gform = data_set$gform))
  
}

# extract fitted model summary (kernels, fitted values, deviance and r2)
extract_outputs <- function(mod_tmp) {
  
  m <- mod_tmp$mod
  samp <- mod_tmp$samples
  gform = mod_tmp$gform
  
  if (gform != "all") {
    kern_samp <- m$kern$rbf$lengthscales$get_samples_df(samp)$values
    kern_samp <- matrix(unlist(kern_samp), ncol = length(kern_samp), byrow = FALSE)
  } else {
    kern_samp <- m$kern$rbf$lengthscales$get_samples_df(samp)$values
    kern_samp <- matrix(unlist(kern_samp), ncol = length(kern_samp), byrow = FALSE)
  }
  
  kern_mean <- apply(kern_samp, 1, mean)
  kern_sd <- apply(kern_samp, 1, sd)
  
  y_est <- m$predict_y(m$X$value)[[1]]
  dev_val <- -2. * m$compute_log_likelihood()
  r2_val <- cor(y_est, m$Y$value, use = "complete") ** 2
  
  return(list(kern_mean = kern_mean,
              kern_sd = kern_sd,
              y_est = y_est,
              devi = dev_val,
              r2 = r2_val))
  
}

# internal cross-validation function
gp_cv_int <- function(i, m_full,
                      n_cv,
                      cv_size,
                      xsp, nsp,
                      max_it_optim,
                      n_iter, eps,
                      lmin, lmax,
                      n_burn) {
  
  print(paste0("Starting cross validation for species ", i))
  y_pred <- NULL
  y_real <- NULL
  if (i < n_cv) {
    sp_set <- ((i - 1) * cv_size + 1):(i * cv_size)
  } else {
    sp_set <- ((i - 1) * cv_size + 1):nsp
  }
  index_set <- which(!is.na(match(xsp, sp_set)))
  y_train <- as.matrix(m_full$mod$Y$value[-index_set, ])
  x_train <- as.matrix(m_full$mod$X$value[-index_set, ])
  z_train <- as.matrix(m_full$mod$Z$value)
  y_test <- m_full$mod$Y$value[index_set, ] 
  x_test <- m_full$mod$X$value[index_set, ]
  mod_cv <- GPflow$sgpmc$SGPMC(X = x_train, Y = y_train,
                               kern = kernel(n_dim = (ncol(x_train) - 2),
                                             gform = m_full$gform),
                               likelihood = m_full$mod$likelihood,
                               mean_function = GPflow$mean_functions$Constant(1),
                               Z = z_train)
  
  mod_cv$kern$rbf$lengthscales$prior <- GPflow$priors$Gamma(2.0, 0.1)
  mod_cv$kern$rbf$lengthscales$fixed <- FALSE
  mod_cv$kern$rbf$variance$fixed <- TRUE
  
  mod_cv$kern$white_1$variance$prior <- GPflow$priors$Gamma(0.1, 0.1)
  mod_cv$kern$white_2$variance$prior <- GPflow$priors$Gamma(0.1, 0.1)
  mod_cv$kern$white_1$variance$fixed <- FALSE
  mod_cv$kern$white_2$variance$fixed <- FALSE
  
  mod_cv$optimize(maxiter = max_it_optim)
  
  samples <- mod_cv$sample(as.integer(n_iter),
                           epsilon = eps,
                           Lmin = as.integer(lmin), Lmax = as.integer(lmax),
                           burn = as.integer(n_burn), verbose = TRUE)
  
  for (s in 1:nrow(samples)) {
    mod_cv$set_state(as.matrix(samples[s, ]))
    y_pred <- cbind(y_pred, mod_cv$predict_y(x_test)[[1]])
  }
  y_pred <- apply(y_pred, 1, mean)
  rm(samples, mod_cv)
  gc()
  return(list(y_pred = y_pred, y_real = y_test))
}

# cross validate fitted model with n_cv folds
gp_cv <- function(m_full,
                  n_cv = 10, n_iter = 2000, n_burn = 1000,
                  eps = 0.005, lmin = 1, lmax = 30,
                  max_it_optim = 1000,
                  par_run = FALSE) {
  xsp <- m_full$mod$X$value[, ncol(m_full$mod$X$value) - 1]
  nsp <- length(unique(xsp))
  if (n_cv == "loo") {
    n_cv <- nsp
  }
  cv_size <- ifelse(floor(nsp / n_cv) < 1, 1, floor(nsp / n_cv))
  cv_out <- mclapply(1:n_cv, gp_cv_int, m_full, n_cv, cv_size, xsp, nsp,
                     max_it_optim,
                     n_iter, eps, lmin, lmax, n_burn,
                     mc.cores = ifelse(par_run, parallel::detectCores(), 1))
  y_pred_all <- unlist(lapply(cv_out, function(x) x$y_pred))
  y_real_all <- unlist(lapply(cv_out, function(x) x$y_real))
  r2_cv <- cor(y_pred_all, y_real_all, use = "complete") ** 2
  return(list(y_pred = y_pred_all,
              y_real = y_real_all,
              r2cv = r2_cv))
}

# predict from fitted GP model (helper function to plot
#   fitted models and calculate eigenvalues)
predict_gp <- function(mod_full,
                       trait_name = "sla",
                       size_res = 100, seq_len = 12,
                       trait_min = -3., trait_max = 3.) {
  m <- mod_full$mod
  size_res2 <- size_res ** 2.
  if (mod_full$gform == "all") {
    if (trait_name != "all") {
      warning("model was fitted to all growth forms; predict function will default to trait_name = 'all'",
              call. = FALSE)
      trait_name <- "all"
    }
    if (seq_len != 4) {
      warning("predict function will set seq_len = 4 for model of all growth forms",
              call. = FALSE)
    }
    trait_seq <- c(1:4)
    col_id <- 2
    seq_len <- length(trait_seq)
  } else {
    trait_seq <- seq(trait_min, trait_max, length = seq_len)
    if (trait_name == "sla") {
      col_id <- 3
    } else {
      if (trait_name == "seed_mass") {
        col_id <- 2
      } else {
        if (trait_name == "ht") {
          col_id <- 4
          if (mod_full$gform == "palm") {
            col_id <- 3
          }
        } else {
          if (trait_name == "wood_den") {
            col_id <- 5 
          }
          else {
            warning(paste0("'", trait_name, "'", " is not a valid trait_name; trait_name = 'seed_mass' by default"),
                    call. = FALSE)
            col_id <- 2
          }
        }
      }
    }
  }
  if (mod_full$param == "growth") {
    zz_mean <- array(NA, dim = c(size_res, size_res, seq_len))
    zz_var <- array(NA, dim = c(size_res, size_res, seq_len))
  } else {
    zz_mean <- array(NA, dim = c(seq_len, size_res))
    zz_var <- array(NA, dim = c(seq_len, size_res))
  }
  for (i in c(1:seq_len)) {
    if (mod_full$param == "growth") {
      if ((mod_full$gform == "tree") | (mod_full$gform == "shrub")) {
        x_new <- cbind(rep(seq(0, 1, length = size_res), times = size_res),
                       rep(seq(0, 1, length = size_res), each = size_res),
                       rep(0, times = size_res2), 
                       rep(0, times = size_res2), 
                       rep(0, times = size_res2), 
                       rep(0, times = size_res2), 
                       rep(0, times = size_res2), 
                       rep(1, times = size_res2), 
                       rep(1, times = size_res2))
      } else {
        if (mod_full$gform == "herb") {
          x_new <- cbind(rep(seq(0, 1, length = size_res), times = size_res),
                         rep(seq(0, 1, length = size_res), each = size_res),
                         rep(0, times = size_res2), 
                         rep(0, times = size_res2), 
                         rep(0, times = size_res2), 
                         rep(1, times = size_res2), 
                         rep(1, times = size_res2))
        } else {
          if (mod_full$gform == "palm") {
            x_new <- cbind(rep(seq(0, 1, length = size_res), times = size_res),
                           rep(seq(0, 1, length = size_res), each = size_res),
                           rep(0, times = size_res2), 
                           rep(0, times = size_res2), 
                           rep(1, times = size_res2), 
                           rep(1, times = size_res2))
          } else {
            x_new <- cbind(rep(seq(0, 1, length = size_res), times = size_res),
                           rep(seq(0, 1, length = size_res), each = size_res),
                           rep(0, times = size_res2), 
                           rep(1, times = size_res2), 
                           rep(1, times = size_res2))
          }
        }
      }
      x_new[, (col_id + 1)] <- rep(trait_seq[i], times = size_res2)
      xx_new <- matrix(rep(seq(0, 1, length = size_res), times = size_res),
                       ncol = size_res,
                       byrow = TRUE)
      yy_new <- matrix(rep(seq(0, 1, length = size_res), times = size_res),
                       ncol = size_res,
                       byrow = TRUE)
    } else {
      if ((mod_full$gform == "tree") | (mod_full$gform == "shrub")) {
        x_new <- cbind(seq(0, 1, length = size_res),
                       rep(0, times = size_res), 
                       rep(0, times = size_res), 
                       rep(0, times = size_res), 
                       rep(0, times = size_res), 
                       rep(1, times = size_res), 
                       rep(1, times = size_res))
      } else {
        if (mod_full$gform == "herb") {
          x_new <- cbind(seq(0, 1, length = size_res),
                         rep(0, times = size_res), 
                         rep(0, times = size_res), 
                         rep(0, times = size_res), 
                         rep(1, times = size_res), 
                         rep(1, times = size_res))
        } else {
          if (mod_full$gform == "palm") {
            x_new <- cbind(seq(0, 1, length = size_res),
                           rep(0, times = size_res), 
                           rep(0, times = size_res), 
                           rep(1, times = size_res), 
                           rep(1, times = size_res))
          } else {
            x_new <- cbind(seq(0, 1, length = size_res),
                           rep(0, times = size_res), 
                           rep(1, times = size_res), 
                           rep(1, times = size_res))
          }
        }
      }
      x_new[, col_id] <- rep(trait_seq[i], times = size_res)
      xx_new <- seq(0, 1, length = size_res)
      yy_new <- NULL
    }
    out_tmp <- m$predict_y(x_new)
    mean <- out_tmp[[1]]
    var <- out_tmp[[2]]
    if (mod_full$param == "growth") {
      zz_mean[, , i] <- matrix(mean, ncol = size_res,
                               nrow = size_res, byrow = TRUE)
      zz_var[, , i] <- matrix(var, ncol = size_res,
                              nrow = size_res, byrow = TRUE)
    } else {
      zz_mean[i, ] <- mean
      zz_var[i, ] <- var
    }
  }
  return(list(z = zz_mean, z_var = zz_var, x = xx_new, y = yy_new))
}

# plot fitted GP model with kernels as a function of traits/predictors
plot_gpmod <- function(m,
                       trait_name = "sla",
                       size_res = 100, seq_len = NULL,
                       trait_min = -3., trait_max = 3.,
                       trait_mean_real = 20, trait_sd_real = 5,
                       col_palette = "inferno") {
  old_mfrow <- par()$mfrow
  old_mar <- par()$mar
  if (m$gform == "all") {
    if (trait_name != "all") {
      warning("model was fitted to all growth forms; predict function will default to trait_name = "all"",
              call. = FALSE)
      trait_name <- "all"
    }
    if (!is.null(seq_len)) {
      if (seq_len != 4) {
        warning("predict function will set seq_len = 4 for model of all growth forms",
                call. = FALSE)
      }
    }
    seq_len <- 4
    n_col <- 2
    n_row <- 2
  } else {
    if (is.na(match(trait_name, c("sla", "seed_mass", "ht", "wood_den")))) {
      warning(paste0(""", trait_name, """, " is not a valid trait_name; trait_name = "seed_mass" by default"),
              call. = FALSE)
      trait_name <- "seed_mass"
    }
    if (is.null(seq_len)) {
      if (m$param == "growth") {
        seq_len <- 4
      } else {
        seq_len <- 100
      }
    }
    n_col <- ceiling(seq_len / 3)
    n_row <- ceiling(seq_len / n_col)
    if (m$gform == "tree") {
      trait_mean_real <- switch(trait_name,
                                "sla" = 13.3,
                                "seed_mass" = 0.55,
                                "ht" = 20.0,
                                "wood_den" = 0.52)
      trait_sd_real <- switch(trait_name,
                              "sla" = 3.4,
                              "seed_mass" = 0.1,
                              "ht" = 6.,
                              "wood_den" = 0.15)
    } else {
      if (m$gform == "shrub") {
        trait_mean_real <- switch(trait_name,
                                  "sla" = 13.6,
                                  "seed_mass" = 2.14,
                                  "ht" = 3.75,
                                  "wood_den" = 0.62)
        trait_sd_real <- switch(trait_name,
                                "sla" = 4.,
                                "seed_mass" = 0.5,
                                "ht" = 0.9,
                                "wood_den" = 0.1)
      } else {
        if (m$gform == "herb") {
          trait_mean_real <- switch(trait_name,
                                    "sla" = 20.4,
                                    "seed_mass" = 0.01,
                                    "ht" = 0.7)
          trait_sd_real <- switch(trait_name,
                                  "sla" = 2.5,
                                  "seed_mass" = 0.001,
                                  "ht" = 0.2)
        } else {
          trait_mean_real <- switch(trait_name,
                                    "seed_mass" = 0.23,
                                    "ht" = 13.3)
          trait_sd_real <- switch(trait_name,
                                  "seed_mass" = 0.05,
                                  "ht" = 3.5)
        }
      }
    }
  }
  plot_vals <- predict_gp(m, trait_name = trait_name, size_res = size_res,
                          seq_len = seq_len, trait_min = trait_min,
                          trait_max = trait_max)
  trait_name_plot <- switch(trait_name,
                            "sla" = "Specific leaf area",
                            "seed_mass" = "Seed mass",
                            "ht" = "Maximum height",
                            "wood_den" = "Wood density")
  gform_names <- c("Trees", "Shrubs", "Herbs", "Palms")
  if (m$param == "growth") {
    par(mfrow = c(n_row, n_col), mar = c(4, 4.5, 2.5, 1.5))
    if (trait_name == "all") {
      trait_seq <- c(0., 1., 2., 3.)
    } else {
      trait_seq <- seq(trait_min, trait_max, length = seq_len)
    }
    for (i in seq(along = trait_seq)) {
      if (trait_name == "all") {
        trait_label <- gform_names[i]
      } else {
        trait_label <- paste0(trait_name_plot, " = ")
      }
      image(plot_vals$z[, , i], col = get(col_palette)(256), las = 1)
      if (!is.na(match(i, seq(1, size_res, by = n_col)))) {
        mtext("Relative size (t + 1)", side = 2, line = 3)
      }
      if (!is.na(match(i, seq(seq_len - n_col + 1, seq_len)))) {
        mtext("Relative size (t)", side = 1, line = 2.5)
      }
      if (trait_name != "all") {
        trait_label <- paste(trait_label, round(trait_mean_real + (trait_sd_real * trait_seq[i]), 2), sep = "")
        mtext(trait_label, side = 3, line = 0.5, adj = 1)
      } else {
        mtext(gform_names[i], side = 3, line = 0.5, adj = 1)
      }
    }
  } else {
    if (trait_name != "all") {
      par(mfrow = c(1, 1))
      trait_label <- trait_name_plot
      trait_seq <- seq(trait_min, trait_max, length = seq_len)
      image(x = plot_vals$x,
            y = (trait_mean_real + (trait_sd_real * trait_seq)),
            z = plogis(t(plot_vals$z)), col = get(col_palette)(256),
            las = 1, xlab = "", ylab = "")
      mtext(trait_label, side = 2, line = 2.5)
      mtext("Relative size", side = 1, line = 2.5)
    } else {
      par(mfrow = c(2, 2))
      if (m$param == "fec") {
        y_label <- "Number of seeds"
      } else {
        y_label <- "Survival probability"
      }
      for (i in c(1:4)) {
        if (m$param == "fec") {
          mid <- exp(plot_vals$z[i, ])
          upper <- exp(plot_vals$z[i, ] + sqrt(plot_vals$z_var[i, ]))
          lower <- exp(plot_vals$z[i, ] - sqrt(plot_vals$z_var[i, ]))
        } else {
          mid <- plogis(plot_vals$z[i, ])
          upper <- plogis(plot_vals$z[i, ] + sqrt(plot_vals$z_var[i, ]))
          lower <- plogis(plot_vals$z[i, ] - sqrt(plot_vals$z_var[i, ]))
        }
        ylim_set <- c(min(c(mid, lower, upper)),
                      max(c(mid, lower, upper)))
        col_set <- get(col_palette)(256)
        plot(mid ~ plot_vals$x,
             type = "n", bty = "l",
             xlab = "", ylab = "",
             ylim = ylim_set, las = 1)
        polygon(c(plot_vals$x, rev(plot_vals$x)),
                c(upper, rev(lower)),
                border = NA,
                col = alpha(col_set[1], 0.5))
        lines(mid ~ plot_vals$x, col = col_set[1], lwd = 2)
        if ((i == 1) | (i == 3)) {
          mtext(y_label, side = 2, line = 2.5)
        }
        if ((i == 3) | (i == 4)) {
          mtext("Relative size", side = 1, line = 2.5)
        }
        mtext(gform_names[i], side = 3, line = 1, adj = 0)
      }
    }
  }
  par(mfrow = old_mfrow, mar = old_mar)
}

# calculate full IPM kernel from fitted growth, fecundity and survival functions
kernel_calc <- function(kern_grow, kern_fec, kern_surv) {
  kerng <- sweep(kern_grow, 1, apply(kern_grow, 2, sum), "/")
  kernf <- matrix(0, nrow = nrow(kern_grow), ncol = ncol(kern_grow))
  kernf[1, ] <- kern_fec
  kerng2 <- kerng * kern_surv
  kern_out <- kerng2 + kernf
  return(kern_out)
}

# define full model run as a function of growth form and parameter
full_gp_run <- function(i, model_opt,
                        n_iter = 1000, n_burn = 1000,
                        n_cv = 10,
                        n_iter_cv = 1000, n_burn_cv = 1000,
                        n_induce_all = n_induce_all,
                        eps_all = eps_all,
                        lmin_all = lmin_all, lmax_all = lmax_all,
                        max_it_optim = 1000,
                        plot_outputs = FALSE,
                        par_run = FALSE) {
  gform <- model_opt$gform[i]
  param <- model_opt$param[i]
  n_induce <- n_induce_all[i]
  eps <- eps_all[i]
  lmin <- lmin_all[i]
  lmax <- lmax_all[i]
  data_set <- load_data(growth_form = gform, param_type = param, n_induce = n_induce)
  if (n_cv == "loo") {
    ncv_set <- length(unique(data_set$X[, ncol(data_set$X) - 1]))
  } else {
    ncv_set <- n_cv
  }
  print(paste0("This model will be cross validated with ", ncv_set, " folds"))
  mod_run <- fit_gp_mod(data_set,
                        n_iter = n_iter, n_burn = n_burn,
                        eps = eps, lmin = lmin, lmax = lmax,
                        max_it_optim = max_it_optim)
  mod_cv <- gp_cv(mod_run, n_cv = ncv_set,
                  n_iter = n_iter_cv, n_burn = n_burn_cv,
                  eps = eps, lmin = lmin, lmax = lmax,
                  max_it_optim = max_it_optim,
                  par_run = par_run)
  if ((gform == "tree") | (gform == "shrub")) {
    trait_list <- c("sla", "seed_mass", "ht", "wood_den")
    if (gform == "tree") {
      trait_mean_set <- c(13.3, 0.55, 20., 0.52)
      trait_sd_set <- c(3.4, 0.1, 6., 0.15)
    } else {
      trait_mean_set <- c(13.6, 2.14, 3.75, 0.62)
      trait_sd_set <- c(4., 0.5, 0.9, 0.1)
    }
  } else {
    if (gform == "herb") {
      trait_list <- c("sla", "seed_mass", "ht")
      trait_mean_set <- c(20.4, 0.01, 0.7)
      trait_sd_set <- c(2.5, 0.001, 0.2)
    } else {
      if (gform == "palm") {
        trait_list <- c("seed_mass", "ht")
        trait_mean_set <- c(0.23, 13.3)
        trait_sd_set <- c(0.05, 3.5)
      } else {
        trait_list <- c("all")
        trait_mean_set <- c(10., 1.)
        trait_sd_set <- c(5., 0.5)
      }
    }
  }
  if (plot_outputs) {
    if (param == "growth") {
      size_res_set <- 25
      seq_len_set <- 9
    } else {
      size_res_set <- 100
      seq_len_set <- 100
    }
    trait_data_min <- -3.
    trait_data_max <- 3.
    if (gform != "all") {
      for (i in seq(along = trait_list)) {
        plot_gpmod(mod_run, trait_name = trait_list[i],
                   size_res = size_res_set, seq_len = seq_len_set,
                   trait_min = trait_data_min, trait_max = trait_data_max,
                   trait_mean_real = trait_mean_set[i], trait_sd_real = trait_sd_set[i])
      }
    } else {
      plot_gpmod(mod_run, trait_name = "all",
                 size_res = size_res_set, seq_len = seq_len_set,
                 trait_min = -5., trait_max = 5.,
                 trait_mean_real = 10., trait_sd_real = 5.)
    }
  }
  return(list(mod = mod_run,
              mod_cv = mod_cv))
}

# summarise fitted GP model (extract fitted and cross-validated values, and 
#    model fit statistics)
gp_mod_summary <- function(mod_sum) {
  mod_summary <- extract_outputs(mod_sum$mod) 
  cv_outputs <- data.frame("pred" = mod_sum$mod_cv$y_pred, "real" = mod_sum$mod_cv$y_real)
  kernel_summary <- data.frame("mean" = mod_summary$kern_mean, "sd" = mod_summary$kern_sd)
  fitted_vals <- mod_summary$y_est
  r2cv_val <- cor(cv_outputs$pred, cv_outputs$real, use = "complete") ** 2
  mod_fit <- c(mod_summary$devi, mod_summary$r2, r2cv_val)
  names(mod_fit) <- c("deviance", "r2", "r2cv")
  return(list(summary = mod_summary,
              cv_out = cv_outputs,
              kernels = kernel_summary,
              fitted = fitted_vals,
              fit_stat = mod_fit))
}

# internal lambda calculation use by lambda_full_calc()
lambda_calc <- function(m_grow, m_fec, m_surv,
                        trait_name = "sla", size_res = 20, seq_len = 10,
                        trait_min = -3, trait_max = 3) {
  pred_vals_grow <- predict_gp(m_grow$mod, trait_name = trait_name, size_res = size_res,
                               seq_len = seq_len, trait_min = trait_min,
                               trait_max = trait_max)
  pred_vals_fecd <- predict_gp(m_fec$mod, trait_name = trait_name, size_res = size_res,
                               seq_len = seq_len, trait_min = trait_min,
                               trait_max = trait_max)
  pred_vals_surv <- predict_gp(m_surv$mod, trait_name = trait_name, size_res = size_res,
                               seq_len = seq_len, trait_min = trait_min,
                               trait_max = trait_max)
  lambda_val <- rep(NA, dim(pred_vals_grow$z)[3])
  for (i in c(1:dim(pred_vals_grow$z)[3])) {
    kernel <- kernel_calc(pred_vals_grow$z[, , i],
                          pred_vals_fecd$z[i, ],
                          pred_vals_surv$z[i, ])
    lambda_val[i] <- Re(eigen(kernel)$values)[1]
  }
  return(lambda_val)
}

# optional helper to calculate long-term population growth rate from 
#   estimated IPM kernel
# calculate lambda from fitted IPM
lambda_full_calc <- function(gform = "tree",
                             mod_grow, mod_fec, mod_surv) {
  if ((gform == "tree") | (gform == "shrub")) {
    trait_list <- c("sla", "seed_mass", "ht", "wood_den")
  } else {
    if (gform == "herb") {
      trait_list <- c("sla", "seed_mass", "ht")
    } else {
      trait_list <- c("seed_mass", "ht")
    }
  }
  size_res_set <- 25
  seq_len_set <- 25
  trait_min_set <- rep(-3, length(trait_list))
  trait_max_set <- rep(3, length(trait_list))
  lambda_out <- matrix(NA, nrow = length(trait_list), ncol = seq_len_set)
  for (i in seq(along = trait_list)) {
    lambda_out[i, ] <- lambda_calc(m_grow = mod_grow, m_fec = mod_fec, m_surv = mod_surv,
                                   trait_name = trait_list[i],
                                   size_res = size_res_set,
                                   seq_len = seq_len_set,
                                   trait_min = trait_min_set[i],
                                   trait_max = trait_max_set[i])
  }
  return(lambda_out)
}
