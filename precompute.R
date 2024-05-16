library(dplyr)
library(tidyr)
library(data.table)
library(monotone)         # for monotone : fast isotonic mean regression

source("data_prep.R")
source("utils.R")

tpath <- "./save_results"

if (!dir.exists(tpath)) {
  dir.create(tpath)
}

################################################################################
# ECDF transform
################################################################################

col_ecdfs <- list()
for (i in 1:length(models)) {
  t <- table(as.vector(models[[i]]))
  col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                               y = c(0, cumsum(as.numeric(t))) / prod(dim(models[[i]])),
                               M = model_names[i]) %>%
    slice(floor(seq(1, nrow(.), length.out = pmin(nrow(.), 10^5))))    # subsample to decrease computational complexity
  write.csv(col_ecdfs[[i]], file.path(tpath, paste0("ecdf_", model_names[i], ".csv")))
}
rm(col_ecdfs)

################################################################################
# Murphy diagramm
################################################################################

S_elem <- compiler::cmpfun(s_theta_da) # compile function to reduce runtime a bit

daily <- FALSE    # for daily = T: aggregate data manually before hand

n_theta <- 100
if (!daily) {
  log_grid <- seq(-24, 4, len = n_theta)
} else {
  log_grid <- seq(-7, 4, len = n_theta)
}
grd <- exp(log_grid)

# use for loops, as otherwise memory demand is too high
# Jonas calculates each row separately, it seems faster? check it!
murphy_df <- matrix(0, nrow = length(grd), ncol = length(models))
for (m in 1:length(models)) {
  print(m)
  for (t in 1:n_theta) {
    cat("*")
    murphy_df[t, m] <- S_elem(models[[m]], obs, grd[t])
  }
}
colnames(murphy_df) <- model_names

write.csv(murphy_df, file.path(tpath, "murphy_df.csv"))

# now look at Murphy diagram of miscalibration and discrimination component
MCB_diag <- DSC_diag <- matrix(0.0, ncol = n_mods, nrow = n_theta)
UNC_diag <- numeric(length(grd))

for (i in 1:n_mods) {
  print(i)
  # recalibrate forecasts x with isotoinc regression from monotone package
  x <- as.vector(models[[i]])
  y <- as.vector(obs)
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)

  for (t in 1:length(grd)) {
    cat("*")
    s <- S_elem(x, y, grd[t])
    s_rc <- S_elem(x_rc, y, grd[t])
    s_mg <- S_elem(mean(y), y, grd[t])

    MCB_diag[t, i] <- s - s_rc
    DSC_diag[t, i] <- s_mg - s_rc
    if (i == 1) {
      UNC_diag[t] <- s_mg
    }
  }
}

df_collect <- rbind(
  cbind(data.frame(MCB_diag), Type = "MCB"),
  cbind(data.frame(DSC_diag), Type = "DSC")
)
colnames(df_collect) <- c(model_names, "Type")

df_collect <- df_collect %>%
  mutate(log_theta = rep(log_grid, 2)) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  mutate(Type = ifelse(Type == "MCB", "Miscalibration", "Discrimination"))

write.csv(df_collect, file.path(tpath, "murphy_mcb_dsc.csv"))
write.csv(data.frame(unc = UNC_diag), file.path(tpath, "murphy_unc.csv"))

rm(df_collect, murphy_df, MCB_diag, DSC_diag, UNC_diag)

################################################################################
# Score component plot
################################################################################

get_score_cmps <- function(x, y, name, scf_list = list("pois" = s_pois_da, "quad" = s_quad_da), daily = F) {
  if (!daily) {
    x <- as.vector(x)
    y <- as.vector(y)
  } else {
    x <- rowSums(x)
    y <- rowSums(y)
  }

  scores <- data.frame()
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)

  for (scf_name in names(scf_list)) {
    scf <- scf_list[[scf_name]]
    s <- scf(x, y)
    s_rc <- scf(x_rc, y)
    s_mg <- scf(mean(y), y)

    scores <- rbind(
      scores,
      data.frame(Model = name, Scoring = scf_name, Type = c("Score", "MCB", "DSC"),
                 value = c(s, s - s_rc, s_mg - s_rc))
    )
  }
  return(scores)
}

get_non_sig_connections <- function(scf, level = 0.1) {
  non_sig_con <- data.frame()

  for (i in 1:(length(models) - 1)) {
    for (j in ((i + 1):length(models))) {
      p <- dm_test(models[[i]], models[[j]], obs, scf)$pval
      if ((p > level / 2) & (p < 1 - level / 2)) {
        non_sig_con <- rbind(non_sig_con, data.frame(m1 = model_names[i], m2 = model_names[j]))
      }
    }
  }
  return(non_sig_con)
}

df_score_cmp <- do.call(
  rbind,
  lapply(1:length(models), function(i) get_score_cmps(models[[i]], obs, model_names[i], daily = F))
)
write.csv(df_score_cmp, file.path(tpath, "score-cmps.csv"))

df_score_cmp <- do.call(
  rbind,
  lapply(1:length(models), function(i) get_score_cmps(models[[i]], obs, model_names[i], daily = T))
)
write.csv(df_score_cmp, file.path(tpath, "score-cmps_daily.csv"))

non_sig_seg_pois <- get_non_sig_connections(s_pois)
write.csv(non_sig_seg_pois, file.path(tpath, "score-cmps_non-sig-seg.csv"))

rm(df_score_cmp, non_sig_seg_pois)

################################################################################
# Reliability diagramm
################################################################################

reldiag <- function(x, y, n_resamples = 99, region_level = 0.9) {
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]

  # compress fit by only storing values around jumps that belong to distinct x-values
  # (if there is a horizontal section of width 0 with more than two data points,
  #  we would have two identical points in the data though we need only one)
  filter_jumps <- function(x, v) {
    n <- length(v)
    jumps <- (v[-1] - v[-n] > 0)
    # filter points around jumps (points immediately before and after a jump)
    next_to_jump <- c(T, jumps) | c(jumps, T)
    x_at_jump <- x[next_to_jump]
    # set to False if corresponding x-values do not change
    next_to_jump[next_to_jump] <- c(T, x_at_jump[-1] - x_at_jump[-length(x_at_jump)] > 0)
    return(next_to_jump)
  }

  score <- s_pois_da

  x_rc <- monotone(y)
  s <- score(x,y)
  s_rc <- score(x_rc, y)
  s_mg <- score(mean(y), y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  jumps <- filter_jumps(x, x_rc)
  rc_fit <- data.table(x = x[jumps], y = x_rc[jumps])

  t <- table(y) / length(y)
  vals <- as.numeric(names(t))
  f_y <- as.numeric(t)
  eps <- sum(f_y[-1] * vals[-1]) / x - 1
  n_split <- 5      # as not everything fits into memory, split
  split <- ceiling(seq(1, length(y), length.out = n_split))
  cmp <- matrix(rep(cumsum(f_y[-length(f_y)]), split[2] - split[1] + 1),
                ncol = length(f_y) - 1, byrow = T)

  collect_vals <- list()
  for (i in 1:n_resamples) {
    for (j in 2:n_split) {
      i_vec <- split[j - 1]:split[j]
      u <- runif(length(i_vec), 0, 1 + eps[i_vec])
      y[i_vec] <- rowSums(u - eps[i_vec] >= cmp[1:pmin(length(i_vec), nrow(cmp)), ])
    }
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x_rc <- monotone(y[ord])
    jumps <- filter_jumps(x, x_rc)
    # x is already sorted, so x[ord] would just resort x on sections on which it
    # is already constant
    collect_vals[[i]] <- data.table(x = x[jumps], y = x_rc[jumps], I = paste0("R", i))
  }
  x_grid <- sort(unique(c(do.call(c, lapply(collect_vals, function(vals) vals$x)), rc_fit$x)))
  # linearly interpolate resampled fits on x values of original fit
  resampl_fits <- do.call(cbind, lapply(collect_vals, function(df) {
    approx(x = df$x, y = df$y, xout = x_grid, method = "linear", yleft = NA, yright = NA)$y
  }))
  val_low_and_up <- apply(resampl_fits, 1, function(row) {
      sorted <- sort(row)
      low <- floor(length(sorted) * (1 - region_level) / 2)   # maybe there are NA, so recalculate it
      up <- length(sorted) - low
      low <- pmax(low, 1)
      c(sorted[low], sorted[up])
    }) %>% t() %>%  # apply writes results in columns
    as.data.table()

  results <- cbind(
    data.table(x = x_grid,
               x_rc = approx(x = rc_fit$x, y = rc_fit$y, xout = x_grid, method = "linear", yleft = NA, yright = NA)$y),
    val_low_and_up
  )
  colnames(results) <- c("x", "x_rc", "lower", "upper")
  stats <- data.frame(Score = s, MCB = mcb, DSC = dsc, UNC = unc)
  return(list(results = results, stats = stats))
}

reldiag_cmp <- compiler::cmpfun(reldiag)

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)   # for resampling fix seed

for (i in 1:length(models)) {
  print(paste(model_names[i], Sys.time()))
  res <- reldiag_cmp(as.vector(models[[i]]), as.vector(obs), n_resamples = 100)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.3f", res$stats[1, ]),
                        collapse = "\n"))
  )
}

write.csv(recal_models, file.path(tpath, "df_rel-Til100.csv"))
write.csv(collect_stats, file.path(tpath, "df_rel-stats.csv"))

rm(recal_models, collect_stats)
