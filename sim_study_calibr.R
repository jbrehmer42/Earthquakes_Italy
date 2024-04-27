library(monotone)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(geomtextpath)     # for geom_labelabline
library(data.table)

fpath <- "./figures"
tpath <- "./save_results"

source("data_prep.R")    # make sure that LM model is loaded
source("utils.R")

################################################################################

scf_list <- list(quad = s_quad_da, pois = s_pois_da)

get_forecaster <- function(fcst, daily = F) {
  if (!daily) {
    y <- as.vector(obs)
    x <- as.vector(models[[which(model_names == "LM")]])
    discretize <- quantile(x, probs = c(0, 0.1, 0.4, 0.5, 0.6, 0.9, 1))
  } else {
    discretize <- 10^((-2:2) / 2)
    y <- rowSums(obs)
    x <- rowSums(models[[1]])
  }

  if (fcst == "lm") {
    x <- x
  } else if (fcst == "lm_rc") {
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x_rc <- monotone(y[ord])
    x <- x_rc[order(ord)]     # restore original sorting
  } else if (fcst == "lm_x5") {
    x <- 4 * x
  } else if (fcst == "lm_x0_2") {
    x <- 0.25 * x
  } else if (fcst == "lm_upped") {
    x <- discretize[sapply(x, function(a) sum(a > discretize) + 1)]
  } else if (fcst == "lm_downed") {
    x <- discretize[sapply(x, function(a) sum(a >= discretize))]
  } else if (fcst == "lm_overconf") {
    x <- ifelse(x > 1e-5, 4 * x, 0.25 * x)
  } else if (fcst == "lm_underconf") {
    x <- ifelse(x > 4e-5, 0.25 * x, ifelse(x < 0.25 * 1e-5, 4 * x, 1e-5))
  }
  return(data.frame(y = y, x = x))
}

filter_jumps <- function(v) {
  n <- length(v)
  # filter points around jumps (points immediately before and after a jump)
  return(c(T, v[-1] - v[-n] > 0) | c(v[-n] - v[-1] < 0, T))
}

# Simulation study run ---------------------------------------------------------

run_simulation_study <- function(vec_fcst = c("lm", "lm_rc", "lm_x5", "lm_x0_2",
                                              "lm_upped", "lm_downed", "lm_underconf",
                                              "lm_overconf"),
                                 daily = FALSE) {
  scores <- data.frame()
  reliability <- data.frame()
  murphy <- data.frame()

  # grid for murphy diagram
  n_theta <- 100
  log_grid <- seq(-15, 4, len = n_theta)
  grd <- 10^log_grid

  for (forecaster in vec_fcst) {
    cat(forecaster)
    df <- get_forecaster(forecaster, daily = daily)
    x <- df$x
    y <- df$y
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x <- x[ord]
    y <- y[ord]
    x_rc <- monotone(y)
    mean_y <- mean(y)

    for (scf_name in names(scf_list)) {
      scf <- scf_list[[scf_name]]
      s <- scf(x, y)
      s_rc <- scf(x_rc, y)
      s_mg <- scf(mean_y, y)

      scores <- rbind(
        scores,
        data.frame(fcst = forecaster, Scoring = scf_name, Type = c("Score", "MCB", "DSC"),
                   value = c(s, s - s_rc, s_mg - s_rc))
      )
    }

    filt_jumps <- filter_jumps(x_rc)

    reliability <- rbind(
      reliability,
      data.frame(fcst = forecaster, x = x[filt_jumps], x_rc = x_rc[filt_jumps])
    )

    elem_scores <- sapply(grd, function(theta) S_theta(x, y, theta))
    elem_scores_rc <- sapply(grd, function(theta) S_theta(x_rc, y, theta))
    elem_scores_mg <- sapply(grd, function(theta) S_theta(mean_y, y, theta))

    murphy <- rbind(
      murphy,
      data.frame(fcst = forecaster, log_x = rep(log_grid, 3),
                 Type = rep(c("Score", "MCB", "DSC"), each = length(grd)),
                 y = c(elem_scores, elem_scores - elem_scores_rc,
                       elem_scores_mg - elem_scores_rc) / n_days)
    )
  }
  return(list(scores = scores, reliability = reliability, murphy = murphy))
}

# Plots ------------------------------------------------------------------------

change_names <- function(fcsts = NULL) {
  mapping <- c("lm" = "LM", "lm_rc" = "LM rc", "lm_x5" = "LM x4", "lm_x0_2" = "LM x0.25",
               "lm_upped" = "LM upped", "lm_downed" = "LM downed",
               "lm_overconf" = "LM overconf", "lm_underconf" = "LM underconf")
  if (is.null(fcsts)) {
    return(mapping)
  } else {
    return(mapping[fcsts])
  }
}

sort_forecasts <- function(fcsts) {
  return(factor(
    fcsts, ordered = T,
    levels = c("lm", "lm_rc", "lm_x5", "lm_x0_2", "lm_upped", "lm_downed",
               "lm_overconf", "lm_underconf")
  ))
}

get_score_cmp_plot <- function(results, non_sig_seg, daily = F) {
  scores_wide <- results %>%
    filter(is.finite(value)) %>%
    pivot_wider(id_cols = c(fcst, Scoring), names_from = Type, values_from = value) %>%
    mutate(Scoring = ifelse(Scoring == "pois", "Poisson", "Quadratic"))
  unc <- mean(with(scores_wide, Score - MCB + DSC), na.rm = T)

  fmt <- paste0("%.",  ifelse(scores_wide$Scoring[1] == "Poisson", 2, 3), "f")

  if (!daily) {
    justs <- data.frame(fcst = c("lm", "lm_rc", "lm_x5", "lm_x0_2", "lm_upped", "lm_downed", "lm_overconf", "lm_underconf"),
                        hjusts = c(-0.3, -0.2, 0.5, -0.1, 0.7, 0.0, 0.7, 0.0),
                        vjusts = c(-0.6, 1.3, 1.8, -0.7, 1.8, 1.8, -1.0, -1.0))
  } else {
    justs <- data.frame(fcst = c("lm", "lm_rc", "lm_x5", "lm_x0_2", "lm_upped", "lm_downed", "lm_overconf", "lm_underconf"),
                        hjusts = c(-0.3, -0.2, 0.5, -0.3, -0.2, 0.0, 0.0, 0.0),
                        vjusts = c(-0.6, 1.3, 1.8, 0.2, 0.2, 1.8, -1.0, -1.0))
  }

  my_hjust <- ifelse(daily, 0.38, 0.55)

  iso <- data.frame(
    intercept = seq(-max(scores_wide$MCB, na.rm = T), max(scores_wide$DSC, na.rm = T),
                    length.out = 10)) %>%
    mutate(label = sprintf(fmt, unc - intercept))    # miscalibration is 0, intercept corresponds to DSC


  mcb_range <- range(scores_wide$MCB)
  dsc_range <- range(scores_wide$DSC)

  pl <- left_join(scores_wide, justs, by = "fcst") %>%
    mutate(fcst = sort_forecasts(fcst)) %>%
    ggplot() +
    facet_wrap(~Scoring) +
    geom_abline(data = iso, aes(intercept = intercept, slope = 1.0), color = "lightgray",
                alpha = 0.5, size = 0.5) +
    geom_labelabline(data = iso, aes(intercept = intercept, slope = 1.0, label = label),
                     color = "gray50", hjust = my_hjust, size = 7 * 0.36, text_only = TRUE,
                     boxcolour = NA, straight = TRUE) +
    geom_segment(data = non_sig_seg, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = "dotted", alpha = 0.5, size = 0.5) +
    geom_point(aes(x = MCB, y = DSC, color = fcst), size = 1.5) +
    geom_text(aes(x = MCB, y = DSC, label = change_names(fcst), hjust = hjusts, vjust = vjusts,
                  color = fcst),
              size = 8 * 0.36) +
    scale_color_discrete() +
    coord_cartesian(xlim = mcb_range + c(0, 1) / 10 * diff(mcb_range),
                    ylim = dsc_range + c(-3, 1) / 10 * diff(dsc_range),
                    expand = F) +
    annotate("label", x = Inf, y = -Inf, label = paste0("UNC = ", sprintf("%.3f", unc)),
             hjust = 1.05, vjust = -0.2) +
    my_theme +
    theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  return(pl)
}

get_non_sig_connections <- function(df_score_cmp, level = 0.1, daily = F) {
  non_sig_con <- data.frame()
  fcsts <- unique(df_score_cmp$fcst)

  for (i in 1:(length(fcsts) - 1)) {
    t <- get_forecaster(fcsts[i], daily = daily)
    obs <- t$y
    model1 <- t$x
    for (j in ((i + 1):length(fcsts))) {
      model2 <- get_forecaster(fcsts[j], daily = daily)$x

      if (!daily) {
        obs <- matrix(obs, nrow = n_days)
        model1 <- matrix(model1, nrow = n_days)
        model2 <- matrix(model2, nrow = n_days)
      }

      p <- dm_test(model1, model2, obs, s_pois, daily = daily)$pval
      if ((p > level / 2) & (p < 1 - level / 2)) {
        non_sig_con <- rbind(
          non_sig_con,
          data.frame(x = filter(df_score_cmp, fcst == fcsts[i], Type == "MCB")$value,
                     y = filter(df_score_cmp, fcst == fcsts[i], Type == "DSC")$value,
                     xend = filter(df_score_cmp, fcst == fcsts[j], Type == "MCB")$value,
                     yend = filter(df_score_cmp, fcst == fcsts[j], Type == "DSC")$value)
        )
      }
    }
  }
  return(non_sig_con)
}

plot_score_components <- function(results, nonsig_seg, daily = F) {
  pois_results <- filter(results, Scoring == "pois")
  pl_pois <- get_score_cmp_plot(pois_results, nonsig_seg, daily = daily)
  return(pl_pois)
}

make_groups <- function(df) {
  df_groups <- rbind(
    filter(df, fcst %in% c("lm", "lm_rc")) %>% mutate(I = "A"),
    filter(df, fcst %in% c("lm", "lm_x5", "lm_x0_2")) %>% mutate(I = "B"),
    filter(df, fcst %in% c("lm", "lm_upped", "lm_downed")) %>% mutate(I = "C"),
    filter(df, fcst %in% c("lm", "lm_overconf", "lm_underconf")) %>% mutate(I = "D")
  ) %>%
    mutate(fcst = sort_forecasts(fcst))
  return(df_groups)
}

get_score_display <- function(df_scores) {
  # first coordinate: upper left, second coordinate: lower right
  text_x <- c(0.02, 0.68)
  text_y <- c(0.68, 0.02)
  push <- 0.03

  df_stats <- df_scores %>%
    filter(Scoring == "pois") %>%
    group_by(fcst) %>%
    # data has to be in the correct order!!!
    summarise(label = paste(c(Type, "UNC"), c("", " ", " ", " "),
                            sprintf("%.3f", c(value, sum(value * c(1, -1, 1)))),
                            collapse = "\n"), .groups = "drop")

  df_groups <- rbind(
    filter(df_stats, fcst %in% c("lm", "lm_rc")) %>% mutate(I = "A", x = rev(text_x), y = rev(text_y)),
    filter(df_stats, fcst %in% c("lm_x5", "lm_x0_2")) %>% mutate(I = "B", x = text_x + c(0, push),
                                                                 y = text_y + c(push, 0)),
    filter(df_stats, fcst %in% c("lm_upped", "lm_downed")) %>% mutate(I = "C", x = text_x, y = text_y),
    filter(df_stats, fcst %in% c("lm_overconf", "lm_underconf")) %>% mutate(I = "D", x = rev(text_x), y = rev(text_y))
  ) %>%
    mutate(fcst = sort_forecasts(fcst))
  return(df_groups)
}

get_ecdf_trans <- function(daily = FALSE) {
  col_ecdfs <- list()
  for (i in 1:length(models)) {
    if (daily) {
      t <- table(rowSums(models[[i]]))
      col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                                   y = c(0, cumsum(as.numeric(t))) / prod(dim(models[[i]])),
                                   M = model_names[i])
    } else {
      col_ecdfs[[i]] <- read.csv(file.path(tpath, paste0("ecdf_", model_names[i], ".csv")), row.names = 1)
    }
  }
  mean_ecdf <- do.call(rbind, col_ecdfs) %>%
    pivot_wider(id_cols = x, names_from = M, values_from = y, values_fill = NA) %>%
    arrange(x) %>%
    setnafill(type = "locf") %>%
    slice(floor(seq(1, nrow(.), length.out = 10^5))) %>%   # sample on grid for plotting
    transmute(x = x, y = rowMeans(select(., all_of(model_names))))

  my_trans <- function(x) {
    # first known values, then values we want to evaluate --> fill with last observation
    rbind(cbind(mean_ecdf, a = -1.0), data.table(x = x, y = NA, a = 1:length(x))) %>% arrange(x) %>%
    setnafill(type = "locf") %>%
    filter(a > 0) %>%
    arrange(a) %>%
    pull(y)
  }
  return(my_trans)
}

plot_rel <- function(reliability, score_cmps, use_ecdf = T, daily = F) {
  pick_panels <- c("A", "B", "D")
  n_rows <- 1
  plot_colors <- scales::hue_pal()(8)[1:6]
  plot_data <- make_groups(reliability) %>%
    filter(I %in% pick_panels)

  use_points <- c("lm_upped", "lm_downed", "lm_rc")
  df_points <- filter(plot_data, fcst %in% use_points)
  df_segments <- filter(plot_data, !(fcst %in% use_points)) %>%
    group_by(fcst, I) %>%
    group_modify(function(df, key) {
      df %>%
        mutate(x_pos = c("x0", ifelse(x_rc[-1] - x_rc[-nrow(df)] > 0, "x0", "x1"))) %>%
        pivot_wider(id_cols = "x_rc", names_from = "x_pos", values_from = "x", values_fill = NA) %>%
        filter(!is.na(x0), !is.na(x1))    # if there is one missing, throw it out
    }, .keep = T)
  df_scores <- get_score_display(score_cmps) %>%
    filter(I %in% pick_panels)

  if (use_ecdf) {
    my_trans <- get_ecdf_trans(daily = daily)
    my_breaks <- c(0, 10^c(-7, -6, -5, -4, 0))
    my_labels <- c("", expression(10^{-7}), "", expression(10^{-5}), expression(10^{-4}), "")
  } else {
    if (daily) {
      my_trans <- function(x) log(x, base = 10)
      my_breaks <- c(0.1, 1, 10)
      my_labels <- my_breaks
    } else {
      d <- 10^6
      my_trans <- function(x) sign(x) * log(abs(x) * d  + 1, base = 10)
      my_breaks <- c(0, 10^c(-6, -4, -2, 0))
      my_labels <- c("0", paste0("1e", c(-6, -4, -2)), "1")
    }
  }
  my_breaks <- my_trans(my_breaks)

  my_plot <- ggplot(plot_data) +
    facet_wrap(~I, nrow = n_rows) +
    geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
                linetype = "dashed") +
    geom_line(aes(x = my_trans(x), y = my_trans(x_rc), color = fcst), size = 0.3) +
    geom_point(data = df_points,
               aes(x = my_trans(x), y = my_trans(x_rc), color = fcst), size = 1) +
    geom_segment(data = df_segments,
                 aes(x = my_trans(x0), xend = my_trans(x1), y = my_trans(x_rc),
                     yend = my_trans(x_rc), color = fcst), size = 0.7) +
    geom_text(data = df_scores, mapping = aes(x = x, y = y, label = label, color = fcst),
              size = 6 * 0.36, hjust = 0, vjust = 0, show.legend = F) +
    xlab("Forecast value") +
    scale_x_continuous(breaks = my_breaks, labels = my_labels) +
    scale_y_continuous(breaks = my_breaks, labels = my_labels) +
    ylab("Conditional mean") +
    scale_color_manual(name = NULL, labels = change_names(), values = plot_colors,
                       guide = guide_legend(nrow = 1)) +
    my_theme +
    theme(strip.text = element_blank(), legend.position = "bottom",
          aspect.ratio = 1, panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(my_plot)
}

# Analysis ---------------------------------------------------------------------

set.seed(111)
cmp_run_sim <- compiler::cmpfun(run_simulation_study)
daily <- FALSE

l_results <- cmp_run_sim(daily = daily)

my_plot <- plot_rel(l_results$reliability, l_results$scores, daily = daily, use_ecdf = T)
file_path <- file.path(fpath, "Fig6_RelDiag-manipulated.pdf")
ggsave(file_path, width = 140, height = 68, unit = "mm", plot = my_plot)

non_sig_seg <- get_non_sig_connections(filter(l_results$scores$Scoring == "pois"))
my_plot <- plot_score_components(l_results$scores, non_sig_seg, daily = daily)
file_path <- file.path(fpath, "Fig7_MCB-DSC-manipulated-seg.pdf")
ggsave(file_path, width = 140, height = 80, unit = "mm", plot = my_plot)

write.csv(l_results$scores, file.path(tpath, "simstudy-rel_scores.csv"))
write.csv(l_results$reliability, file.path(ftath, "simstudy-rel_rel.csv"))
