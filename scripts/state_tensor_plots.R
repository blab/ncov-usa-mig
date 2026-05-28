# File: state_tensor_plots.R
# Author(s): Amin Bemanian
# Date: 5/21/26
# Description: Tensor product GAMs (NB count and RR form) for state-pair
#              identical sequence counts vs. travel covariates and CBSA distance.
#              Can be run standalone or sourced from state_regression.R.
#              When sourced, df_model_data/scenario/results_file are inherited.
#              When run standalone, they are built from data files.
# Arguments: --scenario (default: CAM_1000)
# Output: figs/<scenario>/dist/gam_*_te.jpg
#         figs/<scenario>/dist/gam_fit.jpg
#         figs/<scenario>/dist/nb_gam_dev_decomp.jpg
#         results/<scenario>/tensor_results.txt

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(argparse)
  library(mgcv)
  library(patchwork)
  library(yardstick)
})

# ============================================================================
# ARGUMENTS AND DATA LOADING
# ============================================================================

parser <- ArgumentParser()
parser$add_argument("--scenario", type = "character", default = "CAM_1000")
args     <- parser$parse_args()
scenario <- args$scenario

results_file <- paste0("results/", scenario, "/tensor_results.txt")
if (file.exists(results_file)) file.remove(results_file)

dir.create(paste0("figs/", scenario, "/dist"), recursive = TRUE, showWarnings = FALSE)

normalize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

state_rr <- fread(paste0("results/", scenario, "/df_RR_by_state.tsv")) %>%
  select(x, y, euclid_dist, min_cbsa_dist, RR, N_pairs, N_x, N_y, N_total) %>%
  rename(RR_seq = RR)

df_travel <- fread("data/travel_vars.tsv") %>%
  inner_join(state_rr, by = c("x", "y")) %>%
  filter(x != y, !is.na(RR_seq), RR_seq > 0)

df_model_data <- df_travel %>%
  filter(min_cbsa_dist <= 5000) %>%
  mutate(
    log_RR_trips    = log10(RR_trips),
    log_RR_air_t100 = log10(RR_air_t100),
    log_RR_air_db1b = log10(RR_air_db1b),
    log_RR_move     = log10(RR_move)
  )

message("\n", strrep("=", 80))
message("TENSOR PRODUCT GAMs: DISTANCE x TRAVEL")
message(strrep("=", 80), "\n")

# ============================================================================
# DATA PREPARATION
# ============================================================================

make_model_df <- function(df, rr_col, count_col, floor_val = 0, sd_trim = 3) {
  out <- df %>%
    filter(N_pairs > 0, .data[[count_col]] > floor_val, is.finite(.data[[rr_col]])) %>%
    mutate(log_offset = log(as.numeric(N_x) * as.numeric(N_y) / as.numeric(N_total)))
  mu    <- mean(out[[rr_col]])
  sigma <- sd(out[[rr_col]])
  out %>% filter(abs(.data[[rr_col]] - mu) <= sd_trim * sigma)
}

# Airline floor = 2 (the +1 correction applied to both directions means 0 real flights → value of 2)
df_trips_data <- make_model_df(df_model_data, "log_RR_trips",    "trips_xy",    floor_val = 0)
df_t100_data  <- make_model_df(df_model_data, "log_RR_air_t100", "pass_xy_t100", floor_val = 2)
df_db1b_data  <- make_model_df(df_model_data, "log_RR_air_db1b", "pass_xy_db1b", floor_val = 2)
df_move_data  <- make_model_df(df_model_data, "log_RR_move",     "n_move_avg",   floor_val = 0)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

plot_pred_nb_gam <- function(m, df, model_name) {
  rr_obs <- df$N_pairs / exp(df$log_offset)
  rr_hat <- fitted(m) / exp(df$log_offset)
  eps <- 1e-10
  rmse_val <- as.numeric(yardstick::rmse_vec(log10(rr_obs + eps), log10(rr_hat + eps)))
  message(sprintf("GAM log10 RMSE: %.3f", rmse_val))
  plot_df <- data.frame(rr_obs = rr_obs + eps, rr_hat = rr_hat + eps)
  range_vals <- range(c(plot_df$rr_obs, plot_df$rr_hat), na.rm = TRUE)
  ggplot(plot_df, aes(x = rr_hat, y = rr_obs)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(transform = "log10", name = "Predicted RR") +
    scale_y_continuous(transform = "log10", name = "Observed RR") +
    ggtitle(sprintf("%s - log10 RMSE: %.3f", model_name, rmse_val)) +
    coord_equal(xlim = range_vals, ylim = range_vals) +
    theme_bw()
}

plot_pred_rr_gam <- function(m, df, model_name) {
  rr_obs <- 10^log10(df$RR_seq)
  rr_hat <- 10^fitted(m)
  eps <- 1e-10
  rmse_val <- as.numeric(yardstick::rmse_vec(log10(rr_obs + eps), log10(rr_hat + eps)))
  message(sprintf("GAM RR log10 RMSE: %.3f", rmse_val))
  plot_df <- data.frame(rr_obs = rr_obs + eps, rr_hat = rr_hat + eps)
  range_vals <- range(c(plot_df$rr_obs, plot_df$rr_hat), na.rm = TRUE)
  ggplot(plot_df, aes(x = rr_hat, y = rr_obs)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(transform = "log10", name = "Predicted RR") +
    scale_y_continuous(transform = "log10", name = "Observed RR") +
    ggtitle(sprintf("%s - log10 RMSE: %.3f", model_name, rmse_val)) +
    coord_equal(xlim = range_vals, ylim = range_vals) +
    theme_bw()
}

write_gam_results <- function(m, label) {
  if (is.null(m)) return(invisible(NULL))
  sink(results_file, append = TRUE)
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat(paste0("TENSOR PRODUCT GAM: ", label, "\n"))
  cat(paste0(strrep("=", 80), "\n\n"))
  print(summary(m))
  cat(sprintf("\nDev. explained: %.1f%%\n", summary(m)$dev.expl * 100))
  cat(sprintf("AIC: %.2f\n", AIC(m)))
  sink()
  message(sprintf("GAM %s — Dev. expl.: %.1f%%  AIC: %.2f",
                  label, summary(m)$dev.expl * 100, AIC(m)))
}

plot_te_surface <- function(m, df, xvar, zvar, n = 60, xlabel, zlabel, title, convert_log = TRUE) {
  xseq <- seq(min(df[[xvar]], na.rm = TRUE), max(df[[xvar]], na.rm = TRUE), length.out = n)
  zseq <- seq(min(df[[zvar]], na.rm = TRUE), max(df[[zvar]], na.rm = TRUE), length.out = n)
  grid <- expand.grid(setNames(list(xseq, zseq), c(xvar, zvar)))
  for (col in setdiff(names(df), c(xvar, zvar))) {
    if (is.numeric(df[[col]])) grid[[col]] <- median(df[[col]], na.rm = TRUE)
    else if (is.factor(df[[col]])) grid[[col]] <- levels(df[[col]])[1]
    else grid[[col]] <- df[[col]][1]
  }
  terms_mat <- predict(m, newdata = grid, type = "terms")
  interaction_col <- grep(paste0("(?=.*", xvar, ")(?=.*", zvar, ")"),
                          colnames(terms_mat), perl = TRUE, value = TRUE)
  smooth_name <- if (length(interaction_col) > 0) interaction_col[1] else
                   grep(zvar, colnames(terms_mat), value = TRUE)[1]
  grid$partial <- as.numeric(terms_mat[, smooth_name]) / ifelse(convert_log, log(10), 1)
  too_far <- mgcv::exclude.too.far(grid[[xvar]], grid[[zvar]],
                                    df[[xvar]],   df[[zvar]], dist = 0.3)
  grid$partial[too_far] <- NA

  x_range <- range(xseq)
  z_range <- range(zseq)
  z_breaks <- seq(ceiling(z_range[1]), floor(z_range[2]))

  p_main <- ggplot(grid, aes(x = .data[[xvar]], y = .data[[zvar]], fill = partial)) +
    geom_tile() +
    geom_contour(aes(z = partial), color = "white", alpha = 0.5, linewidth = 0.3) +
    geom_point(data = df, aes(x = .data[[xvar]], y = .data[[zvar]]),
               inherit.aes = FALSE, size = 0.6, alpha = 0.2, color = "black") +
    scale_fill_gradient2(low = "steelblue", mid = "lightyellow", high = "firebrick",
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         labels = formatC(10^c(-1, -0.5, 0, 0.5, 1), digits = 2, format = "g"),
                         name = "Fold\nchange") +
    scale_x_continuous(name = xlabel, limits = x_range, expand = c(0, 0)) +
    scale_y_continuous(name = zlabel, limits = z_range, expand = c(0, 0),
                       breaks = z_breaks,
                       labels = formatC(10^z_breaks, format = "g")) +
    labs(title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_left <- ggplot(df, aes(x = .data[[zvar]])) +
    geom_histogram(bins = 30, fill = "grey40", color = NA) +
    scale_x_continuous(limits = z_range, expand = c(0, 0), breaks = z_breaks) +
    scale_y_reverse(expand = c(0, 0)) +
    coord_flip() +
    theme_void()

  p_bottom <- ggplot(df, aes(x = .data[[xvar]])) +
    geom_histogram(bins = 30, fill = "grey40", color = NA) +
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    theme_void()

  patchwork::wrap_plots(
    p_left, p_main,
    patchwork::plot_spacer(), p_bottom,
    ncol = 2, nrow = 2,
    widths  = c(1, 15),
    heights = c(15, 1)
  )
}

save_te_surface_single <- function(m, df, travel_var, travel_label, out_file) {
  if (is.null(m)) return(invisible(NULL))
  p <- plot_te_surface(m, df,
    xvar = "min_cbsa_dist", zvar = travel_var,
    xlabel = "CBSA distance (km)", zlabel = travel_label,
    title = sprintf("NB GAM: %s × distance", travel_label))
  ggsave(out_file, plot = p, width = 8, height = 6, units = "in", dpi = 300)
  message("Surface saved: ", out_file)
}

save_te_surface_rr_single <- function(m, df, travel_var, travel_label, out_file) {
  if (is.null(m)) return(invisible(NULL))
  p <- plot_te_surface(m, df,
    xvar = "min_cbsa_dist", zvar = travel_var,
    xlabel = "CBSA distance (km)", zlabel = travel_label,
    title = sprintf("RR GAM: %s × distance", travel_label),
    convert_log = FALSE)
  ggsave(out_file, plot = p, width = 8, height = 6, units = "in", dpi = 300)
  message("Surface saved: ", out_file)
}

# ============================================================================
# TENSOR PRODUCT GAMs — NB count form
# ============================================================================

gam_trips <- tryCatch(
  mgcv::gam(N_pairs ~ offset(log_offset) + te(log_RR_trips, min_cbsa_dist),
            family = mgcv::nb(), data = df_trips_data, method = "REML"),
  error = function(e) { message("GAM trips failed: ", conditionMessage(e)); NULL }
)
gam_t100 <- tryCatch(
  mgcv::gam(N_pairs ~ offset(log_offset) + te(log_RR_air_t100, min_cbsa_dist),
            family = mgcv::nb(), data = df_t100_data, method = "REML"),
  error = function(e) { message("GAM T100 failed: ", conditionMessage(e)); NULL }
)
gam_db1b <- tryCatch(
  mgcv::gam(N_pairs ~ offset(log_offset) + te(log_RR_air_db1b, min_cbsa_dist),
            family = mgcv::nb(), data = df_db1b_data, method = "REML"),
  error = function(e) { message("GAM DB1B failed: ", conditionMessage(e)); NULL }
)
gam_move <- tryCatch(
  mgcv::gam(N_pairs ~ offset(log_offset) + te(log_RR_move, min_cbsa_dist),
            family = mgcv::nb(), data = df_move_data, method = "REML"),
  error = function(e) { message("GAM SafeGraph failed: ", conditionMessage(e)); NULL }
)

# ============================================================================
# ti() DECOMPOSED GAMs — NB and RR (REML, for decomposition plots)
# ============================================================================

fit_ti <- function(outcome_formula, travel_var, df, family = gaussian(), label = "") {
  f <- as.formula(sprintf("%s ~ ti(%s) + ti(min_cbsa_dist) + ti(%s, min_cbsa_dist)",
                          outcome_formula, travel_var, travel_var))
  tryCatch(mgcv::gam(f, family = family, data = df, method = "REML"),
           error = function(e) { message(label, " ti failed: ", e$message); NULL })
}

ti_nb_trips <- fit_ti("N_pairs", "log_RR_trips",    df_trips_data, mgcv::nb(), "NB NHTS")
ti_nb_t100  <- fit_ti("N_pairs", "log_RR_air_t100", df_t100_data,  mgcv::nb(), "NB T100")
ti_nb_db1b  <- fit_ti("N_pairs", "log_RR_air_db1b", df_db1b_data,  mgcv::nb(), "NB DB1B")
ti_nb_move  <- fit_ti("N_pairs", "log_RR_move",     df_move_data,  mgcv::nb(), "NB SafeGraph")

ti_rr_trips <- fit_ti("log10(RR_seq)", "log_RR_trips",    df_trips_data, gaussian(), "RR NHTS")
ti_rr_t100  <- fit_ti("log10(RR_seq)", "log_RR_air_t100", df_t100_data,  gaussian(), "RR T100")
ti_rr_db1b  <- fit_ti("log10(RR_seq)", "log_RR_air_db1b", df_db1b_data,  gaussian(), "RR DB1B")
ti_rr_move  <- fit_ti("log10(RR_seq)", "log_RR_move",     df_move_data,  gaussian(), "RR SafeGraph")

# ============================================================================
# TENSOR PRODUCT GAMs — RR form
# ============================================================================

gam_rr_trips <- tryCatch(
  mgcv::gam(log10(RR_seq) ~ te(log_RR_trips, min_cbsa_dist),
            data = df_trips_data, method = "REML"),
  error = function(e) { message("GAM RR trips failed: ", conditionMessage(e)); NULL }
)
gam_rr_t100 <- tryCatch(
  mgcv::gam(log10(RR_seq) ~ te(log_RR_air_t100, min_cbsa_dist),
            data = df_t100_data, method = "REML"),
  error = function(e) { message("GAM RR T100 failed: ", conditionMessage(e)); NULL }
)
gam_rr_db1b <- tryCatch(
  mgcv::gam(log10(RR_seq) ~ te(log_RR_air_db1b, min_cbsa_dist),
            data = df_db1b_data, method = "REML"),
  error = function(e) { message("GAM RR DB1B failed: ", conditionMessage(e)); NULL }
)
gam_rr_move <- tryCatch(
  mgcv::gam(log10(RR_seq) ~ te(log_RR_move, min_cbsa_dist),
            data = df_move_data, method = "REML"),
  error = function(e) { message("GAM RR SafeGraph failed: ", conditionMessage(e)); NULL }
)

# ============================================================================
# RESULTS OUTPUT
# ============================================================================

write_gam_results(gam_trips, "NB: NHTS × distance")
write_gam_results(gam_t100,  "NB: T100 × distance")
write_gam_results(gam_db1b,  "NB: DB1B × distance")
write_gam_results(gam_move,  "NB: SafeGraph × distance")
write_gam_results(gam_rr_trips, "RR: NHTS × distance")
write_gam_results(gam_rr_t100,  "RR: T100 × distance")
write_gam_results(gam_rr_db1b,  "RR: DB1B × distance")
write_gam_results(gam_rr_move,  "RR: SafeGraph × distance")

# ============================================================================
# TENSOR PRODUCT SURFACES
# ============================================================================

save_te_surface_single(gam_trips, df_trips_data, "log_RR_trips",    "RR_NHTS",      paste0("figs/", scenario, "/dist/gam_count_trips_te.jpg"))
save_te_surface_single(gam_t100,  df_t100_data,  "log_RR_air_t100", "RR_T100",      paste0("figs/", scenario, "/dist/gam_count_t100_te.jpg"))
save_te_surface_single(gam_db1b,  df_db1b_data,  "log_RR_air_db1b", "RR_DB1B",      paste0("figs/", scenario, "/dist/gam_count_db1b_te.jpg"))
save_te_surface_single(gam_move,  df_move_data,  "log_RR_move",     "RR_SafeGraph", paste0("figs/", scenario, "/dist/gam_count_move_te.jpg"))

save_te_surface_rr_single(gam_rr_trips, df_trips_data, "log_RR_trips",    "RR_NHTS",      paste0("figs/", scenario, "/dist/gam_rr_trips_te.jpg"))
save_te_surface_rr_single(gam_rr_t100,  df_t100_data,  "log_RR_air_t100", "RR_T100",      paste0("figs/", scenario, "/dist/gam_rr_t100_te.jpg"))
save_te_surface_rr_single(gam_rr_db1b,  df_db1b_data,  "log_RR_air_db1b", "RR_DB1B",      paste0("figs/", scenario, "/dist/gam_rr_db1b_te.jpg"))
save_te_surface_rr_single(gam_rr_move,  df_move_data,  "log_RR_move",     "RR_SafeGraph", paste0("figs/", scenario, "/dist/gam_rr_move_te.jpg"))

# ============================================================================
# LRT: te(X,dist) vs ti(X) + ti(dist) — interaction significance test
# ============================================================================

run_ti_lrt <- function(var, label, df) {
  nb_fam <- mgcv::nb()
  f_add  <- as.formula(sprintf("N_pairs ~ offset(log_offset) + ti(%s) + ti(min_cbsa_dist)", var))
  f_full <- as.formula(sprintf("N_pairs ~ offset(log_offset) + ti(%s) + ti(min_cbsa_dist) + ti(%s, min_cbsa_dist)", var, var))
  m_add  <- tryCatch(mgcv::gam(f_add,  family = nb_fam, data = df, method = "ML"),
                     error = function(e) { message(label, " NB ti additive failed: ", e$message); NULL })
  m_full <- tryCatch(mgcv::gam(f_full, family = nb_fam, data = df, method = "ML"),
                     error = function(e) { message(label, " NB ti full failed: ", e$message); NULL })
  f_rr_add  <- as.formula(sprintf("log10(RR_seq) ~ ti(%s) + ti(min_cbsa_dist)", var))
  f_rr_full <- as.formula(sprintf("log10(RR_seq) ~ ti(%s) + ti(min_cbsa_dist) + ti(%s, min_cbsa_dist)", var, var))
  m_rr_add  <- tryCatch(mgcv::gam(f_rr_add,  data = df, method = "ML"),
                        error = function(e) { message(label, " RR ti additive failed: ", e$message); NULL })
  m_rr_full <- tryCatch(mgcv::gam(f_rr_full, data = df, method = "ML"),
                        error = function(e) { message(label, " RR ti full failed: ", e$message); NULL })
  sink(results_file, append = TRUE)
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat(paste0("LRT: interaction significance — ", label, "\n"))
  cat(paste0(strrep("=", 80), "\n\n"))
  if (!is.null(m_add) && !is.null(m_full)) {
    cat("NB count — ti(RR_travel) + ti(dist)  vs  ti(RR_travel) + ti(dist) + ti(RR_travel,dist):\n")
    print(anova(m_add, m_full, test = "Chisq"))
  }
  if (!is.null(m_rr_add) && !is.null(m_rr_full)) {
    cat("\nRR (Gaussian) — ti(RR_travel) + ti(dist)  vs  ti(RR_travel) + ti(dist) + ti(RR_travel,dist):\n")
    print(anova(m_rr_add, m_rr_full, test = "F"))
  }
  sink()
  message(sprintf("LRT written for %s", label))
}

run_ti_lrt("log_RR_trips",    "NHTS",      df_trips_data)
run_ti_lrt("log_RR_air_t100", "T100",      df_t100_data)
run_ti_lrt("log_RR_air_db1b", "DB1B",      df_db1b_data)
run_ti_lrt("log_RR_move",     "SafeGraph", df_move_data)

# ============================================================================
# DEVIANCE DECOMPOSITION
# ============================================================================

decompose_dev <- function(m_full, travel_var, df, label) {
  m_dist_only   <- mgcv::gam(N_pairs ~ offset(log_offset) + s(min_cbsa_dist),
                              family = mgcv::nb(), data = df, method = "REML")
  f_travel_only <- as.formula(sprintf("N_pairs ~ offset(log_offset) + s(%s)", travel_var))
  m_travel_only <- mgcv::gam(f_travel_only, family = mgcv::nb(), data = df, method = "REML")
  d2_full   <- summary(m_full)$dev.expl
  d2_dist   <- summary(m_dist_only)$dev.expl
  d2_travel <- summary(m_travel_only)$dev.expl
  shared        <- pmax(d2_dist + d2_travel - d2_full, 0)
  unique_dist   <- d2_full - shared - (d2_full - d2_dist)
  unique_travel <- d2_full - d2_dist
  data.frame(source = label, unique_dist = unique_dist, shared = shared,
             unique_travel = unique_travel, unexplained = 1 - d2_full,
             total_explained = d2_full)
}

df_dev_decomp <- bind_rows(
  decompose_dev(gam_trips, "log_RR_trips",    df_trips_data, "NHTS"),
  decompose_dev(gam_t100,  "log_RR_air_t100", df_t100_data,  "T100"),
  decompose_dev(gam_db1b,  "log_RR_air_db1b", df_db1b_data,  "DB1B"),
  decompose_dev(gam_move,  "log_RR_move",     df_move_data,  "SafeGraph")
) %>%
  mutate(
    shared      = pmax(shared, 0),
    unique_dist = total_explained - shared - unique_travel
  ) %>%
  pivot_longer(c(unique_dist, shared, unique_travel, unexplained),
               names_to = "component", values_to = "deviance") %>%
  mutate(
    component = factor(component,
      levels = c("unexplained", "unique_travel", "shared", "unique_dist"),
      labels = c("Unexplained", "Unique to travel", "Shared", "Unique to distance")
    ),
    source = factor(source, levels = c("NHTS", "T100", "DB1B", "SafeGraph"))
  )

p_dev_decomp <- ggplot(df_dev_decomp, aes(x = source, y = deviance * 100, fill = component)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(
    values = c("Unique to distance" = "steelblue", "Shared" = "grey70",
               "Unique to travel" = "firebrick", "Unexplained" = "grey92"),
    name = NULL
  ) +
  scale_y_continuous(name = "Deviance explained (%)", limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(name = "Travel data source") +
  labs(title = "NB GAM deviance decomposition: distance vs. travel variable") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")

ggsave(paste0("figs/", scenario, "/dist/nb_gam_dev_decomp.jpg"),
       plot = p_dev_decomp, width = 7, height = 5, units = "in", dpi = 300)
message("Deviance decomposition plot saved.")

# ============================================================================
# OBSERVED VS PREDICTED PLOTS
# ============================================================================

gam_trips_plot    <- if (!is.null(gam_trips))    plot_pred_nb_gam(gam_trips,    df_trips_data, "GAM NB NHTS (te)")      else patchwork::plot_spacer()
gam_t100_plot     <- if (!is.null(gam_t100))     plot_pred_nb_gam(gam_t100,     df_t100_data,  "GAM NB T100 (te)")      else patchwork::plot_spacer()
gam_db1b_plot     <- if (!is.null(gam_db1b))     plot_pred_nb_gam(gam_db1b,     df_db1b_data,  "GAM NB DB1B (te)")      else patchwork::plot_spacer()
gam_move_plot     <- if (!is.null(gam_move))     plot_pred_nb_gam(gam_move,     df_move_data,  "GAM NB SafeGraph (te)") else patchwork::plot_spacer()
gam_rr_trips_plot <- if (!is.null(gam_rr_trips)) plot_pred_rr_gam(gam_rr_trips, df_trips_data, "GAM RR NHTS (te)")      else patchwork::plot_spacer()
gam_rr_t100_plot  <- if (!is.null(gam_rr_t100))  plot_pred_rr_gam(gam_rr_t100,  df_t100_data,  "GAM RR T100 (te)")      else patchwork::plot_spacer()
gam_rr_db1b_plot  <- if (!is.null(gam_rr_db1b))  plot_pred_rr_gam(gam_rr_db1b,  df_db1b_data,  "GAM RR DB1B (te)")      else patchwork::plot_spacer()
gam_rr_move_plot  <- if (!is.null(gam_rr_move))  plot_pred_rr_gam(gam_rr_move,  df_move_data,  "GAM RR SafeGraph (te)") else patchwork::plot_spacer()

gam_plots <- patchwork::wrap_plots(
  gam_trips_plot,    gam_t100_plot,    gam_db1b_plot,    gam_move_plot,
  gam_rr_trips_plot, gam_rr_t100_plot, gam_rr_db1b_plot, gam_rr_move_plot,
  ncol = 4, nrow = 2
)
ggsave(paste0("figs/", scenario, "/dist/gam_fit.jpg"),
       plot = gam_plots, width = 22, height = 11, units = "in", dpi = 300)
message("GAM fit plots saved to figs/", scenario, "/dist/gam_fit.jpg")

# ============================================================================
# ti() DECOMPOSITION SURFACES
# ============================================================================

plot_ti_decomp <- function(m, df, travel_var, travel_label, convert_log = TRUE, n = 60) {
  if (is.null(m)) return(patchwork::plot_spacer())
  scale_factor <- ifelse(convert_log, log(10), 1)

  travel_seq <- seq(min(df[[travel_var]], na.rm = TRUE), max(df[[travel_var]], na.rm = TRUE), length.out = n)
  dist_seq   <- seq(min(df$min_cbsa_dist, na.rm = TRUE), max(df$min_cbsa_dist, na.rm = TRUE), length.out = n)

  fill_medians <- function(g) {
    for (col in setdiff(names(df), c(travel_var, "min_cbsa_dist"))) {
      if (is.numeric(df[[col]])) g[[col]] <- median(df[[col]], na.rm = TRUE)
      else g[[col]] <- df[[col]][1]
    }
    g
  }

  grid_travel <- fill_medians(data.frame(setNames(list(travel_seq), travel_var),
                                          min_cbsa_dist = median(df$min_cbsa_dist, na.rm = TRUE)))
  grid_dist   <- fill_medians(data.frame(setNames(list(rep(median(df[[travel_var]], na.rm = TRUE), n)), travel_var),
                                          min_cbsa_dist = dist_seq))
  grid_2d     <- fill_medians(expand.grid(setNames(list(travel_seq, dist_seq), c(travel_var, "min_cbsa_dist"))))

  terms_travel <- predict(m, newdata = grid_travel, type = "terms")
  terms_dist   <- predict(m, newdata = grid_dist,   type = "terms")
  terms_2d     <- predict(m, newdata = grid_2d,     type = "terms")

  pick_col <- function(mat, ...) {
    patterns <- c(...)
    for (p in patterns) {
      hit <- grep(p, colnames(mat), value = TRUE, perl = TRUE)
      if (length(hit)) return(hit[1])
    }
    colnames(mat)[1]
  }
  travel_col <- pick_col(terms_travel, paste0("^ti\\(", travel_var, "\\)$"), travel_var)
  dist_col   <- pick_col(terms_dist,   "^ti\\(min_cbsa_dist\\)$", "min_cbsa_dist")
  inter_col  <- pick_col(terms_2d,
                         paste0("(?=.*", travel_var, ")(?=.*min_cbsa_dist)"),
                         travel_var)

  travel_breaks <- seq(ceiling(min(travel_seq)), floor(max(travel_seq)))

  p_travel <- ggplot(data.frame(x = travel_seq,
                                 y = as.numeric(terms_travel[, travel_col]) / scale_factor),
                     aes(x = x, y = y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(color = "steelblue", linewidth = 1) +
    scale_x_continuous(name = travel_label, breaks = travel_breaks,
                       labels = formatC(10^travel_breaks, digits = 2, format = "g")) +
    scale_y_continuous(name = "Partial log10(RR)") +
    ggtitle("Marginal: travel") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_dist <- ggplot(data.frame(x = dist_seq,
                               y = as.numeric(terms_dist[, dist_col]) / scale_factor),
                   aes(x = x, y = y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(color = "steelblue", linewidth = 1) +
    scale_x_continuous(name = "CBSA distance (km)") +
    scale_y_continuous(name = "Partial log10(RR)") +
    ggtitle("Marginal: distance") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  grid_2d$partial <- as.numeric(terms_2d[, inter_col]) / scale_factor
  too_far <- mgcv::exclude.too.far(grid_2d[[travel_var]], grid_2d$min_cbsa_dist,
                                    df[[travel_var]],      df$min_cbsa_dist, dist = 0.3)
  grid_2d$partial[too_far] <- NA

  p_inter <- ggplot(grid_2d, aes(x = min_cbsa_dist, y = .data[[travel_var]], fill = partial)) +
    geom_tile() +
    geom_contour(aes(z = partial), color = "white", alpha = 0.5, linewidth = 0.3, na.rm = TRUE) +
    geom_point(data = df, aes(x = min_cbsa_dist, y = .data[[travel_var]]),
               inherit.aes = FALSE, size = 0.6, alpha = 0.2, color = "black") +
    scale_fill_gradient2(low = "steelblue", mid = "lightyellow", high = "firebrick",
                         midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         labels = formatC(10^c(-1, -0.5, 0, 0.5, 1), digits = 2, format = "g"),
                         name = "Fold\nchange", na.value = "grey90") +
    scale_x_continuous(name = "CBSA distance (km)", expand = c(0, 0)) +
    scale_y_continuous(name = travel_label, expand = c(0, 0),
                       breaks = travel_breaks,
                       labels = formatC(10^travel_breaks, digits = 2, format = "g")) +
    ggtitle("Interaction") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  patchwork::wrap_plots(p_travel, p_dist, p_inter, ncol = 3, widths = c(1, 1, 2))
}

save_ti_decomp <- function(m_nb, m_rr, df, travel_var, travel_label, out_file) {
  p_nb <- plot_ti_decomp(m_nb, df, travel_var, travel_label, convert_log = TRUE)
  p_rr <- plot_ti_decomp(m_rr, df, travel_var, travel_label, convert_log = FALSE)
  p <- patchwork::wrap_plots(p_nb, p_rr, ncol = 1) +
    patchwork::plot_annotation(title = sprintf("ti() decomposition: %s", travel_label),
                               theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  ggsave(out_file, plot = p, width = 14, height = 8, units = "in", dpi = 300)
  message("ti decomp saved: ", out_file)
}

save_ti_decomp(ti_nb_trips, ti_rr_trips, df_trips_data, "log_RR_trips",    "RR_NHTS",      paste0("figs/", scenario, "/dist/ti_decomp_trips.jpg"))
save_ti_decomp(ti_nb_t100,  ti_rr_t100,  df_t100_data,  "log_RR_air_t100", "RR_T100",      paste0("figs/", scenario, "/dist/ti_decomp_t100.jpg"))
save_ti_decomp(ti_nb_db1b,  ti_rr_db1b,  df_db1b_data,  "log_RR_air_db1b", "RR_DB1B",      paste0("figs/", scenario, "/dist/ti_decomp_db1b.jpg"))
save_ti_decomp(ti_nb_move,  ti_rr_move,  df_move_data,  "log_RR_move",     "RR_SafeGraph", paste0("figs/", scenario, "/dist/ti_decomp_move.jpg"))

# ============================================================================
# RR PREDICTION CURVES BY DISTANCE BAND
# ============================================================================

plot_rr_distance_curves <- function(m, df, travel_var, travel_label,
                                     distances = c(100, 500, 2000, 4000),
                                     dist_bins = list(
                                       "100"  = c(0,    250),
                                       "500"  = c(250,  1000),
                                       "2000" = c(1000, 3000),
                                       "4000" = c(3000, 5000)
                                     ),
                                     n = 200, custom_breaks = NULL) {
  global_range <- quantile(df[[travel_var]], c(0.02, 0.98), na.rm = TRUE)

  pred_list <- lapply(distances, function(d) {
    bin    <- dist_bins[[as.character(d)]]
    df_bin <- df[df$min_cbsa_dist >= bin[1] & df$min_cbsa_dist < bin[2], ]
    if (nrow(df_bin) < 5) df_bin <- df

    local_range <- quantile(df_bin[[travel_var]], c(0.02, 0.98), na.rm = TRUE)
    x_lo  <- max(local_range[1], global_range[1])
    x_hi  <- min(local_range[2], global_range[2])
    x_seq <- seq(x_lo, x_hi, length.out = n)

    sub <- data.frame(x_seq); names(sub) <- travel_var
    sub$min_cbsa_dist <- d
    for (col in setdiff(names(df), c(travel_var, "min_cbsa_dist"))) {
      if (is.numeric(df[[col]])) sub[[col]] <- median(df[[col]], na.rm = TRUE)
      else sub[[col]] <- df[[col]][1]
    }
    sub
  })
  pred_df <- do.call(rbind, pred_list)

  preds <- predict(m, newdata = pred_df, se.fit = TRUE)
  pred_df$fit       <- 10^preds$fit
  pred_df$fit_lo    <- 10^(preds$fit - 1.96 * preds$se.fit)
  pred_df$fit_hi    <- 10^(preds$fit + 1.96 * preds$se.fit)
  pred_df$x_rr      <- 10^pred_df[[travel_var]]
  pred_df$dist_label <- factor(pred_df$min_cbsa_dist, levels = distances,
                                labels = paste0(distances, " km"))

  x_breaks_log <- seq(ceiling(global_range[1]), floor(global_range[2]))
  x_breaks     <- if (!is.null(custom_breaks)) custom_breaks else 10^x_breaks_log

  ggplot(pred_df, aes(x = x_rr, y = fit, color = dist_label, fill = dist_label)) +
    geom_ribbon(aes(ymin = fit_lo, ymax = fit_hi), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_x_continuous(transform = "log10", name = travel_label,
                       breaks = x_breaks,
                       labels = formatC(x_breaks, digits = 2, format = "g")) +
    scale_y_continuous(transform = "log10", name = "Predicted Sequence RR",
                       breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
                       labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10")) +
    scale_color_viridis_d(name = "CBSA Distance", direction = -1, option = "plasma") +
    scale_fill_viridis_d(name = "CBSA Distance", direction = -1, option = "plasma") +
    theme_classic(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_curves_db1b <- plot_rr_distance_curves(gam_rr_db1b, df_db1b_data, "log_RR_air_db1b",
                                          "Air Travel RR (DB1B)",
                                          custom_breaks = c(0.1, 0.3, 0.5, 1, 2, 3))
ggsave(paste0("figs/", scenario, "/dist/pred_rr_curves_db1b.jpg"),
       plot = p_curves_db1b, width = 7, height = 5, units = "in", dpi = 300)
message("RR curves saved: pred_rr_curves_db1b.jpg")

p_curves_trips <- plot_rr_distance_curves(gam_rr_trips, df_trips_data, "log_RR_trips",
                                           "Ground Travel RR (NHTS)")
ggsave(paste0("figs/", scenario, "/dist/pred_rr_curves_trips.jpg"),
       plot = p_curves_trips, width = 7, height = 5, units = "in", dpi = 300)
message("RR curves saved: pred_rr_curves_trips.jpg")

# Stitch: NHTS (ground) left, DB1B (air) right, shared legend
p_trips_noleg <- p_curves_trips + theme(legend.position = "none")
p_db1b_noleg  <- p_curves_db1b  + theme(legend.position = "none")
legend_only   <- cowplot::get_legend(p_curves_db1b +
                   theme(legend.position = "right",
                         legend.title = element_text(size = 11),
                         legend.text  = element_text(size = 10)))
p_rr_curves_combined <- cowplot::plot_grid(
  p_trips_noleg, p_db1b_noleg, legend_only,
  nrow = 1, rel_widths = c(1, 1, 0.25)
)
ggsave(paste0("figs/", scenario, "/dist/pred_rr_curves_combined.jpg"),
       plot = p_rr_curves_combined, width = 10, height = 4, units = "in", dpi = 300)
ggsave(paste0("figs/", scenario, "/dist/pred_rr_curves_combined.svg"),
       plot = p_rr_curves_combined, width = 10, height = 4, units = "in")
message("RR curves combined saved: pred_rr_curves_combined.jpg / .svg")

p_curves_move <- plot_rr_distance_curves(gam_rr_move, df_move_data, "log_RR_move",
                                          "Mobility RR (SafeGraph)",
                                          custom_breaks = c(0.05, 0.1, 0.2, 0.5))
ggsave(paste0("figs/", scenario, "/dist/pred_rr_curves_move.jpg"),
       plot = p_curves_move, width = 10, height = 5, units = "in", dpi = 300)
ggsave(paste0("figs/", scenario, "/dist/pred_rr_curves_move.svg"),
       plot = p_curves_move, width = 10, height = 5, units = "in")
message("RR curves saved: pred_rr_curves_move.jpg / .svg")
