# Sweep nuclear-norm penalty strength for the p=30, q=25, n=500, cov-scale=0.2
# diagnostic setting.  This script does not source or modify the full simulation
# driver; it reuses function definitions from test_nuclear_norm_initial_convergence.R
# and generates one fixed data set for all lambda values.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[length(hit)])
}

parse_numeric_grid <- function(value) {
  out <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(out) == 0L || any(!is.finite(out))) stop("Invalid numeric grid: ", value)
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]] %||% NA_character_
script_file <- sub("^--file=", "", file_arg)
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
if (!nzchar(script_dir) || is.na(script_dir) || script_dir == ".") script_dir <- "R"
function_source <- file.path(script_dir, "test_nuclear_norm_initial_convergence.R")
if (!file.exists(function_source)) {
  function_source <- file.path("R", "test_nuclear_norm_initial_convergence.R")
}
if (!file.exists(function_source)) stop("Cannot find test_nuclear_norm_initial_convergence.R")

lines <- readLines(function_source, warn = FALSE)
cfg_start <- grep("^cfg <- list\\(", lines)
if (length(cfg_start) != 1L) stop("Cannot locate cfg block in ", function_source)
eval(parse(text = lines[seq_len(cfg_start - 1L)]), envir = .GlobalEnv)

cfg <- list(
  seed = as.integer(get_arg("seed", "20260527")),
  p = as.integer(get_arg("p", "30")),
  q = as.integer(get_arg("q", "25")),
  rank = as.integer(get_arg("rank", "2")),
  sigma = as.numeric(get_arg("sigma", "0.5")),
  signal = parse_numeric_grid(get_arg("signal", "6.0,4.4")),
  design = get_arg("design", "kronecker"),
  rho_p = as.numeric(get_arg("rho-p", "0.8")),
  rho_q = as.numeric(get_arg("rho-q", "0.7")),
  cov_scale = as.numeric(get_arg("cov-scale", "0.2")),
  sample_size = as.integer(get_arg("sample-size", get_arg("m", "500"))),
  nuclear_lambda_factors = parse_numeric_grid(
    get_arg(
      "nuclear-lambda-factors",
      "0.02,0.05,0.075,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5,2,3,4"
    )
  ),
  nuclear_maxit = as.integer(get_arg("nuclear-maxit", "2000")),
  nuclear_tol = as.numeric(get_arg("nuclear-tol", "1e-8")),
  scaledgd_reference_frob = as.numeric(get_arg("scaledgd-reference-frob", "3.91088544508124")),
  out_root = normalizePath(
    get_arg("out-root", file.path("0608result", "p30_q25_covscale02_m500_nuclear_lambda_sweep")),
    winslash = "/", mustWork = FALSE
  )
)

if (!cfg$design %in% c("isotropic", "kronecker")) stop("--design must be isotropic or kronecker.")
if (cfg$p < 1L || cfg$q < 1L) stop("p and q must be positive integers.")
if (cfg$rank < 1L || cfg$rank > min(cfg$p, cfg$q)) stop("Invalid rank.")
if (length(cfg$signal) != cfg$rank) stop("--signal length must equal --rank.")
if (!is.finite(cfg$sigma) || cfg$sigma < 0) stop("--sigma must be nonnegative.")
if (!is.finite(cfg$cov_scale) || cfg$cov_scale <= 0) stop("--cov-scale must be positive.")
if (cfg$sample_size < 2L) stop("--sample-size must be at least 2.")
if (any(!is.finite(cfg$nuclear_lambda_factors)) || any(cfg$nuclear_lambda_factors <= 0)) {
  stop("--nuclear-lambda-factors must contain positive finite values.")
}
if (cfg$nuclear_maxit < 1L) stop("--nuclear-maxit must be positive.")
if (!is.finite(cfg$nuclear_tol) || cfg$nuclear_tol <= 0) stop("--nuclear-tol must be positive.")

dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$out_root, "figures"), recursive = TRUE, showWarnings = FALSE)

message("Generating one fixed data set...")
set.seed(cfg$seed)
truth <- make_truth(cfg$p, cfg$q, cfg$rank, cfg$signal)
design <- make_design_objects(cfg)
dat <- generate_data(cfg$sample_size, truth, cfg, design$Sigmap, design$Sigmaq)

minimax_rate <- cfg$sigma * sqrt(cfg$rank * (cfg$p + cfg$q) / cfg$sample_size)
adjusted_minimax_rate <- minimax_rate / cfg$cov_scale
base_lambda <- cfg$sigma * sqrt((cfg$p + cfg$q) / cfg$sample_size)

summaries <- vector("list", length(cfg$nuclear_lambda_factors))
histories <- vector("list", length(cfg$nuclear_lambda_factors))

for (ii in seq_along(cfg$nuclear_lambda_factors)) {
  lambda_factor <- cfg$nuclear_lambda_factors[ii]
  message("Running lambda_factor=", lambda_factor, " ...")
  fit <- nuclear_norm_initial_trace(
    dat$xmat, dat$y, p = cfg$p, q = cfg$q, sigma = cfg$sigma,
    M_truth = truth$M0,
    lambda = "auto",
    lambda_factor = lambda_factor,
    maxit = cfg$nuclear_maxit,
    tol = cfg$nuclear_tol
  )

  history <- fit$history
  history$lambda_factor <- lambda_factor
  history$lambda <- fit$lambda
  histories[[ii]] <- history

  final <- history[nrow(history), , drop = FALSE]
  M_rank <- rank_truncate(fit$M, cfg$rank)
  summaries[[ii]] <- data.frame(
    seed = cfg$seed,
    sample_size = cfg$sample_size,
    p = cfg$p,
    q = cfg$q,
    rank = cfg$rank,
    sigma = cfg$sigma,
    signal = paste(cfg$signal, collapse = ","),
    design = cfg$design,
    cov_scale = cfg$cov_scale,
    nuclear_maxit = cfg$nuclear_maxit,
    nuclear_tol = cfg$nuclear_tol,
    lambda_factor = lambda_factor,
    lambda = fit$lambda,
    lambda_over_default = fit$lambda / base_lambda,
    lipschitz = fit$lipschitz,
    iter = fit$iter,
    converged = fit$converged,
    hit_maxit = fit$iter >= cfg$nuclear_maxit && !fit$converged,
    final_objective = fit$objective,
    final_rel_objective_change = final$rel_objective_change,
    final_step_frob = final$step_frob,
    final_rel_step_frob = final$rel_step_frob,
    final_truth_frob = mat_fnorm(fit$M - truth$M0),
    final_rank_truncated_truth_frob = mat_fnorm(M_rank - truth$M0),
    rank_1e8 = fit$rank_1e8,
    rank_1e4 = fit$rank_1e4,
    minimax_rate_standard = minimax_rate,
    ratio_standard = mat_fnorm(fit$M - truth$M0) / minimax_rate,
    rank_truncated_ratio_standard = mat_fnorm(M_rank - truth$M0) / minimax_rate,
    minimax_rate_covscale_adjusted = adjusted_minimax_rate,
    ratio_covscale_adjusted = mat_fnorm(fit$M - truth$M0) / adjusted_minimax_rate,
    rank_truncated_ratio_covscale_adjusted =
      mat_fnorm(M_rank - truth$M0) / adjusted_minimax_rate,
    scaledgd_reference_frob = cfg$scaledgd_reference_frob,
    ratio_to_scaledgd_reference = mat_fnorm(fit$M - truth$M0) / cfg$scaledgd_reference_frob,
    rank_truncated_ratio_to_scaledgd_reference =
      mat_fnorm(M_rank - truth$M0) / cfg$scaledgd_reference_frob
  )
}

summary_df <- do.call(rbind, summaries)
history_df <- do.call(rbind, histories)
summary_df <- summary_df[order(summary_df$lambda_factor), ]
history_df <- history_df[order(history_df$lambda_factor, history_df$iter), ]

best_full <- summary_df[which.min(summary_df$final_truth_frob), , drop = FALSE]
best_rank <- summary_df[which.min(summary_df$final_rank_truncated_truth_frob), , drop = FALSE]

summary_file <- file.path(cfg$out_root, "summary_nuclear_lambda_sweep.csv")
history_file <- file.path(cfg$out_root, "history_nuclear_lambda_sweep.csv")
config_file <- file.path(cfg$out_root, "config.txt")

write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(history_df, history_file, row.names = FALSE)
config_lines <- c(
  vapply(
    names(cfg),
    function(nm) paste(nm, paste(cfg[[nm]], collapse = ","), sep = " = "),
    character(1)
  ),
  paste("base_lambda =", base_lambda),
  paste("minimax_rate_standard =", minimax_rate),
  paste("minimax_rate_covscale_adjusted =", adjusted_minimax_rate),
  paste("best_full_lambda_factor =", best_full$lambda_factor),
  paste("best_full_truth_frob =", best_full$final_truth_frob),
  paste("best_rank_truncated_lambda_factor =", best_rank$lambda_factor),
  paste("best_rank_truncated_truth_frob =", best_rank$final_rank_truncated_truth_frob)
)
writeLines(config_lines, config_file)

plot_error <- function(pdf_file, png_file) {
  draw <- function() {
    op <- par(mar = c(4.8, 4.8, 3.2, 1.0))
    on.exit(par(op), add = TRUE)
    ylim <- range(
      c(summary_df$final_truth_frob, summary_df$final_rank_truncated_truth_frob,
        cfg$scaledgd_reference_frob),
      finite = TRUE
    )
    plot(summary_df$lambda_factor, summary_df$final_truth_frob,
         type = "b", log = "x", pch = 19, lwd = 2, col = "#2166AC",
         xlab = "nuclear lambda factor",
         ylab = expression("Frobenius error to " * M[0]),
         main = "Nuclear lambda sweep",
         ylim = ylim)
    lines(summary_df$lambda_factor, summary_df$final_rank_truncated_truth_frob,
          type = "b", pch = 17, lwd = 2, col = "#1B7837")
    abline(h = cfg$scaledgd_reference_frob, lty = 2, col = "#B2182B", lwd = 2)
    abline(v = 1, lty = 3, col = "#4D4D4D")
    abline(v = best_full$lambda_factor, lty = 4, col = "#2166AC")
    grid()
    legend("topright", bty = "n",
           lwd = c(2, 2, 2, 1), pch = c(19, 17, NA, NA),
           lty = c(1, 1, 2, 3),
           col = c("#2166AC", "#1B7837", "#B2182B", "#4D4D4D"),
           legend = c("full nuclear estimate", "rank-2 truncation",
                      "ScaledGD eta=1 reference", "default factor=1"))
  }
  pdf(pdf_file, width = 7.5, height = 5.5)
  draw()
  dev.off()
  png(png_file, width = 1500, height = 1100, res = 170)
  draw()
  dev.off()
}

plot_selected_histories <- function(pdf_file, png_file) {
  selected_raw <- unique(c(1, best_full$lambda_factor, best_rank$lambda_factor))
  selected <- selected_raw[
    vapply(
      selected_raw,
      function(x) any(abs(history_df$lambda_factor - x) < 1e-12),
      logical(1)
    )
  ]
  if (length(selected) == 0L) return(invisible(NULL))
  draw <- function() {
    op <- par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.0, 1.0))
    on.exit(par(op), add = TRUE)
    cols <- c("#4D4D4D", "#2166AC", "#1B7837")[seq_along(selected)]
    for (panel in seq_len(2)) {
      first <- TRUE
      for (ii in seq_along(selected)) {
        z <- history_df[abs(history_df$lambda_factor - selected[ii]) < 1e-12, ]
        y <- if (panel == 1L) z$truth_frob else pmax(z$rel_objective_change, .Machine$double.eps)
        if (first) {
          plot(z$iter, y, type = "l", lwd = 2, col = cols[ii],
               log = if (panel == 1L) "" else "y",
               xlab = "iteration",
               ylab = if (panel == 1L) expression("||" * M[k] - M[0] * "||"[F]) else "relative objective change",
               main = if (panel == 1L) "Truth Frobenius trajectory" else "Objective-change stopping metric")
          first <- FALSE
        } else {
          lines(z$iter, y, lwd = 2, col = cols[ii])
        }
      }
      if (panel == 1L) abline(h = cfg$scaledgd_reference_frob, lty = 2, col = "#B2182B")
      grid()
      legend("topright", bty = "n", lwd = 2, col = cols,
             legend = paste0("factor=", selected), cex = 0.85)
    }
  }
  pdf(pdf_file, width = 10, height = 4.8)
  draw()
  dev.off()
  png(png_file, width = 2000, height = 960, res = 170)
  draw()
  dev.off()
}

plot_error(
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_sweep_error.pdf"),
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_sweep_error.png")
)
plot_selected_histories(
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_sweep_selected_histories.pdf"),
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_sweep_selected_histories.png")
)

message("Done.")
message("Summary: ", summary_file)
message("History: ", history_file)
message("Best full nuclear factor: ", best_full$lambda_factor,
        ", F-error=", signif(best_full$final_truth_frob, 6))
message("Best rank-truncated factor: ", best_rank$lambda_factor,
        ", F-error=", signif(best_rank$final_rank_truncated_truth_frob, 6))
print(summary_df[, c(
  "lambda_factor", "lambda", "iter", "converged", "final_truth_frob",
  "final_rank_truncated_truth_frob", "rank_1e4", "ratio_covscale_adjusted",
  "rank_truncated_ratio_covscale_adjusted", "ratio_to_scaledgd_reference",
  "rank_truncated_ratio_to_scaledgd_reference"
)])
