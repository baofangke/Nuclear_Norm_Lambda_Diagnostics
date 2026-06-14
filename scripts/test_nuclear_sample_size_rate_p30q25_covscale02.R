# Sample-size sweep for nuclear-norm estimators under p=30, q=25, cov-scale=0.2.
#
# This diagnostic reuses the recommended lambda factors from the n=500 sweep:
#   full nuclear estimator: lambda_factor = 0.04
#   rank-2 truncation:     lambda_factor = 0.03
#
# For each sample size n, lambda is scaled as
#   lambda = lambda_factor * sigma * sqrt((p + q) / n).
#
# Outputs:
#   nuclear_sample_size_rate_summary.csv
#   nuclear_sample_size_rate_fit.txt
#   config.txt
#   figures/fig_nuclear_sample_size_rate.pdf
#   figures/fig_nuclear_sample_size_rate.png

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[length(hit)])
}

parse_numeric_grid_local <- function(value) {
  out <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(out) == 0L || any(!is.finite(out))) stop("Invalid numeric grid: ", value)
  out
}

parse_integer_grid <- function(value) {
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(out) == 0L || any(!is.finite(out))) stop("Invalid integer grid: ", value)
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]] %||% NA_character_
script_file <- sub("^--file=", "", file_arg)
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
if (!nzchar(script_dir) || is.na(script_dir) || script_dir == ".") script_dir <- "scripts"
function_source <- file.path(script_dir, "test_nuclear_norm_initial_convergence.R")
if (!file.exists(function_source)) {
  function_source <- file.path("scripts", "test_nuclear_norm_initial_convergence.R")
}
if (!file.exists(function_source)) stop("Cannot find test_nuclear_norm_initial_convergence.R")

source_lines <- readLines(function_source, warn = FALSE)
cfg_start <- grep("^cfg <- list\\(", source_lines)
if (length(cfg_start) != 1L) stop("Cannot locate cfg block in ", function_source)
eval(parse(text = source_lines[seq_len(cfg_start - 1L)]), envir = .GlobalEnv)

nuclear_norm_initial_fast <- function(xmat, y, p, q, lambda, lipschitz,
                                      maxit = 2000, tol = 1e-8) {
  m <- nrow(xmat)
  L <- lipschitz
  M_old <- matrix(0, p, q)
  Y <- M_old
  t_old <- 1
  obj_old <- nuclear_ls_objective(M_old, xmat, y, lambda)
  converged <- FALSE
  iter <- 0L
  rel_change <- NA_real_
  step_frob <- NA_real_
  rel_step_frob <- NA_real_

  for (iter in seq_len(maxit)) {
    resid <- as.vector(xmat %*% as.vector(Y) - y)
    grad <- matrix(crossprod(xmat, resid) / m, p, q)
    M_new <- svd_soft_threshold(Y - grad / L, lambda / L)
    obj_new <- nuclear_ls_objective(M_new, xmat, y, lambda)

    rel_change <- abs(obj_old - obj_new) / (1 + abs(obj_old))
    step_frob <- mat_fnorm(M_new - M_old)
    rel_step_frob <- step_frob / (1 + mat_fnorm(M_old))

    if (is.finite(rel_change) && rel_change < tol) {
      converged <- TRUE
      M_old <- M_new
      obj_old <- obj_new
      break
    }

    t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
    Y <- M_new + ((t_old - 1) / t_new) * (M_new - M_old)
    M_old <- M_new
    t_old <- t_new
    obj_old <- obj_new
  }

  svals <- svd(M_old, nu = 0, nv = 0)$d
  list(
    M = M_old,
    iter = iter,
    converged = converged,
    objective = obj_old,
    final_rel_objective_change = rel_change,
    final_step_frob = step_frob,
    final_rel_step_frob = rel_step_frob,
    rank_1e8 = sum(svals > 1e-8),
    rank_1e4 = sum(svals > 1e-4),
    singular_values = svals
  )
}

plot_sample_size_rate <- function(summary_df, cfg, slopes, pdf_file, png_file) {
  draw <- function() {
    op <- par(mfrow = c(1, 2), mar = c(5.0, 5.0, 3.2, 1.0))
    on.exit(par(op), add = TRUE)

    n <- summary_df$sample_size
    full <- summary_df$full_truth_frob
    rank <- summary_df$rank_truncated_truth_frob
    anchor_idx <- which.min(abs(n - cfg$anchor_n))
    full_ref <- full[anchor_idx] * sqrt(n[anchor_idx] / n)
    rank_ref <- rank[anchor_idx] * sqrt(n[anchor_idx] / n)
    ylim <- range(c(full, rank, full_ref, rank_ref), finite = TRUE)

    plot(n, full, type = "b", pch = 19, lwd = 2, col = "#2166AC",
         ylim = ylim, xlab = "sample size n",
         ylab = expression("||" * hat(M) - M[0] * "||"[F]),
         main = "Nuclear Frobenius error")
    lines(n, rank, type = "b", pch = 17, lwd = 2, col = "#1B7837")
    lines(n, full_ref, lty = 2, lwd = 1.5, col = "#2166AC")
    lines(n, rank_ref, lty = 2, lwd = 1.5, col = "#1B7837")
    grid()
    legend("topright", bty = "n", lwd = c(2, 2, 1.5), lty = c(1, 1, 2),
           pch = c(19, 17, NA),
           col = c("#2166AC", "#1B7837", "#4D4D4D"),
           legend = c(
             paste0("full, factor=", cfg$full_lambda_factor),
             paste0("rank-", cfg$rank, ", factor=", cfg$rank_lambda_factor),
             paste0("C/sqrt(n), anchored at n=", n[anchor_idx])
           ),
           cex = 0.78)

    plot(n, full, type = "b", log = "xy", pch = 19, lwd = 2, col = "#2166AC",
         ylim = ylim, xlab = "sample size n (log)",
         ylab = expression("||" * hat(M) - M[0] * "||"[F] * " (log)"),
         main = sprintf("Log-log slopes: full %.3f, rank %.3f",
                        slopes$full_slope, slopes$rank_slope))
    lines(n, rank, type = "b", pch = 17, lwd = 2, col = "#1B7837")
    lines(n, full_ref, lty = 2, lwd = 1.5, col = "#2166AC")
    lines(n, rank_ref, lty = 2, lwd = 1.5, col = "#1B7837")
    grid()
    legend("topright", bty = "n", lwd = c(2, 2, 1.5), lty = c(1, 1, 2),
           pch = c(19, 17, NA),
           col = c("#2166AC", "#1B7837", "#4D4D4D"),
           legend = c("full nuclear", paste0("rank-", cfg$rank, " truncation"),
                      "slope -1/2 reference"),
           cex = 0.78)
  }

  pdf(pdf_file, width = 11, height = 5.2)
  draw()
  dev.off()

  png(png_file, width = 2200, height = 1040, res = 180)
  draw()
  dev.off()
}

cfg <- list(
  seed = as.integer(get_arg("seed", "20260527")),
  p = as.integer(get_arg("p", "30")),
  q = as.integer(get_arg("q", "25")),
  rank = as.integer(get_arg("rank", "2")),
  sigma = as.numeric(get_arg("sigma", "0.5")),
  signal = parse_numeric_grid_local(get_arg("signal", "6.0,4.4")),
  design = get_arg("design", "kronecker"),
  rho_p = as.numeric(get_arg("rho-p", "0.8")),
  rho_q = as.numeric(get_arg("rho-q", "0.7")),
  cov_scale = as.numeric(get_arg("cov-scale", "0.2")),
  sample_sizes = parse_integer_grid(get_arg("sample-sizes", "300,400,500,700,1000,1500,2000")),
  full_lambda_factor = as.numeric(get_arg("full-lambda-factor", "0.04")),
  rank_lambda_factor = as.numeric(get_arg("rank-lambda-factor", "0.03")),
  nuclear_maxit = as.integer(get_arg("nuclear-maxit", "2000")),
  nuclear_tol = as.numeric(get_arg("nuclear-tol", "1e-8")),
  anchor_n = as.integer(get_arg("anchor-n", "500")),
  out_root = normalizePath(
    get_arg("out-root", file.path("results", "nuclear_sample_size_rate")),
    winslash = "/", mustWork = FALSE
  )
)

if (!cfg$design %in% c("isotropic", "kronecker")) stop("--design must be isotropic or kronecker.")
if (cfg$p < 1L || cfg$q < 1L) stop("p and q must be positive integers.")
if (cfg$rank < 1L || cfg$rank > min(cfg$p, cfg$q)) stop("Invalid rank.")
if (length(cfg$signal) != cfg$rank) stop("--signal length must equal --rank.")
if (!is.finite(cfg$sigma) || cfg$sigma < 0) stop("--sigma must be nonnegative.")
if (!is.finite(cfg$cov_scale) || cfg$cov_scale <= 0) stop("--cov-scale must be positive.")
if (length(cfg$sample_sizes) < 2L || any(cfg$sample_sizes < 2L)) {
  stop("--sample-sizes must contain at least two integers >= 2.")
}
if (!is.finite(cfg$full_lambda_factor) || cfg$full_lambda_factor <= 0) {
  stop("--full-lambda-factor must be positive.")
}
if (!is.finite(cfg$rank_lambda_factor) || cfg$rank_lambda_factor <= 0) {
  stop("--rank-lambda-factor must be positive.")
}
if (cfg$nuclear_maxit < 1L) stop("--nuclear-maxit must be positive.")
if (!is.finite(cfg$nuclear_tol) || cfg$nuclear_tol <= 0) stop("--nuclear-tol must be positive.")

cfg$sample_sizes <- sort(unique(cfg$sample_sizes))
max_n <- max(cfg$sample_sizes)

dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$out_root, "figures"), recursive = TRUE, showWarnings = FALSE)

message("Generating one data set with max_n=", max_n, " ...")
set.seed(cfg$seed)
truth <- make_truth(cfg$p, cfg$q, cfg$rank, cfg$signal)
design <- make_design_objects(cfg)
dat_full <- generate_data(max_n, truth, cfg, design$Sigmap, design$Sigmaq)

results <- vector("list", length(cfg$sample_sizes))

for (ii in seq_along(cfg$sample_sizes)) {
  n <- cfg$sample_sizes[ii]
  message("Running n=", n, " ...")
  xmat <- dat_full$xmat[seq_len(n), , drop = FALSE]
  y <- dat_full$y[seq_len(n)]
  lipschitz <- estimate_lipschitz_ls(xmat)
  minimax_rate <- cfg$sigma * sqrt(cfg$rank * (cfg$p + cfg$q) / n)
  adjusted_minimax_rate <- minimax_rate / cfg$cov_scale
  base_lambda <- cfg$sigma * sqrt((cfg$p + cfg$q) / n)

  full_lambda <- cfg$full_lambda_factor * base_lambda
  full_fit <- nuclear_norm_initial_fast(
    xmat, y, p = cfg$p, q = cfg$q,
    lambda = full_lambda, lipschitz = lipschitz,
    maxit = cfg$nuclear_maxit, tol = cfg$nuclear_tol
  )
  full_frob <- mat_fnorm(full_fit$M - truth$M0)

  rank_lambda <- cfg$rank_lambda_factor * base_lambda
  rank_fit <- nuclear_norm_initial_fast(
    xmat, y, p = cfg$p, q = cfg$q,
    lambda = rank_lambda, lipschitz = lipschitz,
    maxit = cfg$nuclear_maxit, tol = cfg$nuclear_tol
  )
  M_rank <- rank_truncate(rank_fit$M, cfg$rank)
  raw_rank_factor_frob <- mat_fnorm(rank_fit$M - truth$M0)
  rank_frob <- mat_fnorm(M_rank - truth$M0)

  results[[ii]] <- data.frame(
    seed = cfg$seed,
    sample_size = n,
    p = cfg$p,
    q = cfg$q,
    rank = cfg$rank,
    sigma = cfg$sigma,
    signal = paste(cfg$signal, collapse = ","),
    design = cfg$design,
    cov_scale = cfg$cov_scale,
    full_lambda_factor = cfg$full_lambda_factor,
    full_lambda = full_lambda,
    full_iter = full_fit$iter,
    full_converged = full_fit$converged,
    full_truth_frob = full_frob,
    full_adjusted_minimax_ratio = full_frob / adjusted_minimax_rate,
    full_standard_minimax_ratio = full_frob / minimax_rate,
    full_rank_1e4 = full_fit$rank_1e4,
    full_final_rel_objective_change = full_fit$final_rel_objective_change,
    rank_lambda_factor = cfg$rank_lambda_factor,
    rank_lambda = rank_lambda,
    rank_iter = rank_fit$iter,
    rank_converged = rank_fit$converged,
    rank_factor_raw_truth_frob = raw_rank_factor_frob,
    rank_truncated_truth_frob = rank_frob,
    rank_adjusted_minimax_ratio = rank_frob / adjusted_minimax_rate,
    rank_standard_minimax_ratio = rank_frob / minimax_rate,
    rank_factor_rank_1e4 = rank_fit$rank_1e4,
    rank_final_rel_objective_change = rank_fit$final_rel_objective_change,
    minimax_rate_standard = minimax_rate,
    minimax_rate_covscale_adjusted = adjusted_minimax_rate,
    lipschitz = lipschitz
  )
}

summary_df <- do.call(rbind, results)
full_lm <- stats::lm(log(full_truth_frob) ~ log(sample_size), data = summary_df)
rank_lm <- stats::lm(log(rank_truncated_truth_frob) ~ log(sample_size), data = summary_df)
slopes <- list(
  full_slope = unname(stats::coef(full_lm)[2]),
  rank_slope = unname(stats::coef(rank_lm)[2]),
  full_intercept = unname(stats::coef(full_lm)[1]),
  rank_intercept = unname(stats::coef(rank_lm)[1])
)

summary_file <- file.path(cfg$out_root, "nuclear_sample_size_rate_summary.csv")
fit_file <- file.path(cfg$out_root, "nuclear_sample_size_rate_fit.txt")
config_file <- file.path(cfg$out_root, "config.txt")
fig_pdf <- file.path(cfg$out_root, "figures", "fig_nuclear_sample_size_rate.pdf")
fig_png <- file.path(cfg$out_root, "figures", "fig_nuclear_sample_size_rate.png")

write.csv(summary_df, summary_file, row.names = FALSE)

fit_lines <- c(
  "Nuclear-norm sample-size rate diagnostic",
  paste("sample_sizes =", paste(cfg$sample_sizes, collapse = ",")),
  paste("full_lambda_factor =", cfg$full_lambda_factor),
  paste("rank_lambda_factor =", cfg$rank_lambda_factor),
  paste("full_loglog_slope =", signif(slopes$full_slope, 8)),
  paste("rank_truncated_loglog_slope =", signif(slopes$rank_slope, 8)),
  "reference_slope = -0.5",
  "",
  "Interpretation: factors are fixed from the n=500 lambda diagnostic and are not retuned for each sample size."
)
writeLines(fit_lines, fit_file)

config_lines <- vapply(
  names(cfg),
  function(nm) paste(nm, paste(cfg[[nm]], collapse = ","), sep = " = "),
  character(1)
)
writeLines(config_lines, config_file)

plot_sample_size_rate(summary_df, cfg, slopes, fig_pdf, fig_png)

message("Done.")
message("Summary: ", summary_file)
message("Fit: ", fit_file)
message("Figure: ", fig_png)
print(summary_df[, c(
  "sample_size", "full_iter", "full_truth_frob", "full_adjusted_minimax_ratio",
  "rank_iter", "rank_truncated_truth_frob", "rank_adjusted_minimax_ratio"
)])
message("Full log-log slope: ", signif(slopes$full_slope, 6))
message("Rank-truncated log-log slope: ", signif(slopes$rank_slope, 6))
