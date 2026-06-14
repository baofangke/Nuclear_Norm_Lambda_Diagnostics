# Multi-seed nuclear-norm lambda sweep for p=30, q=25, n=500, cov-scale=0.2.
# The goal is to identify a lambda factor that performs well for most random
# data sets under this fixed simulation setting, rather than tuning to one seed.

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

parse_integer_grid <- function(value) {
  if (grepl(":", value, fixed = TRUE)) {
    endpoints <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    if (length(endpoints) != 2L || any(!is.finite(endpoints))) {
      stop("Invalid integer range: ", value)
    }
    return(seq(endpoints[1], endpoints[2]))
  }
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(out) == 0L || any(!is.finite(out))) stop("Invalid integer grid: ", value)
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

quant <- function(x, prob) {
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

cfg <- list(
  seeds = parse_integer_grid(get_arg("seeds", "20260501:20260530")),
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
      "0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.075,0.1,0.15,0.2"
    )
  ),
  nuclear_maxit = as.integer(get_arg("nuclear-maxit", "2000")),
  nuclear_tol = as.numeric(get_arg("nuclear-tol", "1e-8")),
  scaledgd_reference_frob = as.numeric(get_arg("scaledgd-reference-frob", "3.91088544508124")),
  out_root = normalizePath(
    get_arg("out-root", file.path("0608result", "p30_q25_covscale02_m500_nuclear_lambda_multiseed")),
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

minimax_rate <- cfg$sigma * sqrt(cfg$rank * (cfg$p + cfg$q) / cfg$sample_size)
adjusted_minimax_rate <- minimax_rate / cfg$cov_scale
base_lambda <- cfg$sigma * sqrt((cfg$p + cfg$q) / cfg$sample_size)

results <- vector("list", length(cfg$seeds) * length(cfg$nuclear_lambda_factors))
rr <- 0L

for (seed in cfg$seeds) {
  message("Seed ", seed, " ...")
  set.seed(seed)
  truth <- make_truth(cfg$p, cfg$q, cfg$rank, cfg$signal)
  design <- make_design_objects(cfg)
  dat <- generate_data(cfg$sample_size, truth, cfg, design$Sigmap, design$Sigmaq)
  L <- estimate_lipschitz_ls(dat$xmat)

  for (lambda_factor in cfg$nuclear_lambda_factors) {
    lambda <- lambda_factor * base_lambda
    fit <- nuclear_norm_initial_fast(
      dat$xmat, dat$y,
      p = cfg$p, q = cfg$q,
      lambda = lambda,
      lipschitz = L,
      maxit = cfg$nuclear_maxit,
      tol = cfg$nuclear_tol
    )
    M_rank <- rank_truncate(fit$M, cfg$rank)
    full_frob <- mat_fnorm(fit$M - truth$M0)
    rank_frob <- mat_fnorm(M_rank - truth$M0)
    rr <- rr + 1L
    results[[rr]] <- data.frame(
      seed = seed,
      sample_size = cfg$sample_size,
      p = cfg$p,
      q = cfg$q,
      rank = cfg$rank,
      sigma = cfg$sigma,
      cov_scale = cfg$cov_scale,
      lambda_factor = lambda_factor,
      lambda = lambda,
      lipschitz = L,
      iter = fit$iter,
      converged = fit$converged,
      hit_maxit = fit$iter >= cfg$nuclear_maxit && !fit$converged,
      final_objective = fit$objective,
      final_rel_objective_change = fit$final_rel_objective_change,
      final_step_frob = fit$final_step_frob,
      final_rel_step_frob = fit$final_rel_step_frob,
      final_truth_frob = full_frob,
      final_rank_truncated_truth_frob = rank_frob,
      rank_1e8 = fit$rank_1e8,
      rank_1e4 = fit$rank_1e4,
      minimax_rate_standard = minimax_rate,
      ratio_standard = full_frob / minimax_rate,
      rank_truncated_ratio_standard = rank_frob / minimax_rate,
      minimax_rate_covscale_adjusted = adjusted_minimax_rate,
      ratio_covscale_adjusted = full_frob / adjusted_minimax_rate,
      rank_truncated_ratio_covscale_adjusted = rank_frob / adjusted_minimax_rate,
      scaledgd_reference_frob = cfg$scaledgd_reference_frob,
      ratio_to_scaledgd_reference = full_frob / cfg$scaledgd_reference_frob,
      rank_truncated_ratio_to_scaledgd_reference = rank_frob / cfg$scaledgd_reference_frob
    )
  }
}

result_df <- do.call(rbind, results[seq_len(rr)])

seed_best_full <- stats::aggregate(
  final_truth_frob ~ seed, data = result_df, FUN = min
)
names(seed_best_full)[2] <- "seed_best_full_frob"
seed_best_rank <- stats::aggregate(
  final_rank_truncated_truth_frob ~ seed, data = result_df, FUN = min
)
names(seed_best_rank)[2] <- "seed_best_rank_frob"
result_df <- merge(result_df, seed_best_full, by = "seed")
result_df <- merge(result_df, seed_best_rank, by = "seed")
result_df$full_within_5pct_seed_best <-
  result_df$final_truth_frob <= 1.05 * result_df$seed_best_full_frob
result_df$full_within_10pct_seed_best <-
  result_df$final_truth_frob <= 1.10 * result_df$seed_best_full_frob
result_df$rank_within_5pct_seed_best <-
  result_df$final_rank_truncated_truth_frob <= 1.05 * result_df$seed_best_rank_frob
result_df$rank_within_10pct_seed_best <-
  result_df$final_rank_truncated_truth_frob <= 1.10 * result_df$seed_best_rank_frob

split_by_lambda <- split(result_df, result_df$lambda_factor)
aggregate_df <- do.call(
  rbind,
  lapply(split_by_lambda, function(z) {
    data.frame(
      lambda_factor = z$lambda_factor[1],
      lambda = z$lambda[1],
      n_seeds = length(unique(z$seed)),
      full_mean = mean(z$final_truth_frob),
      full_median = stats::median(z$final_truth_frob),
      full_sd = stats::sd(z$final_truth_frob),
      full_q25 = quant(z$final_truth_frob, 0.25),
      full_q75 = quant(z$final_truth_frob, 0.75),
      full_q90 = quant(z$final_truth_frob, 0.90),
      full_mean_adjusted_ratio = mean(z$ratio_covscale_adjusted),
      full_near_best_5pct = mean(z$full_within_5pct_seed_best),
      full_near_best_10pct = mean(z$full_within_10pct_seed_best),
      full_win_rate = mean(abs(z$final_truth_frob - z$seed_best_full_frob) < 1e-10),
      rank_mean = mean(z$final_rank_truncated_truth_frob),
      rank_median = stats::median(z$final_rank_truncated_truth_frob),
      rank_sd = stats::sd(z$final_rank_truncated_truth_frob),
      rank_q25 = quant(z$final_rank_truncated_truth_frob, 0.25),
      rank_q75 = quant(z$final_rank_truncated_truth_frob, 0.75),
      rank_q90 = quant(z$final_rank_truncated_truth_frob, 0.90),
      rank_mean_adjusted_ratio = mean(z$rank_truncated_ratio_covscale_adjusted),
      rank_near_best_5pct = mean(z$rank_within_5pct_seed_best),
      rank_near_best_10pct = mean(z$rank_within_10pct_seed_best),
      rank_win_rate = mean(abs(z$final_rank_truncated_truth_frob - z$seed_best_rank_frob) < 1e-10),
      mean_iter = mean(z$iter),
      convergence_rate = mean(z$converged),
      median_rank_1e4 = stats::median(z$rank_1e4)
    )
  })
)
aggregate_df <- aggregate_df[order(aggregate_df$lambda_factor), ]
rownames(aggregate_df) <- NULL

aggregate_df$full_robust_score <- aggregate_df$full_median + 0.5 * (aggregate_df$full_q75 - aggregate_df$full_median)
aggregate_df$rank_robust_score <- aggregate_df$rank_median + 0.5 * (aggregate_df$rank_q75 - aggregate_df$rank_median)

recommended_full <- aggregate_df[which.min(aggregate_df$full_robust_score), , drop = FALSE]
recommended_rank <- aggregate_df[which.min(aggregate_df$rank_robust_score), , drop = FALSE]
best_mean_full <- aggregate_df[which.min(aggregate_df$full_mean), , drop = FALSE]
best_mean_rank <- aggregate_df[which.min(aggregate_df$rank_mean), , drop = FALSE]

all_file <- file.path(cfg$out_root, "all_seed_lambda_results.csv")
agg_file <- file.path(cfg$out_root, "aggregate_nuclear_lambda_multiseed.csv")
config_file <- file.path(cfg$out_root, "config.txt")
recommend_file <- file.path(cfg$out_root, "recommendation_summary.txt")

write.csv(result_df[order(result_df$seed, result_df$lambda_factor), ], all_file, row.names = FALSE)
write.csv(aggregate_df, agg_file, row.names = FALSE)

config_lines <- c(
  vapply(
    names(cfg),
    function(nm) paste(nm, paste(cfg[[nm]], collapse = ","), sep = " = "),
    character(1)
  ),
  paste("base_lambda =", base_lambda),
  paste("minimax_rate_standard =", minimax_rate),
  paste("minimax_rate_covscale_adjusted =", adjusted_minimax_rate)
)
writeLines(config_lines, config_file)

recommend_lines <- c(
  "Recommended lambda factors under the fixed p=30,q=25,n=500,cov-scale=0.2 setting",
  paste("n_seeds =", length(cfg$seeds)),
  paste("base_lambda = sigma * sqrt((p+q)/n) =", signif(base_lambda, 8)),
  "",
  paste(
    "Full nuclear robust recommendation: lambda_factor =",
    recommended_full$lambda_factor,
    ", lambda =", signif(recommended_full$lambda, 8),
    ", median F-error =", signif(recommended_full$full_median, 6),
    ", q75 F-error =", signif(recommended_full$full_q75, 6),
    ", near-best-within-10pct =",
    signif(recommended_full$full_near_best_10pct, 4)
  ),
  paste(
    "Rank-2 truncation robust recommendation: lambda_factor =",
    recommended_rank$lambda_factor,
    ", lambda =", signif(recommended_rank$lambda, 8),
    ", median F-error =", signif(recommended_rank$rank_median, 6),
    ", q75 F-error =", signif(recommended_rank$rank_q75, 6),
    ", near-best-within-10pct =",
    signif(recommended_rank$rank_near_best_10pct, 4)
  ),
  "",
  paste(
    "Best mean full nuclear: lambda_factor =",
    best_mean_full$lambda_factor,
    ", mean F-error =", signif(best_mean_full$full_mean, 6)
  ),
  paste(
    "Best mean rank-2 truncation: lambda_factor =",
    best_mean_rank$lambda_factor,
    ", mean F-error =", signif(best_mean_rank$rank_mean, 6)
  )
)
writeLines(recommend_lines, recommend_file)

plot_error_bands <- function(pdf_file, png_file) {
  draw <- function() {
    op <- par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.2, 1.0))
    on.exit(par(op), add = TRUE)

    plot(aggregate_df$lambda_factor, aggregate_df$full_median,
         type = "b", log = "x", pch = 19, lwd = 2, col = "#2166AC",
         ylim = range(c(aggregate_df$full_q25, aggregate_df$full_q75), finite = TRUE),
         xlab = "nuclear lambda factor", ylab = "Frobenius error",
         main = "Full nuclear estimate")
    arrows(aggregate_df$lambda_factor, aggregate_df$full_q25,
           aggregate_df$lambda_factor, aggregate_df$full_q75,
           angle = 90, code = 3, length = 0.03, col = "#2166AC")
    abline(v = recommended_full$lambda_factor, lty = 2, col = "#B2182B")
    abline(h = cfg$scaledgd_reference_frob, lty = 3, col = "#B2182B")
    grid()
    legend("topright", bty = "n", lwd = c(2, 1, 1), lty = c(1, 2, 3),
           pch = c(19, NA, NA),
           col = c("#2166AC", "#B2182B", "#B2182B"),
           legend = c("median with IQR", "recommended", "ScaledGD reference"),
           cex = 0.85)

    plot(aggregate_df$lambda_factor, aggregate_df$rank_median,
         type = "b", log = "x", pch = 17, lwd = 2, col = "#1B7837",
         ylim = range(c(aggregate_df$rank_q25, aggregate_df$rank_q75), finite = TRUE),
         xlab = "nuclear lambda factor", ylab = "Frobenius error",
         main = "After rank-2 truncation")
    arrows(aggregate_df$lambda_factor, aggregate_df$rank_q25,
           aggregate_df$lambda_factor, aggregate_df$rank_q75,
           angle = 90, code = 3, length = 0.03, col = "#1B7837")
    abline(v = recommended_rank$lambda_factor, lty = 2, col = "#B2182B")
    abline(h = cfg$scaledgd_reference_frob, lty = 3, col = "#B2182B")
    grid()
    legend("topright", bty = "n", lwd = c(2, 1, 1), lty = c(1, 2, 3),
           pch = c(17, NA, NA),
           col = c("#1B7837", "#B2182B", "#B2182B"),
           legend = c("median with IQR", "recommended", "ScaledGD reference"),
           cex = 0.85)
  }
  pdf(pdf_file, width = 11, height = 5.2)
  draw()
  dev.off()
  png(png_file, width = 2200, height = 1040, res = 180)
  draw()
  dev.off()
}

plot_near_best <- function(pdf_file, png_file) {
  draw <- function() {
    op <- par(mar = c(4.8, 4.8, 3.2, 1.0))
    on.exit(par(op), add = TRUE)
    ylim <- c(0, 1)
    plot(aggregate_df$lambda_factor, aggregate_df$full_near_best_10pct,
         type = "b", log = "x", pch = 19, lwd = 2, col = "#2166AC",
         ylim = ylim, xlab = "nuclear lambda factor",
         ylab = "fraction within 10% of seed-best",
         main = "Near-best frequency across seeds")
    lines(aggregate_df$lambda_factor, aggregate_df$rank_near_best_10pct,
          type = "b", pch = 17, lwd = 2, col = "#1B7837")
    abline(v = recommended_full$lambda_factor, lty = 2, col = "#2166AC")
    abline(v = recommended_rank$lambda_factor, lty = 2, col = "#1B7837")
    grid()
    legend("topright", bty = "n", lwd = 2, pch = c(19, 17),
           col = c("#2166AC", "#1B7837"),
           legend = c("full nuclear", "rank-2 truncation"))
  }
  pdf(pdf_file, width = 7.5, height = 5.2)
  draw()
  dev.off()
  png(png_file, width = 1500, height = 1040, res = 180)
  draw()
  dev.off()
}

plot_error_bands(
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_multiseed_error_bands.pdf"),
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_multiseed_error_bands.png")
)
plot_near_best(
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_multiseed_near_best.pdf"),
  file.path(cfg$out_root, "figures", "fig_nuclear_lambda_multiseed_near_best.png")
)

message("Done.")
message("All results: ", all_file)
message("Aggregate: ", agg_file)
message("Recommendation: ", recommend_file)
message("Full nuclear robust factor: ", recommended_full$lambda_factor,
        ", lambda=", signif(recommended_full$lambda, 8))
message("Rank-truncated robust factor: ", recommended_rank$lambda_factor,
        ", lambda=", signif(recommended_rank$lambda, 8))
print(aggregate_df[, c(
  "lambda_factor", "lambda", "full_mean", "full_median", "full_q75",
  "full_near_best_10pct", "rank_mean", "rank_median", "rank_q75",
  "rank_near_best_10pct", "convergence_rate", "median_rank_1e4"
)])
