# Trace full nuclear and rank-truncated Frobenius errors at the recommended
# lambda factors for the p=30, q=25, n=500, cov-scale=0.2 diagnostic setting.

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

source_lines <- readLines(function_source, warn = FALSE)
cfg_start <- grep("^cfg <- list\\(", source_lines)
if (length(cfg_start) != 1L) stop("Cannot locate cfg block in ", function_source)
eval(parse(text = source_lines[seq_len(cfg_start - 1L)]), envir = .GlobalEnv)

nuclear_trace_with_rank_error <- function(xmat, y, p, q, rank0, sigma, M_truth,
                                          lambda_factor, maxit = 2000,
                                          tol = 1e-8) {
  m <- nrow(xmat)
  lambda_value <- lambda_factor * sigma * sqrt((p + q) / m)
  L <- estimate_lipschitz_ls(xmat)
  M_old <- matrix(0, p, q)
  Y <- M_old
  t_old <- 1
  obj_old <- nuclear_ls_objective(M_old, xmat, y, lambda_value)
  history <- vector("list", maxit)
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(maxit)) {
    resid <- as.vector(xmat %*% as.vector(Y) - y)
    grad <- matrix(crossprod(xmat, resid) / m, p, q)
    M_new <- svd_soft_threshold(Y - grad / L, lambda_value / L)
    M_rank <- rank_truncate(M_new, rank0)
    obj_new <- nuclear_ls_objective(M_new, xmat, y, lambda_value)

    rel_change <- abs(obj_old - obj_new) / (1 + abs(obj_old))
    step_frob <- mat_fnorm(M_new - M_old)
    rel_step_frob <- step_frob / (1 + mat_fnorm(M_old))
    svals <- svd(M_new, nu = 0, nv = 0)$d

    history[[iter]] <- data.frame(
      iter = iter,
      lambda_factor = lambda_factor,
      lambda = lambda_value,
      objective = obj_new,
      rel_objective_change = rel_change,
      step_frob = step_frob,
      rel_step_frob = rel_step_frob,
      full_truth_frob = mat_fnorm(M_new - M_truth),
      rank_truncated_truth_frob = mat_fnorm(M_rank - M_truth),
      rank_1e8 = sum(svals > 1e-8),
      rank_1e4 = sum(svals > 1e-4),
      converged = FALSE
    )

    if (is.finite(rel_change) && rel_change < tol) {
      converged <- TRUE
      history[[iter]]$converged <- TRUE
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

  history <- do.call(rbind, history[seq_len(iter)])
  list(
    M = M_old,
    lambda = lambda_value,
    lambda_factor = lambda_factor,
    lipschitz = L,
    iter = iter,
    converged = converged,
    objective = obj_old,
    history = history
  )
}

plot_paths <- function(history_df, cfg, summary_df, pdf_file, png_file) {
  draw <- function() {
    op <- par(mfrow = c(1, 3), mar = c(4.6, 4.8, 3.0, 1.0))
    on.exit(par(op), add = TRUE)

    full_factor <- cfg$full_lambda_factor
    rank_factor <- cfg$rank_lambda_factor
    selected <- c(full_factor, rank_factor)
    cols <- c("#2166AC", "#1B7837")

    for (lf in selected) {
      z <- history_df[abs(history_df$lambda_factor - lf) < 1e-12, ]
      ylim <- range(c(z$full_truth_frob, z$rank_truncated_truth_frob), finite = TRUE)
      plot(z$iter, z$full_truth_frob,
           type = "l", lwd = 2, col = cols[1],
           ylim = ylim, xlab = "iteration",
           ylab = expression("Frobenius error to " * M[0]),
           main = paste0("lambda factor = ", lf))
      lines(z$iter, z$rank_truncated_truth_frob, lwd = 2, col = cols[2])
      grid()
      abline(v = z$iter[nrow(z)], lty = 2, col = "#4D4D4D")
      legend("topright", bty = "n", lwd = 2, col = cols,
             legend = c("full nuclear", paste0("rank-", cfg$rank, " truncation")),
             cex = 0.85)
    }

    z_full <- history_df[abs(history_df$lambda_factor - full_factor) < 1e-12, ]
    z_rank <- history_df[abs(history_df$lambda_factor - rank_factor) < 1e-12, ]
    ylim <- range(c(z_full$full_truth_frob, z_rank$rank_truncated_truth_frob), finite = TRUE)
    plot(z_full$iter, z_full$full_truth_frob,
         type = "l", lwd = 2, col = cols[1],
         ylim = ylim, xlab = "iteration",
         ylab = expression("Frobenius error to " * M[0]),
         main = "Recommended trajectories")
    lines(z_rank$iter, z_rank$rank_truncated_truth_frob, lwd = 2, col = cols[2])
    grid()
    legend("topright", bty = "n", lwd = 2, col = cols,
           legend = c(
             paste0("full nuclear, factor=", full_factor),
             paste0("rank-", cfg$rank, ", factor=", rank_factor)
           ),
           cex = 0.8)
  }

  pdf(pdf_file, width = 13, height = 4.8)
  draw()
  dev.off()
  png(png_file, width = 2600, height = 960, res = 180)
  draw()
  dev.off()
}

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
  full_lambda_factor = as.numeric(get_arg("full-lambda-factor", "0.04")),
  rank_lambda_factor = as.numeric(get_arg("rank-lambda-factor", "0.03")),
  nuclear_maxit = as.integer(get_arg("nuclear-maxit", "2000")),
  nuclear_tol = as.numeric(get_arg("nuclear-tol", "1e-8")),
  out_root = normalizePath(
    get_arg("out-root", file.path("0608result", "p30_q25_covscale02_m500_nuclear_recommended_lambda_paths")),
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
if (!is.finite(cfg$full_lambda_factor) || cfg$full_lambda_factor <= 0) {
  stop("--full-lambda-factor must be positive.")
}
if (!is.finite(cfg$rank_lambda_factor) || cfg$rank_lambda_factor <= 0) {
  stop("--rank-lambda-factor must be positive.")
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

lambda_factors <- unique(c(cfg$full_lambda_factor, cfg$rank_lambda_factor))
fits <- vector("list", length(lambda_factors))
for (ii in seq_along(lambda_factors)) {
  message("Tracing lambda_factor=", lambda_factors[ii], " ...")
  fits[[ii]] <- nuclear_trace_with_rank_error(
    dat$xmat, dat$y,
    p = cfg$p, q = cfg$q, rank0 = cfg$rank,
    sigma = cfg$sigma, M_truth = truth$M0,
    lambda_factor = lambda_factors[ii],
    maxit = cfg$nuclear_maxit,
    tol = cfg$nuclear_tol
  )
}

history_df <- do.call(rbind, lapply(fits, `[[`, "history"))
summary_df <- do.call(
  rbind,
  lapply(fits, function(fit) {
    z <- fit$history[nrow(fit$history), ]
    data.frame(
      seed = cfg$seed,
      sample_size = cfg$sample_size,
      p = cfg$p,
      q = cfg$q,
      rank = cfg$rank,
      sigma = cfg$sigma,
      cov_scale = cfg$cov_scale,
      lambda_factor = fit$lambda_factor,
      lambda = fit$lambda,
      iter = fit$iter,
      converged = fit$converged,
      final_full_truth_frob = z$full_truth_frob,
      final_rank_truncated_truth_frob = z$rank_truncated_truth_frob,
      final_rank_1e4 = z$rank_1e4,
      final_rel_objective_change = z$rel_objective_change
    )
  })
)

history_file <- file.path(cfg$out_root, "nuclear_recommended_lambda_iteration_history.csv")
summary_file <- file.path(cfg$out_root, "nuclear_recommended_lambda_iteration_summary.csv")
config_file <- file.path(cfg$out_root, "config.txt")
fig_pdf <- file.path(cfg$out_root, "figures", "fig_nuclear_recommended_lambda_iteration_paths.pdf")
fig_png <- file.path(cfg$out_root, "figures", "fig_nuclear_recommended_lambda_iteration_paths.png")

write.csv(history_df, history_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)
writeLines(
  vapply(
    names(cfg),
    function(nm) paste(nm, paste(cfg[[nm]], collapse = ","), sep = " = "),
    character(1)
  ),
  config_file
)
plot_paths(history_df, cfg, summary_df, fig_pdf, fig_png)

message("Done.")
message("Summary: ", summary_file)
message("History: ", history_file)
message("Figure: ", fig_png)
print(summary_df)
