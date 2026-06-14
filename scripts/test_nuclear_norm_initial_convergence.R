# Trace nuclear_norm_initial convergence for the default maxit=300 setting.
#
# This standalone smoke test mirrors the proximal-gradient/FISTA loop used by
# nuclear_norm_initial() in
# simulation2_rank2_signal64_single_target3_scaledgd_multi_methods_unifdbsp_uniformcache_conditional_parallel_linear.R.
# It records per-iteration Frobenius diagnostics without sourcing the main
# simulation file, because that file runs the full simulation at load time.
#
# Outputs:
#   nuclear_norm_initial_convergence_history.csv
#   nuclear_norm_initial_convergence_summary.csv
#   config.txt
#   figures/fig_nuclear_norm_initial_convergence.pdf
#   figures/fig_nuclear_norm_initial_convergence.png
#
# Example:
#   Rscript R/test_nuclear_norm_initial_convergence.R
#   Rscript R/test_nuclear_norm_initial_convergence.R --sample-size=800 --nuclear-maxit=300

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

mat_fnorm <- function(A) sqrt(sum(A * A))

ar1_cov <- function(d, rho) {
  idx <- seq_len(d)
  outer(idx, idx, function(i, j) rho^abs(i - j))
}

exchangeable_cov <- function(d, rho) {
  (1 - rho) * diag(d) + rho * matrix(1, d, d)
}

rand_orth <- function(n, k) {
  qr.Q(qr(matrix(rnorm(n * k), n, k)))[, seq_len(k), drop = FALSE]
}

make_truth <- function(p, q, rank0, signal) {
  U <- rand_orth(p, rank0)
  V <- rand_orth(q, rank0)
  if (length(signal) < rank0) {
    signal <- seq(signal[1], signal[1] * 0.7, length.out = rank0)
  }
  D <- diag(signal[seq_len(rank0)], nrow = rank0)
  list(M0 = U %*% D %*% t(V), U0 = U, V0 = V, D0 = D)
}

make_design_objects <- function(cfg) {
  Sigmap_raw <- if (cfg$design == "isotropic") diag(cfg$p) else ar1_cov(cfg$p, cfg$rho_p)
  Sigmaq_raw <- if (cfg$design == "isotropic") diag(cfg$q) else exchangeable_cov(cfg$q, cfg$rho_q)
  list(
    Sigmap = cfg$cov_scale * Sigmap_raw,
    Sigmaq = cfg$cov_scale * Sigmaq_raw,
    Sigmap_raw = Sigmap_raw,
    Sigmaq_raw = Sigmaq_raw
  )
}

make_xmat <- function(n, p, q, Sigmap, Sigmaq) {
  z <- matrix(rnorm(n * p * q), n, p * q)
  K <- kronecker(chol(Sigmaq), chol(Sigmap))
  z %*% K
}

generate_data <- function(n, truth, cfg, Sigmap, Sigmaq) {
  beta0 <- as.vector(truth$M0)
  xmat <- make_xmat(n, cfg$p, cfg$q, Sigmap, Sigmaq)
  eps <- cfg$sigma * rnorm(n)
  y <- as.vector(xmat %*% beta0 + eps)
  list(xmat = xmat, y = y, eps = eps)
}

svd_soft_threshold <- function(M, tau) {
  S <- svd(M)
  d <- pmax(S$d - tau, 0)
  keep <- which(d > 0)
  if (length(keep) == 0L) return(matrix(0, nrow(M), ncol(M)))
  S$u[, keep, drop = FALSE] %*% diag(d[keep], nrow = length(keep)) %*%
    t(S$v[, keep, drop = FALSE])
}

nuclear_norm_value <- function(M) {
  sum(svd(M, nu = 0, nv = 0)$d)
}

nuclear_ls_objective <- function(M, xmat, y, lambda) {
  resid <- as.vector(xmat %*% as.vector(M) - y)
  0.5 * mean(resid^2) + lambda * nuclear_norm_value(M)
}

estimate_lipschitz_ls <- function(xmat) {
  G <- crossprod(xmat) / nrow(xmat)
  L <- max(eigen((G + t(G)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(L) || L <= 0) L <- 1
  L
}

rank_truncate <- function(M, rank0) {
  S <- svd(M, nu = rank0, nv = rank0)
  U <- S$u[, seq_len(rank0), drop = FALSE]
  V <- S$v[, seq_len(rank0), drop = FALSE]
  D <- diag(S$d[seq_len(rank0)], nrow = rank0)
  U %*% D %*% t(V)
}

nuclear_norm_initial_trace <- function(xmat, y, p, q, sigma, M_truth = NULL,
                                       lambda = "auto", lambda_factor = 1.0,
                                       maxit = 300, tol = 1e-6) {
  m <- nrow(xmat)
  if (lambda == "auto") {
    lambda_value <- lambda_factor * sigma * sqrt((p + q) / m)
  } else {
    lambda_value <- as.numeric(lambda)
    if (!is.finite(lambda_value) || lambda_value < 0) {
      stop("--nuclear-lambda must be 'auto' or a nonnegative number.")
    }
  }

  L <- estimate_lipschitz_ls(xmat)
  M_old <- matrix(0, p, q)
  Y <- M_old
  t_old <- 1
  obj_old <- nuclear_ls_objective(M_old, xmat, y, lambda_value)
  truth_norm <- if (is.null(M_truth)) NA_real_ else max(mat_fnorm(M_truth), .Machine$double.eps)
  history <- vector("list", maxit)
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(maxit)) {
    resid <- as.vector(xmat %*% as.vector(Y) - y)
    grad <- matrix(crossprod(xmat, resid) / m, p, q)
    M_new <- svd_soft_threshold(Y - grad / L, lambda_value / L)
    obj_new <- nuclear_ls_objective(M_new, xmat, y, lambda_value)

    rel_change <- abs(obj_old - obj_new) / (1 + abs(obj_old))
    step_frob <- mat_fnorm(M_new - M_old)
    rel_step_frob <- step_frob / (1 + mat_fnorm(M_old))
    truth_frob <- if (is.null(M_truth)) NA_real_ else mat_fnorm(M_new - M_truth)
    rel_truth_frob <- if (is.null(M_truth)) NA_real_ else truth_frob / truth_norm
    truth_mse <- if (is.null(M_truth)) NA_real_ else mean((M_new - M_truth)^2)
    svals <- svd(M_new, nu = 0, nv = 0)$d

    history[[iter]] <- data.frame(
      iter = iter,
      objective = obj_new,
      rel_objective_change = rel_change,
      step_frob = step_frob,
      rel_step_frob = rel_step_frob,
      truth_frob = truth_frob,
      rel_truth_frob = rel_truth_frob,
      truth_matrix_mse = truth_mse,
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
  svals <- svd(M_old, nu = 0, nv = 0)$d
  list(
    M = M_old,
    lambda = lambda_value,
    lipschitz = L,
    iter = iter,
    converged = converged,
    objective = obj_old,
    history = history,
    rank_1e8 = sum(svals > 1e-8),
    rank_1e4 = sum(svals > 1e-4),
    singular_values = svals
  )
}

plot_convergence <- function(history, cfg, fit, pdf_file, png_file) {
  draw <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.0, 1.0))
    on.exit(par(op), add = TRUE)

    main_suffix <- paste0("m=", cfg$sample_size, ", maxit=", cfg$nuclear_maxit)
    stop_iter <- if (fit$converged) fit$iter else NA_integer_

    plot(history$iter, pmax(history$step_frob, .Machine$double.eps),
         type = "l", log = "y", lwd = 2, col = "#2166AC",
         xlab = "iteration", ylab = expression("||" * M[k] - M[k-1] * "||"[F]),
         main = paste("Frobenius step,", main_suffix))
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(history$iter, history$truth_frob,
         type = "l", lwd = 2, col = "#1B7837",
         xlab = "iteration", ylab = expression("||" * M[k] - M[0] * "||"[F]),
         main = "Frobenius error to truth")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(history$iter, history$objective,
         type = "l", lwd = 2, col = "#762A83",
         xlab = "iteration", ylab = "objective",
         main = "Penalized objective")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(history$iter, pmax(history$rel_objective_change, .Machine$double.eps),
         type = "l", log = "y", lwd = 2, col = "#B35806",
         xlab = "iteration", ylab = "relative objective change",
         main = "Original stopping criterion")
    abline(h = cfg$nuclear_tol, lty = 3, col = "#B2182B")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")
  }

  pdf(pdf_file, width = 9, height = 7)
  draw()
  dev.off()

  png(png_file, width = 1400, height = 1050, res = 150)
  draw()
  dev.off()
}

cfg <- list(
  seed = as.integer(get_arg("seed", "20260527")),
  p = as.integer(get_arg("p", "10")),
  q = as.integer(get_arg("q", "10")),
  rank = as.integer(get_arg("rank", "2")),
  sigma = as.numeric(get_arg("sigma", "0.5")),
  signal = parse_numeric_grid(get_arg("signal", "6.0,4.4")),
  design = get_arg("design", "kronecker"),
  rho_p = as.numeric(get_arg("rho-p", "0.8")),
  rho_q = as.numeric(get_arg("rho-q", "0.7")),
  cov_scale = as.numeric(get_arg("cov-scale", "1.0")),
  sample_size = as.integer(get_arg("sample-size", get_arg("m", "200"))),
  nuclear_lambda = get_arg("nuclear-lambda", "auto"),
  nuclear_lambda_factor = as.numeric(get_arg("nuclear-lambda-factor", "1.0")),
  nuclear_maxit = as.integer(get_arg("nuclear-maxit", "300")),
  nuclear_tol = as.numeric(get_arg("nuclear-tol", "1e-6")),
  reference_maxit = as.integer(get_arg("reference-maxit", "1000")),
  reference_tol = as.numeric(get_arg("reference-tol", "0")),
  out_root = normalizePath(
    get_arg("out-root", file.path("0608result", "nuclear_norm_initial_convergence_test")),
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
if (!is.finite(cfg$nuclear_lambda_factor) || cfg$nuclear_lambda_factor <= 0) {
  stop("--nuclear-lambda-factor must be positive.")
}
if (cfg$nuclear_maxit < 1L) stop("--nuclear-maxit must be positive.")
if (!is.finite(cfg$nuclear_tol) || cfg$nuclear_tol <= 0) stop("--nuclear-tol must be positive.")
if (cfg$reference_maxit < cfg$nuclear_maxit) {
  stop("--reference-maxit must be at least --nuclear-maxit.")
}
if (!is.finite(cfg$reference_tol) || cfg$reference_tol < 0) {
  stop("--reference-tol must be nonnegative.")
}

dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$out_root, "figures"), recursive = TRUE, showWarnings = FALSE)

message("Generating one test data set...")
set.seed(cfg$seed)
truth <- make_truth(cfg$p, cfg$q, cfg$rank, cfg$signal)
design <- make_design_objects(cfg)
dat <- generate_data(cfg$sample_size, truth, cfg, design$Sigmap, design$Sigmaq)

message("Tracing nuclear_norm_initial with maxit=", cfg$nuclear_maxit, "...")
fit <- nuclear_norm_initial_trace(
  dat$xmat, dat$y, p = cfg$p, q = cfg$q, sigma = cfg$sigma,
  M_truth = truth$M0,
  lambda = cfg$nuclear_lambda,
  lambda_factor = cfg$nuclear_lambda_factor,
  maxit = cfg$nuclear_maxit,
  tol = cfg$nuclear_tol
)

ref_fit <- NULL
if (cfg$reference_maxit > cfg$nuclear_maxit) {
  message("Computing reference run with maxit=", cfg$reference_maxit,
          " and tol=", cfg$reference_tol, "...")
  ref_fit <- nuclear_norm_initial_trace(
    dat$xmat, dat$y, p = cfg$p, q = cfg$q, sigma = cfg$sigma,
    M_truth = truth$M0,
    lambda = cfg$nuclear_lambda,
    lambda_factor = cfg$nuclear_lambda_factor,
    maxit = cfg$reference_maxit,
    tol = cfg$reference_tol
  )
}

history <- fit$history
final <- history[nrow(history), , drop = FALSE]
M_rank <- rank_truncate(fit$M, cfg$rank)
truth_norm <- max(mat_fnorm(truth$M0), .Machine$double.eps)

summary_df <- data.frame(
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
  nuclear_lambda = fit$lambda,
  lipschitz = fit$lipschitz,
  iter = fit$iter,
  converged = fit$converged,
  hit_maxit = fit$iter >= cfg$nuclear_maxit && !fit$converged,
  final_objective = fit$objective,
  final_rel_objective_change = final$rel_objective_change,
  final_step_frob = final$step_frob,
  final_rel_step_frob = final$rel_step_frob,
  final_truth_frob = mat_fnorm(fit$M - truth$M0),
  final_rel_truth_frob = mat_fnorm(fit$M - truth$M0) / truth_norm,
  final_rank_truncated_truth_frob = mat_fnorm(M_rank - truth$M0),
  rank_1e8 = fit$rank_1e8,
  rank_1e4 = fit$rank_1e4,
  reference_maxit = cfg$reference_maxit,
  reference_tol = cfg$reference_tol,
  reference_iter = if (is.null(ref_fit)) NA_integer_ else ref_fit$iter,
  reference_converged = if (is.null(ref_fit)) NA else ref_fit$converged,
  frob_diff_vs_reference = if (is.null(ref_fit)) NA_real_ else mat_fnorm(fit$M - ref_fit$M),
  rel_frob_diff_vs_reference = if (is.null(ref_fit)) NA_real_ else mat_fnorm(fit$M - ref_fit$M) / truth_norm,
  objective_gap_vs_reference = if (is.null(ref_fit)) NA_real_ else fit$objective - ref_fit$objective
)

history_file <- file.path(cfg$out_root, "nuclear_norm_initial_convergence_history.csv")
summary_file <- file.path(cfg$out_root, "nuclear_norm_initial_convergence_summary.csv")
config_file <- file.path(cfg$out_root, "config.txt")
fig_pdf <- file.path(cfg$out_root, "figures", "fig_nuclear_norm_initial_convergence.pdf")
fig_png <- file.path(cfg$out_root, "figures", "fig_nuclear_norm_initial_convergence.png")

write.csv(history, history_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)
config_lines <- vapply(
  names(cfg),
  function(nm) paste(nm, paste(cfg[[nm]], collapse = ","), sep = " = "),
  character(1)
)
writeLines(config_lines, config_file)
plot_convergence(history, cfg, fit, fig_pdf, fig_png)

message("Done.")
message("Summary: ", summary_file)
message("History: ", history_file)
message("Figure: ", fig_pdf)
print(summary_df)
