# Benchmark: kd-tree ROBIN init (.robinden_data) vs matrix ROBIN init (robinden)
#
# Both functions live in this package, so this is a clean A/B: same binary, same
# compiler flags, same session, only the initialization function differs.
#
#   robinden(as.matrix(dist(X)), K, mp)  # upstream path: O(n^2) distance matrix
#   ktaucenters:::.robinden_data(X, K, mp)             # new path: kd-tree on the data
#
# It reports (1) equivalence of results and (2) the speed gain.
# Run with:  Rscript benchmark/robin_kdtree_benchmark.R

suppressMessages(library(ktaucenters))

make_data <- function(n, p = 2, K = 3, contam = 0.2, seed = 6) {
  set.seed(seed)
  per <- ceiling(n / K)
  centers <- matrix(seq(0, by = 6, length.out = K), nrow = K, ncol = p)
  lab <- rep(seq_len(K), each = per)[seq_len(n)]
  X <- matrix(rnorm(n * p), ncol = p) + centers[lab, , drop = FALSE]
  n_out <- floor(n * contam)
  if (n_out > 0) {
    idx <- sample.int(n, n_out)
    X[idx, ] <- matrix(runif(n_out * p, 2 * min(X), 2 * max(X)), ncol = p)
  }
  list(X = X, labels = lab)
}

med_time <- function(expr, reps = 7) {
  ts <- numeric(reps)
  for (i in seq_len(reps)) ts[i] <- system.time(force(eval.parent(expr)))[["elapsed"]]
  median(ts)
}

## ---------------------------------------------------------------------------
## 1. Equivalence: identical ROBIN seeds across sizes, dimensions, K and seeds
## ---------------------------------------------------------------------------
cat("== Equivalence: .robinden_data vs robinden ==\n")
grid <- expand.grid(n = c(300, 1500, 4500), p = c(2, 5), K = c(3, 5),
                    seed = c(1, 6, 42))
n_match <- 0L
max_idp_diff <- 0
mismatches <- list()
for (r in seq_len(nrow(grid))) {
  g <- grid[r, ]
  X <- make_data(g$n, g$p, g$K, seed = g$seed)$X
  old <- robinden(as.matrix(dist(X)), g$K, 10)
  new <- ktaucenters:::.robinden_data(X, g$K, 10)
  same <- identical(sort(old$centers), sort(new$centers))
  max_idp_diff <- max(max_idp_diff, max(abs(old$idpoints - new$idpoints)))
  if (same) n_match <- n_match + 1L
  else mismatches[[length(mismatches) + 1]] <-
      sprintf("n=%d p=%d K=%d seed=%d: old=%s new=%s", g$n, g$p, g$K, g$seed,
              paste(sort(old$centers), collapse = ","),
              paste(sort(new$centers), collapse = ","))
}
cat(sprintf("  identical center seeds: %d / %d configs\n", n_match, nrow(grid)))
cat(sprintf("  max |idpoints diff|:    %.3e\n", max_idp_diff))
if (length(mismatches)) { cat("  MISMATCHES:\n"); for (m in mismatches) cat("   ", m, "\n") }

## ---------------------------------------------------------------------------
## 2. Speed: init step, matrix path vs kd-tree path
## ---------------------------------------------------------------------------
cat("\n== Init speed (p=2, K=3, mp=10), median of reps ==\n")
cat(sprintf("%-8s %12s %12s %9s %12s\n",
            "n", "matrix(s)", "kd-tree(s)", "speedup", "matrix RAM"))
for (n in c(1000, 3000, 6000, 9000, 12000)) {
  X <- make_data(n)$X
  t_mat <- med_time(quote({ D <- as.matrix(dist(X)); robinden(D, 3, 10) }), reps = 5)
  t_kd  <- med_time(quote(ktaucenters:::.robinden_data(X, 3, 10)), reps = 5)
  ram <- n * n * 8 / 1024^2
  cat(sprintf("%-8d %12.4f %12.4f %8.1fx %9.0f MB\n", n, t_mat, t_kd, t_mat / t_kd, ram))
}

## ---------------------------------------------------------------------------
## 3. kd-tree scaling where the matrix cannot allocate
## ---------------------------------------------------------------------------
cat("\n== kd-tree scaling beyond the matrix ceiling ==\n")
cat(sprintf("%-8s %12s %14s\n", "n", "kd-tree(s)", "matrix RAM"))
for (n in c(30000, 100000, 300000)) {
  X <- make_data(n)$X
  t_kd <- med_time(quote(ktaucenters:::.robinden_data(X, 3, 10)), reps = 3)
  ram <- n * n * 8 / 1024^3
  cat(sprintf("%-8d %12.4f %11.0f GB\n", n, t_kd, ram))
}

cat("\nDone.\n")
