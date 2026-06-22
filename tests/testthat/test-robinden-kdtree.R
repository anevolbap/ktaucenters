test_that(".robinden_data matches robinden on the distance matrix", {
  make_data <- function(n, p, K, seed) {
    set.seed(seed)
    per <- ceiling(n / K)
    centers <- matrix(seq(0, by = 6, length.out = K), nrow = K, ncol = p)
    lab <- rep(seq_len(K), each = per)[seq_len(n)]
    matrix(rnorm(n * p), ncol = p) + centers[lab, , drop = FALSE]
  }

  grid <- expand.grid(n = c(300, 1500), p = c(2, 5), K = c(3, 5),
                      seed = c(1, 6, 42))
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    X <- make_data(g$n, g$p, g$K, g$seed)
    old <- robinden(as.matrix(dist(X)), g$K, 10)
    new <- ktaucenters:::.robinden_data(X, g$K, 10)
    expect_equal(sort(new$centers), sort(old$centers))
    expect_equal(new$idpoints, old$idpoints, tolerance = 1e-10)
  }
})
