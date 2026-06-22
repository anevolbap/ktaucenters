# kd-tree ROBIN init: benchmark results

Hardware/software: single machine, R 4.2.2, g++ -O2. Reproduce with
`Rscript benchmark/robin_kdtree_benchmark.R` (times vary by machine; ratios are
the point). Data: K well-separated Gaussian clusters with 20% uniform
contamination.

## 1. Same results (equivalence)

`.robinden_data` (kd-tree) vs `robinden(as.matrix(dist(X)))` (matrix), swept over
n in {300, 1500, 4500}, p in {2, 5}, K in {3, 5}, seed in {1, 6, 42}:

- identical center seeds: **36 / 36** configurations
- max |idpoints difference|: **7.1e-15** (floating-point noise)

Full `ktaucentersfast` output vs CRAN 1.0.0 (matrix path), n in {1500, 4500, 9000}:

| n    | centers diff | clusters | outliers | di diff | tau diff | ARI (kd) | ARI (CRAN) |
|------|--------------|----------|----------|---------|----------|----------|------------|
| 1500 | 0            | equal    | equal    | 0       | 0        | 1.0000   | 1.0000     |
| 4500 | 0            | equal    | equal    | 0       | 0        | 1.0000   | 1.0000     |
| 9000 | 0            | equal    | equal    | 0       | 0        | 1.0000   | 1.0000     |

Output is bit-identical and clustering accuracy (Adjusted Rand Index vs ground
truth) is unchanged.

## 2. Speed: init step (matrix vs kd-tree)

| n     | matrix (s) | kd-tree (s) | speedup | matrix RAM |
|-------|------------|-------------|---------|------------|
| 1000  | 0.019      | 0.001       | 19x     | 8 MB       |
| 3000  | 0.424      | 0.004       | 106x    | 69 MB      |
| 6000  | 1.517      | 0.009       | 169x    | 275 MB     |
| 9000  | 3.459      | 0.013       | 266x    | 618 MB     |
| 12000 | 6.213      | 0.017       | 366x    | 1099 MB    |

The speedup grows with n because the matrix path is O(n^2) and the kd-tree path
is ~O(n log n).

## 3. Scaling beyond the matrix ceiling

kd-tree init at sizes where the n x n matrix cannot allocate:

| n      | kd-tree (s) | matrix RAM needed |
|--------|-------------|-------------------|
| 30000  | 0.048       | 7 GB              |
| 100000 | 0.195       | 75 GB             |
| 300000 | 0.675       | 671 GB            |

## 4. Full pipeline (ktaucentersfast vs CRAN 1.0.0)

| n    | speedup |
|------|---------|
| 1500 | 1.8x    |
| 4500 | 6.4x    |
| 9000 | 6.7x    |

## 5. Tests

`tests/testthat`: `rho-opt` plus a new `robinden-kdtree` test asserting the
kd-tree path equals the matrix path across 24 configurations. All pass.
