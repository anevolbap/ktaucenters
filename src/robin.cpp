#include "knn.h"
#include "utils.h"
#include "nanoflann.hpp"
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

NumericVector point_density(NumericMatrix D, const std::size_t k) {

  //' Estimates the local points density.
  //'
  //' @param D a distance matrix, which contains the distances between the rows
  //' of a matrix.
  //' @param k number of neighbors to calculate local point
  //' density.
  //'
  //' @return
  //' A vector containing the density values for each point.

  const std::size_t n = D.nrow();

  List knn = dist_to_kNN(D, k);
  IntegerMatrix id = knn["id"];
  NumericMatrix distances = knn["dist"];

  NumericVector out(no_init(n));

  for (std::size_t i = 0; i < n; ++i) {
    NumericVector max_distance(k);
    for (std::size_t j = 0; j < k; ++j) {
      max_distance[j] = std::max(distances(id(i, j), k - 1), distances(i, j));
    }
    out[i] = k / sum(max_distance);
  }
  return out;
}

std::size_t robin_center(NumericVector idp, IntegerVector indices,
                         const double crit_robin) {
  //' Utility function to estimate robinden center
  //'
  //' @param idp a vector with containing the inverse density each point.
  //' @param indices vector with sorted indices.
  //' @param crit_robin critical robin value.
  //'
  //' @return
  //' Index of the cluster center

  NumericVector idp_sorted_points = idp[indices];

  // Sometimes all idp_sorted_points are greater than the
  // crit_robin value, then we take the nearest point to crit_robin
  LogicalVector comp = idp_sorted_points <= crit_robin;

  if (is_true(any(comp))) {
    IntegerVector tmp = indices[comp];
    return tmp[0];
  } else {
    return indices[which_min(idp_sorted_points - crit_robin)];
  }
}

//' Robust Initialization based on Inverse Density estimator (ROBINDEN)
//'
//' Searches for k initial cluster seeds for k-means based clustering methods.
//'
//' @param D a distance matrix, which contains the distances between the rows of
//' a matrix.
//' @param n_clusters number of cluster centers to find.
//' @param mp number of nearest neighbors to compute point density.
//'
//' @return A list with the following components:
//' \item{\code{centers}}{: A numeric vector with the initial cluster centers
//' indices.}
//' \item{\code{idpoints}}{: A real vector containing the inverse of point
//' density estimation.}
//'
//' @details
//' The centers are the observations located in the most dense region
//' and far away from each other at the same time.
//' In order to find the observations in the highly dense region, this function
//' uses point density estimation (instead of Local Outlier Factor, Breunig et
//' al (2000)), see more details.
//'
//' @note This is a slightly modified version of ROBIN algorithm
//' implementation done by Sarka Brodinova <sarka.brodinova@tuwien.ac.at>.
//' @author Juan Domingo Gonzalez <juanrst@hotmail.com>
//'
//' @examples
//' # Generate synthetic data (7 cluster well separated)
//' K <- 5
//' nk <- 100
//' Z <- rnorm(2 * K * nk)
//' mues <- rep(5 * -floor(K/2):floor(K/2), 2 * nk * K)
//' X <-  matrix(Z + mues, ncol = 2)
//'
//' # Generate synthetic outliers (contamination level 20%)
//' X[sample(1:(nk * K), (nk * K) * 0.2), ] <-
//'   matrix(runif((nk * K) * 0.2 * 2, 3 * min(X), 3 * max(X)),
//'          ncol = 2,
//'          nrow = (nk * K)* 0.2)
//' res <- robinden(D = as.matrix(dist(X)), n_clusters = K, mp = 10);
//' # plot the Initial centers found
//' plot(X)
//' points(X[res$centers, ], pch = 19, col = 4, cex = 2)
//'
//' @references Hasan AM, et al. Robust partitional clustering by
//' outlier and density insensitive seeding. Pattern Recognition Letters,
//' 30(11), 994-1002, 2009.
//'
//'@export
// [[Rcpp::export]]
List robinden(NumericMatrix D, const std::size_t n_clusters,
              const std::size_t mp) {

  const std::size_t n = D.nrow();

  // Compute the inverse density points.
  NumericVector idp = 1.0 / point_density(D, mp);

  // Outliers have a high idp value. In unbalanced cases and when the number of
  // clusters increases, all the observations from a group might be above the
  // crit_robin. So we need to increase the crit_robin in order to
  // avoid two initials centers from the same group.

  // Minus 1 to get 0 based index
  const std::size_t position =
      trunc(std::max(0.5, 0.96 * (1 - (1.5 / n_clusters))) * n) - 1;
  NumericVector sorted_idp = clone(idp).sort(false);
  const double crit_robin = sorted_idp[position];

  // Start with a point with maximum density
  std::size_t r = which_min(idp);
  IntegerVector sorted_indices = top_index(D.column(r), n, true);

  IntegerVector centers(n_clusters);
  centers[0] = robin_center(idp, sorted_indices, crit_robin);

  for (std ::size_t iter = 1; iter < n_clusters; ++iter) {
    NumericVector minimum_values(no_init(n));

    for (std::size_t column = 0; column < n; ++column) {
      NumericVector c = D.column(column);
      NumericVector tmp = c[centers[Range(0, iter - 1)]];
      minimum_values[column] = min(tmp);
    }

    sorted_indices = top_index(minimum_values, n, true);
    centers[iter] = robin_center(idp, sorted_indices, crit_robin);
  }
  return List::create(_["centers"] = centers, _["idpoints"] = idp);
}

// ---- kd-tree fast path (internal) ----
// Same algorithm and tie behavior as robinden(), but sources the kNN from a
// kd-tree on the data and computes center distances on the fly, avoiding the
// O(n^2) distance matrix. Used by ktaucenters()/ktaucentersfast().

struct MatAdaptor {
  const NumericMatrix &mat;
  MatAdaptor(const NumericMatrix &m) : mat(m) {}
  inline std::size_t kdtree_get_point_count() const { return mat.nrow(); }
  inline double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const {
    return mat(idx, dim);
  }
  template <class BBOX> bool kdtree_get_bbox(BBOX &) const { return false; }
};
typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, MatAdaptor>, MatAdaptor>
    KDTree;

static inline double row_dist(const NumericMatrix &x, std::size_t a,
                              std::size_t b) {
  const std::size_t p = x.ncol();
  double s = 0.0, d;
  for (std::size_t j = 0; j < p; ++j) {
    d = x(a, j) - x(b, j);
    s += d * d;
  }
  return std::sqrt(s);
}

// [[Rcpp::export(".robinden_data")]]
List robinden_data(NumericMatrix x, const std::size_t n_clusters,
                   const std::size_t mp) {
  const std::size_t n = x.nrow();
  const std::size_t p = x.ncol();

  // Build the kd-tree and query the mp nearest neighbors of each point.
  MatAdaptor adaptor(x);
  KDTree index(p, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  NumericMatrix knn_dist(n, mp);
  IntegerMatrix knn_id(n, mp);
  std::vector<double> kdist(n);
  std::vector<unsigned int> ret_idx(mp + 1);
  std::vector<double> ret_d2(mp + 1);
  std::vector<double> query(p);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < p; ++j)
      query[j] = x(i, j);
    index.knnSearch(&query[0], mp + 1, &ret_idx[0], &ret_d2[0]);
    std::size_t out = 0;
    for (std::size_t r = 0; r < mp + 1 && out < mp; ++r) {
      if ((std::size_t)ret_idx[r] == i)
        continue; // drop self
      knn_id(i, out) = (int)ret_idx[r];
      knn_dist(i, out) = std::sqrt(ret_d2[r]);
      ++out;
    }
    kdist[i] = knn_dist(i, mp - 1);
  }

  // Inverse point density (same definition as point_density()).
  NumericVector idp(n);
  for (std::size_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (std::size_t j = 0; j < mp; ++j)
      s += std::max(kdist[knn_id(i, j)], knn_dist(i, j));
    idp[i] = s / (double)mp;
  }

  const std::size_t position =
      (std::size_t)(std::trunc(
          std::max(0.5, 0.96 * (1 - (1.5 / n_clusters))) * n)) -
      1;
  NumericVector sorted_idp = clone(idp).sort(false);
  const double crit_robin = sorted_idp[position];

  std::size_t r = which_min(idp);
  NumericVector dist_r(no_init(n));
  for (std::size_t c = 0; c < n; ++c)
    dist_r[c] = row_dist(x, c, r);
  IntegerVector sorted_indices = top_index(dist_r, n, true);

  IntegerVector centers(n_clusters);
  centers[0] = robin_center(idp, sorted_indices, crit_robin);

  // Running minimum distance to the centers chosen so far.
  NumericVector mindist(no_init(n));
  for (std::size_t c = 0; c < n; ++c)
    mindist[c] = row_dist(x, c, centers[0]);

  for (std::size_t iter = 1; iter < n_clusters; ++iter) {
    sorted_indices = top_index(mindist, n, true);
    centers[iter] = robin_center(idp, sorted_indices, crit_robin);
    std::size_t cc = centers[iter];
    for (std::size_t c = 0; c < n; ++c) {
      double dd = row_dist(x, c, cc);
      if (dd < mindist[c])
        mindist[c] = dd;
    }
  }

  return List::create(_["centers"] = centers, _["idpoints"] = idp);
}
