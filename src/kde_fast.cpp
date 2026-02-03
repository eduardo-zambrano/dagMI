// kde_fast.cpp - Fast kernel density estimation using Rcpp
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// Gaussian kernel
inline double gaussian_kernel(double u) {
  return exp(-0.5 * u * u) / sqrt(2.0 * M_PI);
}

// [[Rcpp::export]]
NumericVector kde_eval_cpp(NumericVector x_eval,
                            NumericVector x_data,
                            double bandwidth) {
  int n_eval = x_eval.size();
  int n_data = x_data.size();
  NumericVector density(n_eval);

  double h = bandwidth;
  double norm_const = 1.0 / (n_data * h);

  for (int i = 0; i < n_eval; i++) {
    double sum = 0.0;
    for (int j = 0; j < n_data; j++) {
      double u = (x_eval[i] - x_data[j]) / h;
      sum += gaussian_kernel(u);
    }
    density[i] = sum * norm_const;
  }

  return density;
}

// [[Rcpp::export]]
NumericMatrix kde_bivariate_eval_cpp(NumericMatrix xy_eval,
                                      NumericMatrix xy_data,
                                      NumericVector bandwidth) {
  int n_eval = xy_eval.nrow();
  int n_data = xy_data.nrow();
  NumericVector density(n_eval);

  double h1 = bandwidth[0];
  double h2 = bandwidth[1];
  double norm_const = 1.0 / (n_data * h1 * h2 * 2.0 * M_PI);

  for (int i = 0; i < n_eval; i++) {
    double sum = 0.0;
    for (int j = 0; j < n_data; j++) {
      double u1 = (xy_eval(i, 0) - xy_data(j, 0)) / h1;
      double u2 = (xy_eval(i, 1) - xy_data(j, 1)) / h2;
      sum += exp(-0.5 * (u1 * u1 + u2 * u2));
    }
    density[i] = sum * norm_const;
  }

  // Return as matrix for consistency
  NumericMatrix result(n_eval, 1);
  result(_, 0) = density;
  return result;
}

// Efficient leave-one-out KDE for cross-validation
// [[Rcpp::export]]
NumericVector kde_loo_cpp(NumericVector x_data, double bandwidth) {
  int n = x_data.size();
  NumericVector density(n);

  double h = bandwidth;
  double norm_const = 1.0 / ((n - 1) * h);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      if (i != j) {
        double u = (x_data[i] - x_data[j]) / h;
        sum += gaussian_kernel(u);
      }
    }
    density[i] = sum * norm_const;
  }

  return density;
}

// Cross-validation score for bandwidth selection
// [[Rcpp::export]]
double kde_cv_score_cpp(NumericVector x_data, double bandwidth) {
  int n = x_data.size();
  double h = bandwidth;

  // Leave-one-out log-likelihood
  double log_lik = 0.0;
  double norm_const = 1.0 / ((n - 1) * h);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      if (i != j) {
        double u = (x_data[i] - x_data[j]) / h;
        sum += gaussian_kernel(u);
      }
    }
    double f_i = sum * norm_const;
    if (f_i > 1e-300) {
      log_lik += log(f_i);
    } else {
      log_lik += log(1e-300);
    }
  }

  return log_lik;
}

// KDE with density truncation (for numerical stability)
// [[Rcpp::export]]
NumericVector kde_eval_truncated_cpp(NumericVector x_eval,
                                      NumericVector x_data,
                                      double bandwidth,
                                      double min_density = 1e-10) {
  NumericVector density = kde_eval_cpp(x_eval, x_data, bandwidth);

  int n = density.size();
  for (int i = 0; i < n; i++) {
    if (density[i] < min_density) {
      density[i] = min_density;
    }
  }

  return density;
}

// Conditional density estimation p(y|x)
// Uses Nadaraya-Watson kernel regression approach
// [[Rcpp::export]]
NumericVector kde_conditional_cpp(NumericVector y_eval,
                                   NumericVector x_condition,
                                   NumericVector y_data,
                                   NumericVector x_data,
                                   double bandwidth_x,
                                   double bandwidth_y,
                                   double min_density = 1e-10) {
  int n_eval = y_eval.size();
  int n_data = y_data.size();
  NumericVector density(n_eval);

  double hx = bandwidth_x;
  double hy = bandwidth_y;

  for (int i = 0; i < n_eval; i++) {
    double num = 0.0;
    double denom = 0.0;

    for (int j = 0; j < n_data; j++) {
      double ux = (x_condition[i] - x_data[j]) / hx;
      double uy = (y_eval[i] - y_data[j]) / hy;

      double kx = gaussian_kernel(ux);
      double ky = gaussian_kernel(uy);

      num += kx * ky;
      denom += kx;
    }

    if (denom > 1e-15) {
      density[i] = num / (denom * hy);
    } else {
      density[i] = min_density;
    }

    if (density[i] < min_density) {
      density[i] = min_density;
    }
  }

  return density;
}

// Product of marginal densities raised to power 1/(n+1)
// for the density correction factor
// [[Rcpp::export]]
NumericVector density_correction_factor_cpp(NumericMatrix x_data,
                                             List kde_objects,
                                             double power) {
  int n_obs = x_data.nrow();
  int n_vars = x_data.ncol();
  NumericVector result(n_obs, 1.0);

  for (int j = 0; j < n_vars; j++) {
    List kde_j = kde_objects[j];
    NumericVector data_j = kde_j["data"];
    double bw_j = kde_j["bandwidth"];

    // Evaluate density at each observation
    NumericVector x_col = x_data(_, j);
    NumericVector dens_j = kde_eval_truncated_cpp(x_col, data_j, bw_j, 1e-10);

    // Multiply into result with power
    for (int i = 0; i < n_obs; i++) {
      result[i] *= pow(dens_j[i], power);
    }
  }

  return result;
}

// Fast computation of multivariate KDE using product kernel
// [[Rcpp::export]]
NumericVector kde_multivariate_cpp(NumericMatrix x_eval,
                                    NumericMatrix x_data,
                                    NumericVector bandwidth) {
  int n_eval = x_eval.nrow();
  int n_data = x_data.nrow();
  int d = x_data.ncol();
  NumericVector density(n_eval);

  double norm_const = 1.0 / n_data;
  for (int k = 0; k < d; k++) {
    norm_const /= bandwidth[k];
  }
  norm_const /= pow(2.0 * M_PI, d / 2.0);

  for (int i = 0; i < n_eval; i++) {
    double sum = 0.0;
    for (int j = 0; j < n_data; j++) {
      double kernel_prod = 1.0;
      for (int k = 0; k < d; k++) {
        double u = (x_eval(i, k) - x_data(j, k)) / bandwidth[k];
        kernel_prod *= exp(-0.5 * u * u);
      }
      sum += kernel_prod;
    }
    density[i] = sum * norm_const;
  }

  return density;
}

// Weighted KDE for bootstrap resampling
// [[Rcpp::export]]
NumericVector kde_weighted_cpp(NumericVector x_eval,
                                NumericVector x_data,
                                NumericVector weights,
                                double bandwidth) {
  int n_eval = x_eval.size();
  int n_data = x_data.size();
  NumericVector density(n_eval);

  double h = bandwidth;
  double sum_weights = sum(weights);
  double norm_const = 1.0 / (sum_weights * h);

  for (int i = 0; i < n_eval; i++) {
    double sum = 0.0;
    for (int j = 0; j < n_data; j++) {
      double u = (x_eval[i] - x_data[j]) / h;
      sum += weights[j] * gaussian_kernel(u);
    }
    density[i] = sum * norm_const;
  }

  return density;
}
