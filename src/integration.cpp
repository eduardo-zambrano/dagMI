// integration.cpp - Numerical integration routines for Q_n computation
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// Gaussian kernel for KDE evaluation
inline double gaussian_kernel(double u) {
  return exp(-0.5 * u * u) / sqrt(2.0 * M_PI);
}

// Evaluate marginal density at a single point (for integration)
double eval_marginal_density(double x,
                              NumericVector data,
                              double bandwidth,
                              double min_dens = 1e-10) {
  int n = data.size();
  double h = bandwidth;
  double sum = 0.0;

  for (int j = 0; j < n; j++) {
    double u = (x - data[j]) / h;
    sum += gaussian_kernel(u);
  }

  double dens = sum / (n * h);
  return std::max(dens, min_dens);
}

// Evaluate bivariate density f(x_{j-1}, x_j) for consecutive pair
double eval_bivariate_density(double x1, double x2,
                               NumericVector data1,
                               NumericVector data2,
                               double bw1, double bw2,
                               double min_dens = 1e-10) {
  int n = data1.size();
  double sum = 0.0;

  for (int i = 0; i < n; i++) {
    double u1 = (x1 - data1[i]) / bw1;
    double u2 = (x2 - data2[i]) / bw2;
    sum += exp(-0.5 * (u1 * u1 + u2 * u2));
  }

  double dens = sum / (n * bw1 * bw2 * 2.0 * M_PI);
  return std::max(dens, min_dens);
}

// [[Rcpp::export]]
double compute_qn_quadrature_cpp(List kde_list,
                                  NumericMatrix quad_grid,
                                  NumericVector quad_weights,
                                  int n_vars) {
  // Q_n^G = integral of product of bivariate densities raised to 1/(n+1)
  // over marginal densities raised to (n-1)/(n+1)

  int n_points = quad_grid.nrow();
  double integral = 0.0;

  double power_bivar = 1.0 / (n_vars + 1);
  double power_marg = (n_vars - 1.0) / (n_vars + 1);

  for (int i = 0; i < n_points; i++) {
    double integrand = 1.0;

    // Product over bivariate densities (consecutive pairs)
    for (int j = 1; j < n_vars; j++) {
      List kde_j_prev = kde_list[j - 1];
      List kde_j = kde_list[j];

      NumericVector data_prev = kde_j_prev["data"];
      NumericVector data_j = kde_j["data"];
      double bw_prev = kde_j_prev["bandwidth"];
      double bw_j = kde_j["bandwidth"];

      double f_bivar = eval_bivariate_density(
        quad_grid(i, j - 1), quad_grid(i, j),
        data_prev, data_j,
        bw_prev, bw_j
      );

      integrand *= pow(f_bivar, power_bivar);
    }

    // Product over marginal densities
    for (int j = 0; j < n_vars; j++) {
      List kde_j = kde_list[j];
      NumericVector data_j = kde_j["data"];
      double bw_j = kde_j["bandwidth"];

      double f_marg = eval_marginal_density(
        quad_grid(i, j), data_j, bw_j
      );

      integrand /= pow(f_marg, power_marg);
    }

    integral += quad_weights[i] * integrand;
  }

  return integral;
}

// Monte Carlo integration using importance sampling
// [[Rcpp::export]]
List compute_qn_montecarlo_cpp(List kde_list,
                                NumericMatrix samples,
                                int n_vars,
                                int n_samples) {
  double power_bivar = 1.0 / (n_vars + 1);
  double power_marg = (n_vars - 1.0) / (n_vars + 1);

  NumericVector values(n_samples);

  for (int i = 0; i < n_samples; i++) {
    double integrand = 1.0;

    // Product over bivariate densities
    for (int j = 1; j < n_vars; j++) {
      List kde_j_prev = kde_list[j - 1];
      List kde_j = kde_list[j];

      NumericVector data_prev = kde_j_prev["data"];
      NumericVector data_j = kde_j["data"];
      double bw_prev = kde_j_prev["bandwidth"];
      double bw_j = kde_j["bandwidth"];

      double f_bivar = eval_bivariate_density(
        samples(i, j - 1), samples(i, j),
        data_prev, data_j,
        bw_prev, bw_j
      );

      integrand *= pow(f_bivar, power_bivar);
    }

    // Product over marginal densities
    for (int j = 0; j < n_vars; j++) {
      List kde_j = kde_list[j];
      NumericVector data_j = kde_j["data"];
      double bw_j = kde_j["bandwidth"];

      double f_marg = eval_marginal_density(
        samples(i, j), data_j, bw_j
      );

      integrand /= pow(f_marg, power_marg);
    }

    values[i] = integrand;
  }

  // Compute mean and standard error
  double qn = mean(values);
  double se = sd(values) / sqrt((double)n_samples);

  return List::create(
    Named("qn") = qn,
    Named("se") = se,
    Named("values") = values
  );
}

// Compute T_h statistic for a given test function
// [[Rcpp::export]]
double compute_test_stat_cpp(NumericMatrix data,
                              List kde_list,
                              NumericVector h_values,
                              double qn,
                              int n_vars) {
  int n_obs = data.nrow();

  double power_dens = 1.0 / (n_vars + 1);
  double power_h = 1.0 / (n_vars + 1);

  // Compute LHS: E[prod h_i(X_i) * prod p_i(X_i)^{1/(n+1)}]
  double lhs = 0.0;

  for (int i = 0; i < n_obs; i++) {
    double term = 1.0;

    for (int j = 0; j < n_vars; j++) {
      // h_j(X_j) value
      int idx = i * n_vars + j;
      term *= h_values[idx];

      // p_j(X_j)^{1/(n+1)}
      List kde_j = kde_list[j];
      NumericVector data_j = kde_j["data"];
      double bw_j = kde_j["bandwidth"];

      double dens_j = eval_marginal_density(data(i, j), data_j, bw_j);
      term *= pow(dens_j, power_dens);
    }

    lhs += term;
  }
  lhs /= n_obs;

  // Compute RHS: Q_n * prod (E[h_i(X_i)^{n+1}])^{1/(n+1)}
  double rhs = qn;

  for (int j = 0; j < n_vars; j++) {
    double moment = 0.0;
    for (int i = 0; i < n_obs; i++) {
      int idx = i * n_vars + j;
      moment += pow(h_values[idx], n_vars + 1);
    }
    moment /= n_obs;
    rhs *= pow(moment, power_h);
  }

  // T_h = RHS - LHS (should be >= 0 under H0)
  return rhs - lhs;
}

// Compute density correction factor for all observations
// [[Rcpp::export]]
NumericVector compute_density_correction_cpp(NumericMatrix data,
                                              List kde_list,
                                              int n_vars) {
  int n_obs = data.nrow();
  double power = 1.0 / (n_vars + 1);

  NumericVector correction(n_obs);

  for (int i = 0; i < n_obs; i++) {
    double prod = 1.0;

    for (int j = 0; j < n_vars; j++) {
      List kde_j = kde_list[j];
      NumericVector data_j = kde_j["data"];
      double bw_j = kde_j["bandwidth"];

      double dens_j = eval_marginal_density(data(i, j), data_j, bw_j);
      prod *= pow(dens_j, power);
    }

    correction[i] = prod;
  }

  return correction;
}

// Compute (n+1)-th moments for test statistic denominator
// [[Rcpp::export]]
NumericVector compute_moments_cpp(NumericMatrix h_values,
                                   int n_vars,
                                   int n_obs) {
  NumericVector moments(n_vars);
  double power = 1.0 / (n_vars + 1);

  for (int j = 0; j < n_vars; j++) {
    double sum = 0.0;
    for (int i = 0; i < n_obs; i++) {
      sum += pow(h_values(i, j), n_vars + 1);
    }
    moments[j] = pow(sum / n_obs, power);
  }

  return moments;
}

// Gauss-Legendre quadrature nodes and weights for [-1, 1]
// [[Rcpp::export]]
List gauss_legendre_cpp(int n) {
  NumericVector nodes(n);
  NumericVector weights(n);

  int m = (n + 1) / 2;

  for (int i = 0; i < m; i++) {
    // Initial guess
    double z = cos(M_PI * (i + 0.75) / (n + 0.5));
    double p1 = 1.0;
    double p2 = 0.0;

    // Newton iteration
    for (int iter = 0; iter < 100; iter++) {
      p1 = 1.0;
      p2 = 0.0;

      for (int j = 1; j <= n; j++) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }

      double pp = n * (z * p1 - p2) / (z * z - 1.0);
      double z_new = z - p1 / pp;

      if (fabs(z_new - z) < 1e-15) break;
      z = z_new;
    }

    // Final Legendre polynomial evaluation at converged z
    p1 = 1.0;
    p2 = 0.0;
    for (int j = 1; j <= n; j++) {
      double p3 = p2;
      p2 = p1;
      p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
    }

    nodes[i] = -z;
    nodes[n - 1 - i] = z;
    double pp = n * (z * p1 - p2) / (z * z - 1.0);
    double w = 2.0 / ((1.0 - z * z) * pp * pp);
    weights[i] = w;
    weights[n - 1 - i] = w;
  }

  return List::create(
    Named("nodes") = nodes,
    Named("weights") = weights
  );
}

// Transform quadrature from [-1,1] to [a,b]
// [[Rcpp::export]]
List transform_quadrature_cpp(NumericVector nodes,
                               NumericVector weights,
                               double a, double b) {
  double mid = (a + b) / 2.0;
  double half_width = (b - a) / 2.0;

  int n = nodes.size();
  NumericVector new_nodes(n);
  NumericVector new_weights(n);

  for (int i = 0; i < n; i++) {
    new_nodes[i] = mid + half_width * nodes[i];
    new_weights[i] = half_width * weights[i];
  }

  return List::create(
    Named("nodes") = new_nodes,
    Named("weights") = new_weights
  );
}
