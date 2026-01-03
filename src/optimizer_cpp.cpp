// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace Numer;

// Constants for numerical stability
static const double LOG_2PI = std::log(2.0 * M_PI);
static const double SQRT_2 = std::sqrt(2.0);
static const double INV_SQRT_2 = 1.0 / SQRT_2;

//------------------------------------------------------------------------------
// Standard normal PDF (log scale)
//------------------------------------------------------------------------------
inline double dnorm_log(double x) {
  return -0.5 * (LOG_2PI + x * x);
}

//------------------------------------------------------------------------------
// Standard normal CDF (log scale, upper tail)
//------------------------------------------------------------------------------
inline double pnorm_log_upper(double x) {
  return std::log(0.5 * std::erfc(x * INV_SQRT_2));
}

//------------------------------------------------------------------------------
// Normal Inverse Mills Ratio
//------------------------------------------------------------------------------
inline double nimr_scalar(double x) {
  return std::exp(dnorm_log(x) - pnorm_log_upper(x));
}

//------------------------------------------------------------------------------
// Class for optimization problem - computes negative log-likelihood and gradient
//------------------------------------------------------------------------------
class NegLogLikelihood: public MFuncGrad {
private:
  arma::mat X;
  arma::mat Xt;
  arma::vec times;
  arma::ivec status;
  arma::uvec theta_keep;
  arma::uvec eta_keep;
  int n, p, pt, ptheta, peta;
  bool is_aft;

  // Pre-computed values
  arma::vec yi;
  arma::vec delta;
  int n_obs;

public:
  NegLogLikelihood(const arma::mat& X_, const arma::mat& Xt_,
                   const arma::vec& times_, const arma::ivec& status_,
                   const arma::uvec& theta_keep_, const arma::uvec& eta_keep_,
                   bool is_aft_) :
    X(X_), Xt(Xt_), times(times_), status(status_),
    theta_keep(theta_keep_), eta_keep(eta_keep_), is_aft(is_aft_) {

    n = times.n_elem;
    p = X.n_cols;
    pt = Xt.n_cols;
    ptheta = arma::sum(theta_keep);
    peta = arma::sum(eta_keep);

    // Pre-compute
    yi = arma::log(times);
    delta = arma::conv_to<arma::vec>::from(status);
    n_obs = arma::sum(status);
  }

  // Compute negative log-likelihood and gradient
  double f_grad(Constvec& par, Refvec grad) {
    // Parse parameters from optimization vector
    double nu = par[0];
    double theta0 = par[1];

    // Reconstruct full theta and eta vectors
    arma::vec theta(pt, arma::fill::zeros);
    arma::vec eta(p, arma::fill::zeros);

    int idx = 2;
    for (int j = 0; j < pt; j++) {
      if (theta_keep(j) == 1) {
        theta(j) = par[idx++];
      }
    }
    for (int j = 0; j < p; j++) {
      if (eta_keep(j) == 1) {
        eta(j) = par[idx++];
      }
    }

    // Common quantities
    double enu = std::exp(nu);
    double e_nu = std::exp(-nu);

    arma::vec X_theta = Xt * theta;
    arma::vec v_i(n);
    arma::vec X_eta(n);

    if (is_aft) {
      v_i.ones();
      X_eta.zeros();
    } else {
      X_eta = X * eta;
      v_i = arma::exp(X_theta * e_nu - X_eta);
    }

    arma::vec u_i = enu * yi - X_theta - theta0;
    arma::vec vi_di = v_i - delta;

    // Compute NIMR and related quantities
    arma::vec f_F(n);
    arma::vec log_Fv(n);
    for (int i = 0; i < n; i++) {
      f_F(i) = nimr_scalar(u_i(i));
      log_Fv(i) = pnorm_log_upper(u_i(i)) * v_i(i);
    }
    arma::vec IMR_vi_di = f_F % vi_di;

    // === Compute log-likelihood ===
    double sum_yi_obs = 0.0;
    double sum_log_vi_obs = 0.0;
    double sum_lgf_obs = 0.0;
    double sum_lgF_vi_di = 0.0;

    for (int i = 0; i < n; i++) {
      double lgf_t = dnorm_log(u_i(i));
      double lgF_t = pnorm_log_upper(u_i(i));

      if (status(i) == 1) {
        sum_yi_obs += yi(i);
        sum_lgf_obs += lgf_t;
        if (!is_aft) {
          sum_log_vi_obs += std::log(v_i(i));
        }
      }
      sum_lgF_vi_di += lgF_t * vi_di(i);
    }

    double ll = n_obs * nu - sum_yi_obs + sum_log_vi_obs + sum_lgf_obs + sum_lgF_vi_di;

    // === Compute gradients ===
    double grad_theta0 = 0.0;
    for (int i = 0; i < n; i++) {
      if (status(i) == 1) {
        grad_theta0 += u_i(i);
      }
    }
    grad_theta0 += arma::sum(IMR_vi_di);

    // Initialize gradient output
    grad.setZero();

    if (is_aft) {
      // AFT model gradients
      arma::vec grad_theta_full(pt, arma::fill::zeros);
      for (int j = 0; j < pt; j++) {
        double sum_obs = 0.0;
        double sum_all = 0.0;
        for (int i = 0; i < n; i++) {
          if (status(i) == 1) {
            sum_obs += u_i(i) * Xt(i, j);
          }
          sum_all += IMR_vi_di(i) * Xt(i, j);
        }
        grad_theta_full(j) = sum_obs + sum_all;
      }

      double sum1 = 0.0, sum2 = 0.0;
      for (int i = 0; i < n; i++) {
        if (status(i) == 1) {
          sum1 += u_i(i) * yi(i);
        }
        sum2 += IMR_vi_di(i) * yi(i);
      }
      double grad_nu = n_obs - sum1 * enu - sum2 * enu;

      grad[0] = -grad_nu;
      grad[1] = -grad_theta0;
      idx = 2;
      for (int j = 0; j < pt; j++) {
        if (theta_keep(j) == 1) {
          grad[idx++] = -grad_theta_full(j);
        }
      }

    } else {
      // General model gradients
      arma::vec grad_theta_full(pt, arma::fill::zeros);
      arma::vec grad_eta_full(p, arma::fill::zeros);

      for (int j = 0; j < pt; j++) {
        double sum_status = 0.0;
        double sum_obs = 0.0;
        double sum_IMR = 0.0;
        double sum_logFv = 0.0;
        for (int i = 0; i < n; i++) {
          if (status(i) == 1) {
            sum_obs += u_i(i) * Xt(i, j);
            sum_status += Xt(i, j);
          }
          sum_IMR += IMR_vi_di(i) * Xt(i, j);
          sum_logFv += log_Fv(i) * Xt(i, j);
        }
        grad_theta_full(j) = sum_status * e_nu + sum_obs + sum_IMR + sum_logFv * e_nu;
      }

      for (int j = 0; j < p; j++) {
        double sum_status = 0.0;
        double sum_logFv = 0.0;
        for (int i = 0; i < n; i++) {
          if (status(i) == 1) {
            sum_status += X(i, j);
          }
          sum_logFv += log_Fv(i) * X(i, j);
        }
        grad_eta_full(j) = -sum_status - sum_logFv;
      }

      double sum_Xtheta_obs = 0.0;
      double sum_ui_yi_obs = 0.0;
      double sum_IMR_yi = 0.0;
      double sum_logFv_Xtheta = 0.0;
      for (int i = 0; i < n; i++) {
        if (status(i) == 1) {
          sum_Xtheta_obs += X_theta(i);
          sum_ui_yi_obs += u_i(i) * yi(i);
        }
        sum_IMR_yi += IMR_vi_di(i) * yi(i);
        sum_logFv_Xtheta += log_Fv(i) * X_theta(i);
      }
      double grad_nu = n_obs - sum_Xtheta_obs * e_nu - sum_ui_yi_obs * enu -
                       sum_IMR_yi * enu - sum_logFv_Xtheta * e_nu;

      grad[0] = -grad_nu;
      grad[1] = -grad_theta0;
      idx = 2;
      for (int j = 0; j < pt; j++) {
        if (theta_keep(j) == 1) {
          grad[idx++] = -grad_theta_full(j);
        }
      }
      for (int j = 0; j < p; j++) {
        if (eta_keep(j) == 1) {
          grad[idx++] = -grad_eta_full(j);
        }
      }
    }

    return -ll;
  }
};

//------------------------------------------------------------------------------
// Main optimization function - replaces mle.gh in R
//------------------------------------------------------------------------------
// [[Rcpp::export]]
List mle_gh_cpp(const arma::mat& X, const arma::mat& Xt,
                const arma::vec& times, const arma::ivec& status,
                const arma::ivec& gam,
                int max_iter = 100,
                double eps_f = 1e-6,
                double eps_g = 1e-5) {

  int n = times.n_elem;
  int p = X.n_cols;
  int pt = Xt.n_cols;
  int gam_len = gam.n_elem;

  // Determine which parameters to include
  arma::uvec theta_keep(pt, arma::fill::zeros);
  arma::uvec eta_keep(p, arma::fill::zeros);

  bool is_aft = true;
  for (int j = 0; j < gam_len; j++) {
    if (gam(j) != 4 && gam(j) != 0) {
      is_aft = false;
      break;
    }
  }

  for (int j = 0; j < std::min(pt, gam_len); j++) {
    theta_keep(j) = (gam(j) == 1 || gam(j) == 3 || gam(j) == 4) ? 1 : 0;
  }
  for (int j = 0; j < std::min(p, gam_len); j++) {
    eta_keep(j) = (gam(j) == 2 || gam(j) == 3) ? 1 : 0;
  }

  int ptheta = arma::sum(theta_keep);
  int peta = arma::sum(eta_keep);
  int D = 2 + ptheta + peta;

  // Initialize parameter vector
  Eigen::VectorXd par(D);
  par.setZero();

  // Create optimization problem
  NegLogLikelihood nll(X, Xt, times, status, theta_keep, eta_keep, is_aft);

  // Run L-BFGS optimization
  double fopt;
  int niter = optim_lbfgs(nll, par, fopt, max_iter, eps_f, eps_g);

  // Reconstruct full parameter vectors
  arma::vec theta_out(pt, arma::fill::zeros);
  arma::vec eta_out(p, arma::fill::zeros);

  double nu_out = par[0];
  double theta0_out = par[1];

  int idx = 2;
  for (int j = 0; j < pt; j++) {
    if (theta_keep(j) == 1) {
      theta_out(j) = par[idx++];
    }
  }
  for (int j = 0; j < p; j++) {
    if (eta_keep(j) == 1) {
      eta_out(j) = par[idx++];
    }
  }

  return List::create(
    Named("nu") = nu_out,
    Named("theta0") = theta0_out,
    Named("theta") = theta_out,
    Named("eta") = eta_out,
    Named("convergence") = (niter < max_iter) ? 0 : 1,
    Named("iterations") = niter,
    Named("value") = -fopt
  );
}
