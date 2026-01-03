// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// Constants for numerical stability
const double LOG_2PI = std::log(2.0 * M_PI);
const double SQRT_2 = std::sqrt(2.0);
const double INV_SQRT_2 = 1.0 / SQRT_2;

//------------------------------------------------------------------------------
// Standard normal PDF (log scale)
//------------------------------------------------------------------------------
inline double dnorm_log(double x) {
  return -0.5 * (LOG_2PI + x * x);
}

//------------------------------------------------------------------------------
// Standard normal CDF (log scale, upper tail)
// log(1 - Phi(x)) = log(Phi(-x))
// Uses erfc for numerical stability
//------------------------------------------------------------------------------
inline double pnorm_log_upper(double x) {
  // log(1 - Phi(x)) = log(0.5 * erfc(x / sqrt(2)))
  return std::log(0.5 * std::erfc(x * INV_SQRT_2));
}

//------------------------------------------------------------------------------
// Normal Inverse Mills Ratio: phi(x) / (1 - Phi(x))
// Numerically stable implementation
//------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec NIMR_cpp(const arma::vec& x) {
  int n = x.n_elem;
  arma::vec result(n);

  for (int i = 0; i < n; i++) {
    double xi = x(i);
    // IMR = exp(log(phi(x)) - log(1 - Phi(x)))
    double log_phi = dnorm_log(xi);
    double log_surv = pnorm_log_upper(xi);
    result(i) = std::exp(log_phi - log_surv);
  }

  return result;
}

//------------------------------------------------------------------------------
// Log-likelihood function for General Hazard model
// Returns: scalar log-likelihood value
//------------------------------------------------------------------------------
// [[Rcpp::export]]
double lglik_cpp(double nu, double theta0,
                 const arma::vec& theta, const arma::vec& eta,
                 const arma::mat& X, const arma::mat& Xt,
                 const arma::vec& times, const arma::ivec& status,
                 const arma::ivec& gam) {

  int n = times.n_elem;

  // Precompute common quantities
  double enu = std::exp(nu);
  double e_nu = std::exp(-nu);

  // Log times
  arma::vec yi = arma::log(times);

  // X.theta = Xt %*% theta
  arma::vec X_theta = Xt * theta;

  // Check if AFT model (all gam == 4)
  bool is_aft = arma::all(gam == 4);

  arma::vec v_i(n);
  if (is_aft) {
    v_i.ones();
  } else {
    // X.eta = X %*% eta
    arma::vec X_eta = X * eta;
    // v.i = exp(X.theta * e_nu - X.eta)
    v_i = arma::exp(X_theta * e_nu - X_eta);
  }

  // u.i = enu * yi - X.theta - theta0
  arma::vec u_i = enu * yi - X_theta - theta0;

  // vi.di = v.i - delta
  arma::vec delta = arma::conv_to<arma::vec>::from(status);
  arma::vec vi_di = v_i - delta;

  // Compute log-likelihood components
  double log_lik = 0.0;
  int n_obs = 0;
  double sum_log_vi_obs = 0.0;
  double sum_lgf_obs = 0.0;
  double sum_lgF_vi_di = 0.0;

  for (int i = 0; i < n; i++) {
    double lgf_t = dnorm_log(u_i(i));           // log(f(u_i))
    double lgF_t = pnorm_log_upper(u_i(i));     // log(F(-u_i)) = log(1 - Phi(u_i))

    if (status(i) == 1) {
      n_obs++;
      sum_lgf_obs += lgf_t;
      if (!is_aft) {
        sum_log_vi_obs += std::log(v_i(i));
      }
    }
    sum_lgF_vi_di += lgF_t * vi_di(i);
  }

  // Sum of log(times) for observed
  double sum_yi_obs = 0.0;
  for (int i = 0; i < n; i++) {
    if (status(i) == 1) {
      sum_yi_obs += yi(i);
    }
  }

  // Final log-likelihood
  log_lik = n_obs * nu - sum_yi_obs + sum_log_vi_obs + sum_lgf_obs + sum_lgF_vi_di;

  return log_lik;
}

//------------------------------------------------------------------------------
// Gradient of log-likelihood (lnorm distribution only)
// Returns: vector of gradients [grad_nu, grad_theta0, grad_theta, grad_eta]
//------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec gradLLS_cpp(double nu, double theta0,
                      const arma::vec& theta, const arma::vec& eta,
                      const arma::mat& X, const arma::mat& Xt,
                      const arma::vec& times, const arma::ivec& status,
                      const arma::ivec& gam) {

  int n = times.n_elem;
  int p = X.n_cols;
  int pt = Xt.n_cols;

  // Precompute common quantities
  double enu = std::exp(nu);
  double e_nu = std::exp(-nu);

  arma::vec yi = arma::log(times);
  arma::vec X_theta = Xt * theta;

  bool is_aft = arma::all(gam == 4);

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
  arma::vec delta = arma::conv_to<arma::vec>::from(status);
  arma::vec vi_di = v_i - delta;

  // Compute NIMR
  arma::vec f_F = NIMR_cpp(u_i);

  // log(F(-u_i)) * v_i
  arma::vec log_Fv(n);
  for (int i = 0; i < n; i++) {
    log_Fv(i) = pnorm_log_upper(u_i(i)) * v_i(i);
  }

  // IMR * (v_i - delta)
  arma::vec IMR_vi_di = f_F % vi_di;

  // Initialize gradients
  double grad_nu = 0.0;
  double grad_theta0 = 0.0;
  arma::vec grad_theta(pt, arma::fill::zeros);
  arma::vec grad_eta(p, arma::fill::zeros);

  // Compute observed sums
  int n_obs = arma::sum(status);

  // grad_theta0 = sum(u_i[obs]) + sum(IMR * vi_di)
  for (int i = 0; i < n; i++) {
    if (status(i) == 1) {
      grad_theta0 += u_i(i);  // -fp.fu.i.obs = u_i for observed
    }
  }
  grad_theta0 += arma::sum(IMR_vi_di);

  if (is_aft) {
    // AFT model gradients
    for (int j = 0; j < pt; j++) {
      double sum_obs = 0.0;
      double sum_all = 0.0;
      for (int i = 0; i < n; i++) {
        if (status(i) == 1) {
          sum_obs += u_i(i) * Xt(i, j);
        }
        sum_all += IMR_vi_di(i) * Xt(i, j);
      }
      grad_theta(j) = sum_obs + sum_all;
    }

    // grad_nu for AFT
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; i++) {
      if (status(i) == 1) {
        sum1 += u_i(i) * yi(i);
      }
      sum2 += IMR_vi_di(i) * yi(i);
    }
    grad_nu = n_obs - sum1 * enu - sum2 * enu;

  } else {
    // General model gradients
    // grad_theta
    for (int j = 0; j < pt; j++) {
      double sum_obs = 0.0;
      double sum_IMR = 0.0;
      double sum_logFv = 0.0;
      double sum_status = 0.0;
      for (int i = 0; i < n; i++) {
        if (status(i) == 1) {
          sum_obs += u_i(i) * Xt(i, j);
          sum_status += Xt(i, j);
        }
        sum_IMR += IMR_vi_di(i) * Xt(i, j);
        sum_logFv += log_Fv(i) * Xt(i, j);
      }
      grad_theta(j) = sum_status * e_nu + sum_obs + sum_IMR + sum_logFv * e_nu;
    }

    // grad_eta
    for (int j = 0; j < p; j++) {
      double sum_status = 0.0;
      double sum_logFv = 0.0;
      for (int i = 0; i < n; i++) {
        if (status(i) == 1) {
          sum_status += X(i, j);
        }
        sum_logFv += log_Fv(i) * X(i, j);
      }
      grad_eta(j) = -sum_status - sum_logFv;
    }

    // grad_nu
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
    grad_nu = n_obs - sum_Xtheta_obs * e_nu - sum_ui_yi_obs * enu -
              sum_IMR_yi * enu - sum_logFv_Xtheta * e_nu;
  }

  // Combine into output vector: [nu, theta0, theta, eta]
  arma::vec result(2 + pt + p);
  result(0) = grad_nu;
  result(1) = grad_theta0;
  result.subvec(2, 1 + pt) = grad_theta;
  result.subvec(2 + pt, 1 + pt + p) = grad_eta;

  return result;
}

//------------------------------------------------------------------------------
// Hessian of log-likelihood (lnorm distribution only)
// Returns: (2 + pt + p) x (2 + pt + p) matrix
//------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat hessLLS_cpp(double nu, double theta0,
                      const arma::vec& theta, const arma::vec& eta,
                      const arma::mat& X, const arma::mat& Xt,
                      const arma::vec& times, const arma::ivec& status,
                      const arma::ivec& gam) {

  int n = times.n_elem;
  int p = X.n_cols;
  int pt = Xt.n_cols;
  int D = 2 + pt + p;

  // Precompute common quantities
  double enu = std::exp(nu);
  double e_nu = std::exp(-nu);
  double e2nu = enu * enu;
  double e_2nu = e_nu * e_nu;

  arma::vec yi = arma::log(times);
  arma::vec yi2 = yi % yi;
  arma::vec yi_enu = yi * enu;
  arma::vec yi_enu2 = yi_enu % yi_enu;

  arma::vec X_theta = Xt * theta;
  arma::vec X_theta_e_nu = X_theta * e_nu;

  bool is_aft = arma::all(gam == 4);

  arma::vec v_i(n);
  if (is_aft) {
    v_i.ones();
  } else {
    arma::vec X_eta = X * eta;
    v_i = arma::exp(X_theta * e_nu - X_eta);
  }

  arma::vec u_i = enu * yi - X_theta - theta0;
  arma::vec u_i2 = u_i % u_i;
  arma::vec delta = arma::conv_to<arma::vec>::from(status);
  arma::vec vi_di = v_i - delta;

  // Compute NIMR and related quantities
  arma::vec f_F = NIMR_cpp(u_i);
  arma::vec f_F2 = f_F % f_F;
  arma::vec fp_F = -u_i % f_F;  // f'(u)/F(-u) = -u * f(u)/F(-u)

  // log(F(-u_i)) * v_i
  arma::vec log_Fv(n);
  for (int i = 0; i < n; i++) {
    log_Fv(i) = pnorm_log_upper(u_i(i)) * v_i(i);
  }

  // For observed: fpp/f = -1 + u^2, (f'/f)^2 = u^2
  // fpp.f.obs - fp.f.obs2 = -1 + u^2 - u^2 = -1  (but we compute per element)

  // Initialize Hessian matrix
  arma::mat HESS(D, D, arma::fill::zeros);

  // Separate observed and all data computations
  arma::uvec obs_idx = arma::find(status == 1);
  int n_obs = obs_idx.n_elem;

  arma::vec u_i_obs = u_i.elem(obs_idx);
  arma::vec u_i2_obs = u_i2.elem(obs_idx);
  arma::vec yi_enu_obs = yi_enu.elem(obs_idx);
  arma::vec yi_enu2_obs = yi_enu2.elem(obs_idx);
  arma::mat Xt_obs = Xt.rows(obs_idx);

  // fpp.f.obs - fp.f.obs2 for observed = -1 + u_i2 - u_i2 = -1
  // But actually: fpp/f = -1 + u^2, fp/f = -u, so (fp/f)^2 = u^2
  // fpp.f.obs - fp.f.obs2 = (-1 + u^2) - u^2 = -1
  arma::vec fpp_f_obs = -1.0 + u_i2_obs;
  arma::vec fp_f_obs2 = u_i2_obs;
  arma::vec w_obs = fpp_f_obs - fp_f_obs2;  // = -1 for all observed

  // (fp.F + f.F2) * vi.di
  arma::vec fpF_fF2_vidi = (fp_F + f_F2) % vi_di;

  if (is_aft) {
    // AFT model Hessian

    // hess_theta0
    double hess_theta0 = arma::sum(w_obs) - arma::sum(fpF_fF2_vidi);

    // hess_theta: weighted crossprod
    arma::vec w_theta_obs = w_obs;
    arma::vec w_theta = -fpF_fF2_vidi;

    arma::mat hess_theta = Xt_obs.t() * (Xt_obs.each_col() % w_theta_obs) +
                           Xt.t() * (Xt.each_col() % w_theta);

    // hess_nu
    double hess_nu = arma::sum(w_obs % yi_enu2_obs) - arma::sum(u_i_obs % yi_enu_obs) -
                     arma::sum(fpF_fF2_vidi % yi_enu2) - arma::sum(f_F % vi_di % yi_enu);

    // hess_theta0_theta
    arma::vec hess_theta0_theta = Xt_obs.t() * w_obs - Xt.t() * fpF_fF2_vidi;

    // hess_theta0_nu
    double hess_theta0_nu = -arma::sum(w_obs % yi_enu_obs) + arma::sum(fpF_fF2_vidi % yi_enu);

    // hess_theta_nu
    arma::vec hess_theta_nu = -Xt_obs.t() * (w_obs % yi_enu_obs) +
                               Xt.t() * (fpF_fF2_vidi % yi_enu);

    // Fill Hessian matrix
    HESS(0, 0) = hess_nu;
    HESS(0, 1) = hess_theta0_nu;
    HESS(1, 0) = hess_theta0_nu;
    HESS(1, 1) = hess_theta0;
    HESS.submat(0, 2, 0, 1 + pt) = hess_theta_nu.t();
    HESS.submat(2, 0, 1 + pt, 0) = hess_theta_nu;
    HESS.submat(1, 2, 1, 1 + pt) = hess_theta0_theta.t();
    HESS.submat(2, 1, 1 + pt, 1) = hess_theta0_theta;
    HESS.submat(2, 2, 1 + pt, 1 + pt) = hess_theta;
    // eta block is zero for AFT

  } else {
    // General model Hessian

    // hess_theta0
    double hess_theta0 = arma::sum(w_obs) - arma::sum(fpF_fF2_vidi);

    // Weights for theta
    arma::vec w_theta = -fpF_fF2_vidi + 2.0 * f_F % v_i * e_nu + log_Fv * e_2nu;

    arma::mat hess_theta = Xt_obs.t() * (Xt_obs.each_col() % w_obs) +
                           Xt.t() * (Xt.each_col() % w_theta);

    // hess_eta
    arma::mat hess_eta = X.t() * (X.each_col() % log_Fv);

    // hess_nu
    arma::vec X_theta_obs = X_theta.elem(obs_idx);
    double hess_nu = arma::sum(X_theta_e_nu.elem(obs_idx)) +
                     arma::sum(w_obs % yi_enu2_obs) - arma::sum(u_i_obs % yi_enu_obs) -
                     arma::sum(fpF_fF2_vidi % yi_enu2) +
                     arma::sum(log_Fv % X_theta_e_nu % (1.0 + X_theta_e_nu)) +
                     2.0 * arma::sum(f_F % yi % v_i % X_theta) -
                     arma::sum(f_F % vi_di % yi_enu);

    // hess_theta0_theta
    arma::vec hess_theta0_theta = Xt_obs.t() * w_obs +
                                  Xt.t() * (f_F % v_i) * e_nu -
                                  Xt.t() * fpF_fF2_vidi;

    // hess_theta0_eta
    arma::vec hess_theta0_eta = -X.t() * (f_F % v_i);

    // hess_theta0_nu
    double hess_theta0_nu = -arma::sum(w_obs % yi_enu_obs) +
                            arma::sum(fpF_fF2_vidi % yi_enu) -
                            arma::sum(f_F % v_i % X_theta_e_nu);

    // hess_theta_eta
    arma::vec w_theta_eta = -f_F % v_i - log_Fv * e_nu;
    arma::mat hess_theta_eta = Xt.t() * (X.each_col() % w_theta_eta);

    // hess_theta_nu
    arma::vec Xt_ones_obs = arma::sum(Xt_obs, 0).t();
    arma::vec hess_theta_nu = -Xt_ones_obs * e_nu -
                               Xt_obs.t() * (w_obs % yi_enu_obs) -
                               Xt.t() * (f_F % v_i % X_theta_e_nu) -
                               Xt.t() * (log_Fv % (1.0 + X_theta_e_nu)) * e_nu +
                               Xt.t() * (fpF_fF2_vidi % yi_enu) -
                               Xt.t() * (f_F % yi % v_i);

    // hess_eta_nu
    arma::vec hess_eta_nu = X.t() * (f_F % v_i % yi_enu) +
                            X.t() * (log_Fv % X_theta_e_nu);

    // Fill Hessian matrix
    HESS(0, 0) = hess_nu;
    HESS(0, 1) = hess_theta0_nu;
    HESS(1, 0) = hess_theta0_nu;
    HESS(1, 1) = hess_theta0;

    HESS.submat(0, 2, 0, 1 + pt) = hess_theta_nu.t();
    HESS.submat(2, 0, 1 + pt, 0) = hess_theta_nu;

    HESS.submat(0, 2 + pt, 0, 1 + pt + p) = hess_eta_nu.t();
    HESS.submat(2 + pt, 0, 1 + pt + p, 0) = hess_eta_nu;

    HESS.submat(1, 2, 1, 1 + pt) = hess_theta0_theta.t();
    HESS.submat(2, 1, 1 + pt, 1) = hess_theta0_theta;

    HESS.submat(1, 2 + pt, 1, 1 + pt + p) = hess_theta0_eta.t();
    HESS.submat(2 + pt, 1, 1 + pt + p, 1) = hess_theta0_eta;

    HESS.submat(2, 2, 1 + pt, 1 + pt) = hess_theta;
    HESS.submat(2, 2 + pt, 1 + pt, 1 + pt + p) = hess_theta_eta;
    HESS.submat(2 + pt, 2, 1 + pt + p, 1 + pt) = hess_theta_eta.t();
    HESS.submat(2 + pt, 2 + pt, 1 + pt + p, 1 + pt + p) = hess_eta;
  }

  return HESS;
}
