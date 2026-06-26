// MCMC inner loop in C++ using RcppArmadillo.
//
// Invoked from R when both densityFunction and proposalFunction passed to
// mcmcSampling() carry an "rcpp_spec" attribute attached by
// mclustDensityFunction() and addHighDimGaussian(). The R side hands the raw
// arrays to mcmc_loop_cpp() and gets a matrix of samples back. No callbacks
// into R during the loop.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Minimum proposal-scale floor; keep in sync with R/mcmcSampling.R.
static const double MIN_COV_CORRECTION = 1e-10;

// Precomputed parameters of a single Gaussian mixture.
struct GMMSpec {
  arma::cube inv_sigma;  // d x d x G
  arma::vec  log_norm;   // G  (log(pro_k) - d/2 log(2pi) - 0.5 log|sigma_k|)
  arma::mat  means;      // d x G
};

// Full target-density configuration: env GMM gates the species GMM and the
// returned value is floored at floor_value.
struct DensityConfig {
  GMMSpec env;
  GMMSpec sp;
  double env_threshold;
  double sp_cutoff;
  double floor_value;
};

static inline double gmm_density(const arma::vec& x, const GMMSpec& g) {
  const arma::uword G = g.log_norm.n_elem;
  arma::vec log_d(G);
  for (arma::uword k = 0; k < G; ++k) {
    arma::vec diff = x - g.means.col(k);
    log_d[k] = g.log_norm[k] - 0.5 * arma::as_scalar(diff.t() * g.inv_sigma.slice(k) * diff);
  }
  const double m = log_d.max();
  return std::exp(m + std::log(arma::sum(arma::exp(log_d - m))));
}

static inline double combined_density(const arma::vec& x, const DensityConfig& cfg) {
  const double env_d = gmm_density(x, cfg.env);
  if (env_d < cfg.env_threshold) return cfg.floor_value;
  // Uniform-sampling contract: mclustDensityFunction(species.model = NULL) sets
  // cfg.sp_cutoff = Inf as the sentinel for "no presence exclusion". Detect it and
  // return the in-support target 1.0 directly; the species GMM passed in that mode
  // is dummy marshalling ballast and is never evaluated here.
  if (std::isinf(cfg.sp_cutoff)) return 1.0;
  const double sp_d = gmm_density(x, cfg.sp);
  const double val = 1.0 - sp_d / cfg.sp_cutoff;
  return val > cfg.floor_value ? val : cfg.floor_value;
}

// Parameter groups (no inline comments inside the signature so
// Rcpp::compileAttributes() can parse it):
//   env GMM:     env_inv_sigma, env_log_norm, env_means, env_threshold
//   species GMM: sp_inv_sigma, sp_log_norm, sp_means, sp_cutoff, floor_value
//   proposal:    proposal_mean, proposal_cov
// [[Rcpp::export]]
List mcmc_loop_cpp(NumericVector start_point,
                   int n_sample_points,
                   int burn_in,
                   double covariance_correction,
                   arma::cube env_inv_sigma,
                   arma::vec env_log_norm,
                   arma::mat env_means,
                   double env_threshold,
                   arma::cube sp_inv_sigma,
                   arma::vec sp_log_norm,
                   arma::mat sp_means,
                   double sp_cutoff,
                   double floor_value,
                   arma::vec proposal_mean,
                   arma::mat proposal_cov) {

  RNGScope rngScope;  // honors set.seed() from the caller

  DensityConfig cfg;
  cfg.env = {env_inv_sigma, env_log_norm, env_means};
  cfg.sp  = {sp_inv_sigma,  sp_log_norm,  sp_means};
  cfg.env_threshold = env_threshold;
  cfg.sp_cutoff     = sp_cutoff;
  cfg.floor_value   = floor_value;

  const int d = start_point.size();
  arma::vec current(d);
  for (int j = 0; j < d; ++j) current[j] = start_point[j];

  // Lower-triangular Cholesky factor for the proposal covariance. Combined
  // with sqrt(cov_correction), this gives proposals from N(0, c * Sigma).
  // The R caller already validates positive-definiteness, so this won't throw.
  arma::mat L = arma::chol(proposal_cov, "lower");

  double current_density = combined_density(current, cfg);

  // Pre-allocate per-step buffers once so the inner loops have zero heap
  // traffic on every iteration.
  arma::vec z(d, arma::fill::zeros);
  arma::vec proposed(d, arma::fill::zeros);

  // ---- Burn-in: Robbins-Monro on log(covariance_correction) ----
  const double target_accept = 0.234;
  const double rm_exponent = 0.6;
  double log_cc = std::log(covariance_correction);
  int burn_accepted = 0;
  for (int t = 1; t <= burn_in; ++t) {
    const double scale = std::sqrt(std::exp(log_cc));
    for (int j = 0; j < d; ++j) z[j] = R::rnorm(0.0, 1.0);
    proposed = current + proposal_mean + scale * (L * z);
    const double proposed_density = combined_density(proposed, cfg);
    bool accept = false;
    if (!ISNAN(proposed_density) && proposed_density > 0.0 && current_density > 0.0) {
      const double ratio = proposed_density / current_density;
      if (ratio >= 1.0) {
        accept = true;
      } else if (R::runif(0.0, 1.0) < ratio) {
        accept = true;
      }
    }
    if (accept) {
      current = proposed;
      current_density = proposed_density;
      burn_accepted += 1;
    }
    const double gamma_t = 1.0 / std::pow(static_cast<double>(t) + 1.0, rm_exponent);
    log_cc += gamma_t * ((accept ? 1.0 : 0.0) - target_accept);
  }
  double cc_final = std::exp(log_cc);
  if (cc_final < MIN_COV_CORRECTION) cc_final = MIN_COV_CORRECTION;

  // ---- Sampling at the adapted, frozen step size ----
  arma::mat samples(n_sample_points, d + 1, arma::fill::zeros);
  const double scale = std::sqrt(cc_final);
  int rejected = 0;
  for (int i = 0; i < n_sample_points; ++i) {
    for (int j = 0; j < d; ++j) z[j] = R::rnorm(0.0, 1.0);
    proposed = current + proposal_mean + scale * (L * z);
    const double proposed_density = combined_density(proposed, cfg);
    bool accept = false;
    if (!ISNAN(proposed_density) && proposed_density > 0.0 && current_density > 0.0) {
      const double ratio = proposed_density / current_density;
      if (ratio >= 1.0) {
        accept = true;
      } else if (R::runif(0.0, 1.0) < ratio) {
        accept = true;
      }
    }
    if (accept) {
      current = proposed;
      current_density = proposed_density;
    } else {
      rejected += 1;
    }
    for (int j = 0; j < d; ++j) samples(i, j) = current[j];
    samples(i, d) = current_density;
  }

  return List::create(
    _["samples"] = samples,
    _["covariance_correction"] = cc_final,
    _["burnin_accepted"] = burn_accepted,
    _["sampling_rejected"] = rejected
  );
}
