/*
  Copyright (C) 2019 Andrew D. Smith
  Author: Andrew D. Smith

  This is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This software is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
*/

#include "TwoStateHMM.hpp"

#include <memory>
#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#include "smithlab_utils.hpp"

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::min;
using std::abs;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;
using std::isfinite;
using std::make_pair;

using smithlab::log_sum_log_vec;

struct betabin {
  betabin(const double a, const double b) :
    alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}
  double operator()(const pair<double, double> &val) const;
  void fit(const vector<double> &vals_a, const vector<double> &vals_b,
           const vector<double> &p);
  string tostring() const;
  double alpha;
  double beta;
  double lnbeta_helper;
};

string
betabin::tostring() const {
  std::ostringstream os;
  os << std::fixed << setprecision(3) << alpha << " "
     << std::fixed << setprecision(3) << beta;
  return os.str();
}

double
betabin::operator()(const pair<double, double> &val) const {
  const size_t x = static_cast<size_t>(val.first);
  const size_t n = static_cast<size_t>(x + val.second);
  return gsl_sf_lnchoose(n, x) +
    gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
}

inline static double
sign(double x) {
  return (x >= 0) ? 1.0 : -1.0;
}
static const double tolerance = 1e-10;

inline static double
invpsi(const double tolerance, const double x) {
  double L = 1.0, Y = std::exp(x);
  while (L > tolerance) {
    Y += L*sign(x - gsl_sf_psi(Y));
    L /= 2.0;
  }
  return Y;
}

inline static double
movement(const double curr, const double prev) {
  return std::abs(curr - prev)/max(std::abs(curr), std::abs(prev));
}

void
betabin::fit(const vector<double> &vals_a, const vector<double> &vals_b,
             const vector<double> &p) {
  const double p_total = std::accumulate(begin(p), end(p), 0.0);

  const double alpha_rhs =
    inner_product(begin(vals_a), end(vals_a), begin(p), 0.0)/p_total;
  const double beta_rhs =
    inner_product(begin(vals_b), end(vals_b), begin(p), 0.0)/p_total;

  double prev_alpha = 0.0, prev_beta = 0.0;
  alpha = beta = 0.01;
  while (movement(alpha, prev_alpha) > tolerance &&
         movement(beta, prev_beta) > tolerance) {
    prev_alpha = alpha;
    prev_beta = beta;
    alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
    beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
  }
  lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

inline double
log_sum_log(const double p, const double q) {
  if (p == 0.0) {return q;}
  else if (q == 0.0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

static double
log_sum_log_vec(const vector<double> &vals, size_t a, size_t b) {
  auto x = std::max_element(begin(vals) + a, begin(vals) + b);
  const double max_val = *x;
  const size_t max_idx = x - begin(vals);
  double sum = 1.0;
  for (size_t i = a; i < b; ++i) {
    if (i != max_idx)
      sum += exp(vals[i] - max_val);
  }
  return max_val + log(sum);
}

static double
log_sum_log_vec(const vector<double> &vals, const vector<size_t> &resets) {
  vector<double> w(resets.size() - 1);
  for (size_t i = 0; i < resets.size() - 1; ++i)
    w[i] = log_sum_log_vec(vals, resets[i], resets[i+1] - 1);
  return log_sum_log_vec(w, w.size());
}


////////////////////////////////////////////////////////////////////////
/////////// WRAPPER FUNCTIONS
////////////////////////////////////////////////////////////////////////

double
TwoStateHMM::ViterbiDecoding(const vector<pair<double, double> > &values,
                             const vector<size_t> &reset_points,
                             const double p_fb, const double p_bf,
                             const double fg_alpha, const double fg_beta,
                             const double bg_alpha, const double bg_beta,
                             vector<bool> &classes) const {

  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  return ViterbiDecoding(values, reset_points, p_fb, p_bf,
                         fg_distro, bg_distro, classes);
}


double
TwoStateHMM::PosteriorDecoding(const vector<pair<double, double> > &values,
                               const vector<size_t> &reset_points,
                               const double p_fb, const double p_bf,
                               const double fg_alpha, const double fg_beta,
                               const double bg_alpha, const double bg_beta,
                               vector<bool> &classes,
                               vector<double> &posteriors) const {

  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  return PosteriorDecoding(values, reset_points, p_fb, p_bf,
                           fg_distro, bg_distro, classes, posteriors);
}

void
TwoStateHMM::PosteriorScores(const vector<pair<double, double> > &values,
                             const vector<size_t> &reset_points,
                             const double p_fb, const double p_bf,
                             const double fg_alpha, const double fg_beta,
                             const double bg_alpha, const double bg_beta,
                             const bool fg_class,
                             vector<double> &posteriors) const {

  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  PosteriorScores(values, reset_points, p_fb, p_bf,
                  fg_distro, bg_distro, fg_class, posteriors);
}

double
TwoStateHMM::BaumWelchTraining(const std::vector<pair<double, double> > &values,
                               const std::vector<size_t> &reset_points,
                               double &p_fb, double &p_bf,
                               double &fg_alpha, double &fg_beta,
                               double &bg_alpha, double &bg_beta) const {

  betabin fg_distro(fg_alpha, fg_beta);
  betabin bg_distro(bg_alpha, bg_beta);

  const double score = BaumWelchTraining(values, reset_points, p_fb, p_bf,
                                         fg_distro, bg_distro);
  fg_alpha = fg_distro.alpha;
  fg_beta = fg_distro.beta;
  bg_alpha = bg_distro.alpha;
  bg_beta = bg_distro.beta;

  return score;
}

////////////////////////////////////////////////////////////////////////
/////////// INTERNAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

static void
get_emissions(vector<pair<double, double> >::const_iterator v,
              const vector<pair<double, double> >::const_iterator v_end,
              vector<double>::iterator emit, const betabin &distr) {
  while (v != v_end)
    *emit++ = distr(*v++);
}


double
forward_algorithm(const size_t start, const size_t end,
                  const double lp_sf, const double lp_sb,
                  const double lp_ff, const double lp_fb,
                  const double lp_bf, const double lp_bb,
                  const vector<double> &fg_emit, const vector<double> &bg_emit,
                  vector<pair<double, double> > &f) {
  f[start].first = fg_emit[start] + lp_sf;
  f[start].second = bg_emit[start] + lp_sb;
  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    f[i].first = fg_emit[i] + log_sum_log(f[k].first + lp_ff,
                                          f[k].second + lp_bf);
    f[i].second = bg_emit[i] + log_sum_log(f[k].first + lp_fb,
                                           f[k].second + lp_bb);
  }
  return log_sum_log(f[end - 1].first, f[end - 1].second);
}

double
backward_algorithm(const size_t start, const size_t end,
                   const double lp_sf, const double lp_sb,
                   const double lp_ff, const double lp_fb,
                   const double lp_bf, const double lp_bb,
                   const vector<double> &fg_emit,
                   const vector<double> &bg_emit,
                   vector<pair<double, double> > &b) {
  b[end - 1].first = 0;
  b[end - 1].second = 0;
  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    const double fg_a = fg_emit[k] + b[k].first;
    const double bg_a = bg_emit[k] + b[k].second;
    b[i].first = log_sum_log(fg_a + lp_ff, bg_a + lp_fb);
    b[i].second = log_sum_log(fg_a + lp_bf, bg_a + lp_bb);
  }
  return log_sum_log(b[start].first + fg_emit[start] + lp_sf,
                     b[start].second + bg_emit[start] + lp_sb);
}


template <class T> void
one_minus(T a, const T a_end, T b) {
  while (a != a_end)
    *b++ = 1.0 - *a++;
}

inline static double
get_posterior(const pair<double, double> &f, const pair<double, double> &b) {
  const double fg = f.first + b.first;
  return exp(fg - log_sum_log(fg, f.second + b.second));
}

inline static void
get_posteriors(const vector<pair<double, double> > &forward,
               const vector<pair<double, double> > &backward,
               vector<double> &posteriors) {
  posteriors.resize(forward.size());
  for (size_t i = 0; i < forward.size(); ++i)
    posteriors[i] = get_posterior(forward[i], backward[i]);
}


void
summarize_transitions(const size_t start, const size_t end,
                      const vector<pair<double, double> > &f,
                      const vector<pair<double, double> > &b,
                      const double total,
                      const vector<double> &fg_emit, const vector<double> &bg_emit,
                      const double lp_ff, const double lp_fb,
                      const double lp_bf, const double lp_bb,
                      vector<double> &ff_vals, vector<double> &fb_vals,
                      vector<double> &bf_vals, vector<double> &bb_vals) {

  for (size_t i = start + 1; i < end; ++i) {
    const double b_first = b[i].first + fg_emit[i] - total;
    const double b_second = b[i].second + bg_emit[i] - total;

    const size_t k = i - 1;
    const double ff = f[k].first;
    const double bb = f[k].second;

    ff_vals[k] = ff + lp_ff + b_first;
    fb_vals[k] = ff + lp_fb + b_second;

    bf_vals[k] = bb + lp_bf + b_first;
    bb_vals[k] = bb + lp_bb + b_second;
  }
}


double
single_iteration(const vector<pair<double, double> > &values,
                 const vector<double> &vals_a, const vector<double> &vals_b,
                 const vector<size_t> &reset_points,
                 vector<pair<double, double> > &forward,
                 vector<pair<double, double> > &backward,
                 double &p_fb, double &p_bf,
                 betabin &fg_distro, betabin &bg_distro,
                 vector<double> &fg_emit, vector<double> &bg_emit,
                 vector<double> &ff_vals, vector<double> &fb_vals,
                 vector<double> &bf_vals, vector<double> &bb_vals) {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);
  assert(isfinite(lp_sf) && isfinite(lp_sb) && isfinite(lp_ff) &&
         isfinite(lp_fb) && isfinite(lp_bf) && isfinite(lp_bb));

  get_emissions(begin(values), end(values), begin(fg_emit), fg_distro);
  get_emissions(begin(values), end(values), begin(bg_emit), bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    summarize_transitions(reset_points[i], reset_points[i + 1],
                          forward, backward, score, fg_emit, bg_emit,
                          lp_ff, lp_fb, lp_bf, lp_bb,
                          ff_vals, fb_vals, bf_vals, bb_vals);

    total_loglik += score;
  }

  const double p_ff_update = exp(log_sum_log_vec(ff_vals, reset_points));
  const double p_fb_update = exp(log_sum_log_vec(fb_vals, reset_points));
  const double f_denom = p_ff_update + p_fb_update;
  assert(p_fb_update/f_denom > tolerance);
  p_fb = p_fb_update/f_denom;

  const double p_bf_update = exp(log_sum_log_vec(bf_vals, reset_points));
  const double p_bb_update = exp(log_sum_log_vec(bb_vals, reset_points));
  const double b_denom = p_bb_update + p_bf_update;
  assert(p_bf_update/b_denom > tolerance);
  p_bf = p_bf_update/b_denom;

  vector<double> posteriors;
  get_posteriors(forward, backward, posteriors);
  fg_distro.fit(vals_a, vals_b, posteriors);

  one_minus(begin(posteriors), end(posteriors), begin(posteriors));
  bg_distro.fit(vals_a, vals_b, posteriors);

  return total_loglik;
}


static void
report_param_header_for_verbose() {
  cerr << setw(3) << "ITR"
       << setw(8) << "F size"
       << setw(8) << "B size"
       << setw(14) << "F PARAMS"
       << setw(14) << "B PARAMS"
       << setw(11) << "DELTA"
       << endl;
}

inline double
get_delta(const double a, const double b) {
  return (b - a)/max(abs(a), abs(b));
}


static void
report_params_for_verbose(const size_t i,
                          const double p_fb_est,
                          const double p_bf_est,
                          const betabin &fg_distro,
                          const betabin &bg_distro,
                          const double total,
                          const double prev_total) {
  std::ios_base::fmtflags orig_flags(cerr.flags());
  cerr.precision(2);
  cerr << setw(3) << i + 1
       << setw(8) << std::fixed << 1/p_fb_est
       << setw(8) << std::fixed << 1/p_bf_est
       << setw(14) << fg_distro.tostring()
       << setw(14) << bg_distro.tostring()
       << setw(11) << std::scientific
       << abs(get_delta(prev_total, total))
       << endl;
  cerr.flags(orig_flags);
}


static void
extract_fractional_values(const vector<pair<double, double> > &values,
                          vector<double> &vals_a, vector<double> &vals_b) {

  const size_t n_vals = values.size();
  vals_a.resize(n_vals);
  vals_b.resize(n_vals);

  static const double epsilon = 1e-2;
  // const double pseudocount = 1.0;
  for (size_t i = 0; i < n_vals; ++i) {
    // const double a = values[i].first + pseudocount;
    // const double b = values[i].second + pseudocount;
    // const double val = a/(a + b);
    const double val = values[i].first/(values[i].first + values[i].second);
    const double adjusted_val = min(max(val, epsilon), 1 - epsilon);
    vals_a[i] = log(adjusted_val);
    vals_b[i] = log(1.0 - adjusted_val);
  }
}


double
TwoStateHMM::BaumWelchTraining(const vector<pair<double, double> > &values,
                               const vector<size_t> &reset_points,
                               double &p_fb, double &p_bf,
                               betabin &fg_distro, betabin &bg_distro) const {

  vector<double> vals_a, vals_b;
  extract_fractional_values(values, vals_a, vals_b);

  const size_t n_vals = values.size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> ff_vals(n_vals), fb_vals(n_vals); // for estimating transitions
  vector<double> bf_vals(n_vals), bb_vals(n_vals);
  vector<double> fg_emit(n_vals), bg_emit(n_vals); // avoid recomp of emissions

  if (VERBOSE)
    report_param_header_for_verbose();

  double prev_total = -std::numeric_limits<double>::max();
  bool converged = false;
  for (size_t i = 0; i < max_iterations && !converged; ++i) {

    double p_fb_est = p_fb, p_bf_est = p_bf;

    const double total =
      single_iteration(values, vals_a, vals_b, reset_points, forward, backward,
                       p_fb_est, p_bf_est, fg_distro, bg_distro,
                       ff_vals, fb_vals, bf_vals, bb_vals, fg_emit, bg_emit);

    if (VERBOSE)
      report_params_for_verbose(i, p_fb_est, p_bf_est,
                                fg_distro, bg_distro, total, prev_total);

    // ADS: removing the check based on expected log likelihood from
    // forward/backward as these seem to have some problem...
    converged = ((get_delta(p_fb_est, p_fb) < tolerance) &&
		 (get_delta(p_bf_est, p_bf) < tolerance));
    // converged = (get_delta(prev_total, total) < tolerance);

    if (converged) {
      if (VERBOSE)
        cerr << "CONVERGED" << endl;
    }
    else {
      p_fb = p_fb_est;
      p_bf = p_bf_est;
      prev_total = total;
    }
  }
  return prev_total;
}


void
TwoStateHMM::PosteriorScores(const vector<pair<double, double> > &values,
                             const vector<size_t> &reset_points,
                             const double p_fb, const double p_bf,
                             const betabin &fg_distro, const betabin &bg_distro,
                             const bool fg_class,
                             vector<double> &posteriors) const {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  assert(isfinite(lp_sf) && isfinite(lp_sb) && isfinite(lp_ff) &&
         isfinite(lp_fb) && isfinite(lp_bf) && isfinite(lp_bb));

  const size_t n_vals = values.size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> fg_emit(n_vals), bg_emit(n_vals);
  get_emissions(begin(values), end(values), begin(fg_emit), fg_distro);
  get_emissions(begin(values), end(values), begin(bg_emit), bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    total_loglik += score;
  }

  get_posteriors(forward, backward, posteriors);
  if (!fg_class)
    one_minus(begin(posteriors), end(posteriors), begin(posteriors));
}


void
TwoStateHMM::TransitionPosteriors(const vector<pair<double, double> > &values,
                                  const vector<size_t> &reset_points,
                                  const double p_fb, const double p_bf,
                                  const double fg_alpha, const double fg_beta,
                                  const double bg_alpha, const double bg_beta,
                                  const size_t transition,
                                  vector<double> &posteriors) const {

  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  return TransitionPosteriors(values, reset_points, p_fb, p_bf,
                              fg_distro, bg_distro, transition, posteriors);
}


static void
get_joint_posteriors(const pair<double, double> &forward,
                     const pair<double, double> &backward,
                     const double fg_emit, const double bg_emit,
                     const double lp_ff, const double lp_fb,
                     const double lp_bf, const double lp_bb,
                     double &ff_c, double &fb_c, double &bf_c, double &bb_c) {
  // (forward val) + transition + emission + (backward val offset by 1)
  ff_c = forward.first + lp_ff + fg_emit + backward.first;
  fb_c = forward.first + lp_fb + bg_emit + backward.second;
  bf_c = forward.second + lp_bf + fg_emit + backward.first;
  bb_c = forward.second + lp_bb + bg_emit + backward.second;
}


void
TwoStateHMM::TransitionPosteriors(const vector<pair<double, double> > &values,
                                  const vector<size_t> &reset_points,
                                  const double p_fb, const double p_bf,
                                  const betabin &fg_distro,
                                  const betabin &bg_distro,
                                  const size_t transition,
                                  vector<double> &scores) const {

  assert(transition < 4);

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  const size_t n_vals = values.size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> fg_emit(n_vals), bg_emit(n_vals);
  get_emissions(begin(values), end(values), begin(fg_emit), fg_distro);
  get_emissions(begin(values), end(values), begin(bg_emit), bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    total_loglik += score;
  }

  scores.clear();
  scores.resize(n_vals, 0.0);
  size_t j = 0;
  for (size_t i = 0; i < n_vals; ++i) {
    if (i == reset_points[j])
      ++j;
    else {
      double ff_c, fb_c, bf_c, bb_c;
      get_joint_posteriors(forward[i - 1], backward[i], fg_emit[i], bg_emit[i],
                           lp_ff, lp_fb, lp_bf, lp_bb, ff_c, fb_c, bf_c, bb_c);
      double numerator = ff_c;
      if (transition == 1) numerator = fb_c;
      if (transition == 2) numerator = bf_c;
      if (transition == 3) numerator = bb_c;
      const double denom = log_sum_log(log_sum_log(ff_c, fb_c),
                                       log_sum_log(bf_c, bb_c));
      scores[i] = exp(numerator - denom);
    }
  }
}


double
TwoStateHMM::PosteriorDecoding(const vector<pair<double, double> > &values,
                               const vector<size_t> &reset_points,
                               const double p_fb, const double p_bf,
                               const betabin &fg_distro,
                               const betabin &bg_distro,
                               vector<bool> &classes,
                               vector<double> &posteriors) const {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  const size_t n_vals = values.size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> fg_emit(n_vals), bg_emit(n_vals);
  get_emissions(begin(values), end(values), begin(fg_emit), fg_distro);
  get_emissions(begin(values), end(values), begin(bg_emit), bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    total_loglik += score;
  }

  get_posteriors(forward, backward, posteriors);

  classes.resize(n_vals);
  for (size_t i = 0; i < n_vals; ++i)
    classes[i] = (posteriors[i] > 0.5);

  return total_loglik;
}


double
TwoStateHMM::ViterbiDecoding(const vector<pair<double, double> > &values,
                             const vector<size_t> &reset_points,
                             const double p_fb, const double p_bf,
                             const betabin &fg_distro, const betabin &bg_distro,
                             vector<bool> &ml_classes) const {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  double total = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const size_t lim = reset_points[i + 1] - reset_points[i];

    vector<pair<double, double> > v(lim, make_pair(0.0, 0.0));
    vector<pair<size_t, size_t> > trace(lim, make_pair(0ul, 0ul));

    v[0].first = fg_distro(values[reset_points[i]]) + lp_sf;
    v[0].second = bg_distro(values[reset_points[i]]) + lp_sb;

    for (size_t j = 1; j < lim; ++j) {

      const double ff = v[j - 1].first + lp_ff;
      const double bf = v[j - 1].second + lp_bf;
      const double fg_log_emmit = fg_distro(values[reset_points[i] + j]);
      if (ff > bf) {
        v[j].first = fg_log_emmit + ff;
        trace[j].first = 0;
      }
      else {
        v[j].first = fg_log_emmit + bf;
        trace[j].first = 1;
      }

      const double fb = v[j - 1].first + lp_fb;
      const double bb = v[j - 1].second + lp_bb;
      const double bg_log_emmit = bg_distro(values[reset_points[i] + j]);
      if (fb > bb) {
        v[j].second = bg_log_emmit + fb;
        trace[j].second = 0;
      }
      else {
        v[j].second = bg_log_emmit + bb;
        trace[j].second = 1;
      }
    }

    vector<bool> inner_ml_classes;

    // do the traceback
    size_t prev = 0;
    if (v.back().first > v.back().second) {
      inner_ml_classes.push_back(true);
      prev = trace.back().first;
    }
    else {
      inner_ml_classes.push_back(false);
      prev = trace.back().second;
    }

    for (size_t j = trace.size() - 1; j > 0; --j) {
      const size_t k = j - 1;
      if (prev == 0) {
        inner_ml_classes.push_back(true);
        prev = trace[k].first;
      }
      else {
        inner_ml_classes.push_back(false);
        prev = trace[k].second;
      }
    }
    reverse(begin(inner_ml_classes), end(inner_ml_classes));
    ml_classes.insert(end(ml_classes),
                      begin(inner_ml_classes), end(inner_ml_classes));

    total += max(v.back().first, v.back().second);
  }
  return total;
}


////////////////////////////////////////////////////////////////////////////////
///////////////  FOR MULTIPLE REPLICATES

// WRAPPER FUNCTIONS

double
TwoStateHMM::BaumWelchTraining(const vector<vector<pair<double, double> > > &values,
                               const vector<size_t> &reset_points,
                               double &p_fb, double &p_bf,
                               vector<double> &fg_alpha,
                               vector<double> &fg_beta,
                               vector<double> &bg_alpha,
                               vector<double> &bg_beta) const {
  vector<betabin> fg_distro, bg_distro;
  for (size_t i = 0; i < values.size(); ++i) {
    fg_distro.push_back(betabin(fg_alpha[i], fg_beta[i]));
    bg_distro.push_back(betabin(bg_alpha[i], bg_beta[i]));
  }

  const double score = BaumWelchTraining(values, reset_points,
                                         p_fb, p_bf, fg_distro, bg_distro);
  for (size_t r = 0; r < values.size(); ++r) {
    fg_alpha[r] = fg_distro[r].alpha;
    fg_beta[r] = fg_distro[r].beta;
    bg_alpha[r] = bg_distro[r].alpha;
    bg_beta[r] = bg_distro[r].beta;
  }
  return score;
}

double
TwoStateHMM::PosteriorDecoding(const vector<vector<pair<double, double> > > &values,
                               const vector<size_t> &reset_points,
                               const double p_fb, const double p_bf,
                               const vector<double> &fg_alpha,
                               const vector<double> &fg_beta,
                               const vector<double> &bg_alpha,
                               const vector<double> &bg_beta,
                               vector<bool> &classes,
                               vector<double> &posteriors) const {

  vector<betabin> fg_distro, bg_distro;
  for (size_t i = 0; i < values.size(); ++i) {
    fg_distro.push_back(betabin(fg_alpha[i], fg_beta[i]));
    bg_distro.push_back(betabin(bg_alpha[i], bg_beta[i]));
  }
  return PosteriorDecoding(values, reset_points, p_fb, p_bf,
                           fg_distro, bg_distro, classes, posteriors);
}


void
TwoStateHMM::PosteriorScores(const vector<vector<pair<double, double> > > &values,
                             const vector<size_t> &reset_points,
                             const double p_fb, const double p_bf,
                             const vector<double> &fg_alpha,
                             const vector<double> &fg_beta,
                             const vector<double> &bg_alpha,
                             const vector<double> &bg_beta,
                             const bool &fg_class,
                             vector<double> &posteriors) const {

  vector<betabin> fg_distro, bg_distro;
  for (size_t i = 0; i < values.size(); ++i) {
    fg_distro.push_back(betabin(fg_alpha[i], fg_beta[i]));
    bg_distro.push_back(betabin(bg_alpha[i], bg_beta[i]));
  }
  return PosteriorScores(values, reset_points, p_fb, p_bf,
                         fg_distro, bg_distro, fg_class, posteriors);
}

// INTERNAL FUNCTIONS (FOR REPLICATES)

inline bool
has_data(const pair<double, double> &p) {
  return p.first + p.second >= 1.0;
}

static void
get_emissions_rep(const vector<vector<pair<double, double> > > &v,
                  vector<double> &emit, const vector<betabin> &distr) {
  fill(begin(emit), end(emit), 0.0);
  for (size_t r = 0; r < v.size(); ++r)
    for (size_t i = 0; i < v[r].size(); ++i)
      if (has_data(v[r][i]))
        emit[i] += distr[r](v[r][i]);
}


static void
fit_distro_rep(betabin &distro, const vector<pair<double, double> > &values,
               const vector<double> &vals_a, const vector<double> &vals_b,
               const vector<double> &posteriors,
               vector<double> &tmp_a, vector<double> &tmp_b,
               vector<double> &tmp_p) {
  tmp_a.clear();
  tmp_b.clear();
  tmp_p.clear();
  for (size_t i = 0; i < values.size(); ++i)
    if (has_data(values[i])) {
      tmp_a.push_back(vals_a[i]);
      tmp_b.push_back(vals_b[i]);
      tmp_p.push_back(posteriors[i]);
    }
  distro.fit(tmp_a, tmp_b, tmp_p);
}


double
single_iteration_rep(const vector<vector<pair<double, double> > > &values,
                     const vector<vector<double> > &vals_a,
                     const vector<vector<double> > &vals_b,
                     const vector<size_t> &reset_points,
                     vector<pair<double, double> > &forward,
                     vector<pair<double, double> > &backward,
                     double &p_fb, double &p_bf,
                     vector<betabin> &fg_distro, vector<betabin> &bg_distro,
                     vector<double> &fg_emit, vector<double> &bg_emit,
                     vector<double> &ff_vals, vector<double> &fb_vals,
                     vector<double> &bf_vals, vector<double> &bb_vals) {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  get_emissions_rep(values, fg_emit, fg_distro);
  get_emissions_rep(values, bg_emit, bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    summarize_transitions(reset_points[i], reset_points[i + 1],
                          forward, backward, score, fg_emit, bg_emit,
                          lp_ff, lp_fb, lp_bf, lp_bb,
                          ff_vals, fb_vals, bf_vals, bb_vals);

    total_loglik += score;
  }

  const double p_ff_update = exp(log_sum_log_vec(ff_vals, reset_points));
  const double p_fb_update = exp(log_sum_log_vec(fb_vals, reset_points));
  const double f_denom = p_ff_update + p_fb_update;
  assert(p_fb_update/f_denom > tolerance);
  p_fb = p_fb_update/f_denom;

  const double p_bf_update = exp(log_sum_log_vec(bf_vals, reset_points));
  const double p_bb_update = exp(log_sum_log_vec(bb_vals, reset_points));
  const double b_denom = p_bb_update + p_bf_update;
  assert(p_bf_update/b_denom > tolerance);
  p_bf = p_bf_update/b_denom;

  vector<double> posteriors;
  get_posteriors(forward, backward, posteriors);

  const size_t n_reps = values.size();
  const size_t n_vals = values[0].size();

  vector<double> tmp_a, tmp_b, tmp_p;
  tmp_a.reserve(n_vals);
  tmp_b.reserve(n_vals);
  tmp_p.reserve(n_vals);
  for (size_t r = 0; r < n_reps; ++r)
    fit_distro_rep(fg_distro[r], values[r],
                   vals_a[r], vals_b[r], posteriors, tmp_a, tmp_b, tmp_p);

  one_minus(begin(posteriors), end(posteriors), begin(posteriors));
  for (size_t r = 0; r < n_reps; ++r)
    fit_distro_rep(bg_distro[r], values[r],
                   vals_a[r], vals_b[r], posteriors, tmp_a, tmp_b, tmp_p);

  return total_loglik;
}



double
TwoStateHMM::BaumWelchTraining(const vector<vector<pair<double, double> > > &values,
                               const vector<size_t> &reset_points,
                               double &p_fb, double &p_bf,
                               vector<betabin> &fg_distro,
                               vector<betabin> &bg_distro) const {

  // extract the fractional values (both fraction meth and unmeth)
  const size_t n_reps = values.size();
  vector<vector<double> > vals_a(n_reps), vals_b(n_reps);
  for (size_t r = 0; r < n_reps; ++r)
    extract_fractional_values(values[r], vals_a[r], vals_b[r]);

  const size_t n_vals = values[0].size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> ff_vals(n_vals), fb_vals(n_vals); // for estimating transitions
  vector<double> bf_vals(n_vals), bb_vals(n_vals);
  vector<double> fg_emit(n_vals), bg_emit(n_vals); // avoid recomp of emissions

  if (VERBOSE)
    report_param_header_for_verbose();

  double prev_total = -std::numeric_limits<double>::max();
  bool converged = false;
  for (size_t i = 0; i < max_iterations && !converged; ++i) {

    double p_fb_est = p_fb, p_bf_est = p_bf;

    const double total =
      single_iteration_rep(values, vals_a, vals_b, reset_points, forward, backward,
                           p_fb_est, p_bf_est, fg_distro, bg_distro,
                           ff_vals, fb_vals, bf_vals, bb_vals, fg_emit, bg_emit);

    if (VERBOSE) // reporting for first replicate
      report_params_for_verbose(i, p_fb_est, p_bf_est,
                                fg_distro[0], bg_distro[0], total, prev_total);

    // ADS: removing the check based on expected log likelihood from
    // forward/backward as these seem to have some problem...
    converged = ((get_delta(p_fb_est, p_fb) < tolerance) &&
		 (get_delta(p_bf_est, p_bf) < tolerance));
    // converged = (get_delta(prev_total, total) < tolerance);

    if (converged) {
      if (VERBOSE)
        cerr << "CONVERGED" << endl;
    }
    else {
      p_fb = p_fb_est;
      p_bf = p_bf_est;
      prev_total = total;
    }
  }
  return prev_total;
}


void
TwoStateHMM::PosteriorScores(const vector<vector<pair<double, double> > > &values,
                             const vector<size_t> &reset_points,
                             const double p_fb, const double p_bf,
                             const vector<betabin> &fg_distro,
                             const vector<betabin> &bg_distro,
                             const bool fg_class,
                             vector<double> &posteriors) const {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  const size_t n_vals = values[0].size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> fg_emit(n_vals), bg_emit(n_vals);
  get_emissions_rep(values, fg_emit, fg_distro);
  get_emissions_rep(values, bg_emit, bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    total_loglik += score;
  }

  get_posteriors(forward, backward, posteriors);
  if (!fg_class)
    one_minus(begin(posteriors), end(posteriors), begin(posteriors));
}


double
TwoStateHMM::PosteriorDecoding(const vector<vector<pair<double, double> > > &values,
                               const vector<size_t> &reset_points,
                               const double p_fb, const double p_bf,
                               const vector<betabin> &fg_distro,
                               const vector<betabin> &bg_distro,
                               vector<bool> &classes,
                               vector<double> &posteriors) const {

  const double lp_sf = log(p_bf/(p_bf + p_fb));
  const double lp_sb = log(p_fb/(p_bf + p_fb));
  const double lp_ff = log(1.0 - p_fb);
  const double lp_fb = log(p_fb);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(1.0 - p_bf);

  const size_t n_vals = values[0].size();
  vector<pair<double, double> > forward(n_vals, make_pair(0.0, 0.0));
  vector<pair<double, double> > backward(n_vals, make_pair(0.0, 0.0));

  vector<double> fg_emit(n_vals), bg_emit(n_vals);
  get_emissions_rep(values, fg_emit, fg_distro);
  get_emissions_rep(values, bg_emit, bg_distro);

  double total_loglik = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const double score =
      forward_algorithm(reset_points[i], reset_points[i + 1],
                        lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                        fg_emit, bg_emit, forward);

    const double backward_score =
      backward_algorithm(reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb, lp_ff, lp_fb, lp_bf, lp_bb,
                         fg_emit, bg_emit, backward);

    assert(fabs(score - backward_score)/max(score, backward_score) < tolerance);

    total_loglik += score;
  }

  get_posteriors(forward, backward, posteriors);

  classes.resize(n_vals);
  for (size_t i = 0; i < n_vals; ++i)
    classes[i] = (posteriors[i] > 0.5);

  return total_loglik;
}
