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

#include "TwoStateCTHMM.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>
#include <iomanip>

#include <unordered_map>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;
using std::isfinite;
using std::numeric_limits;

struct prob_mat {
  double bb;
  double bf;
  double fb;
  double ff;
  prob_mat() {}
  prob_mat(const double mu0, const double mu1, const double t) {
    const double mu_sum = mu0 + mu1;
    const double frac = mu0/mu_sum;
    const double c = exp(-mu_sum*t);
    ff = (1.0 - frac)*c + frac;
    fb = 1.0 - ff;
    bb = frac*c + (1.0 - frac);
    bf = 1.0 - bb;
  }
  void make_logs() {
    bb = log(bb);
    bf = log(bf);
    fb = log(fb);
    ff = log(ff);
  }
};

std::ostream &
operator<<(std::ostream &os, const prob_mat &pm) {
  return os << pm.bb << '\t' << pm.bf << '\t'
            << pm.fb << '\t' << pm.ff;
}

struct two_by_two {
  double bb;
  double bf;
  double fb;
  double ff;
  two_by_two() : bb(0.0), bf(0.0), fb(0.0), ff(0.0) {}
  double total() const {return bb + bf + fb + ff;}
};

std::ostream &
operator<<(std::ostream &os, const two_by_two &x) {
  return os << x.bb << '\t' << x.bf << '\t'
            << x.fb << '\t' << x.ff;
}

struct betabin {
  betabin(const double a, const double b) :
    alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}
  double operator()(const pair<double, double> &val) const;
  void fit(const vector<double> &vals_a,
           const vector<double> &vals_b,
           const vector<double> &p);
  string tostring() const;
  double alpha;
  double beta;
  double lnbeta_helper;
};

string
betabin::tostring() const {
  std::ostringstream os;
  os << setprecision(4) << alpha << " " << setprecision(4) << beta;
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

static double
movement(const double curr, const double prev) {
  return std::abs(curr - prev)/max(std::abs(curr), std::abs(prev));
}

void
betabin::fit(const vector<double> &vals_a, const vector<double> &vals_b,
             const vector<double> &p) {
  const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
  const double alpha_rhs = inner_product(vals_a.begin(), vals_a.end(),
                                         p.begin(), 0.0)/p_total;
  const double beta_rhs = inner_product(vals_b.begin(), vals_b.end(),
                                        p.begin(), 0.0)/p_total;
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
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


inline double
TwoStateCTHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


double
TwoStateCTHMM::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
  const vector<double>::const_iterator x =
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
    }
  }
  return max_val + log(sum);
}


double
TwoStateCTHMM::forward_algorithm(const vector<uint32_t> &pos,
                                 const vector<pair<double, double> > &vals,
                                 const size_t start, const size_t end,
                                 const double lp_sf, const double lp_sb,
                                 const vector<prob_mat> &lm,
                                 const double lp_ft, const double lp_bt,
                                 const betabin &fg_distro,
                                 const betabin &bg_distro,
                                 vector<pair<double, double> > &f) const {
  f[start].first = fg_distro(vals[start]) + lp_sf;
  f[start].second = bg_distro(vals[start]) + lp_sb;
  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    f[i].first = (fg_distro(vals[i]) +
                  log_sum_log(f[k].first + lm[pos[i]].ff,
                              f[k].second + lm[pos[i]].bf));
    f[i].second = (bg_distro(vals[i]) +
                   log_sum_log(f[k].first + lm[pos[i]].fb,
                               f[k].second + lm[pos[i]].bb));
  }
  return log_sum_log(f[end - 1].first + lp_ft, f[end - 1].second + lp_bt);
}


double
TwoStateCTHMM::backward_algorithm(const vector<uint32_t> &pos,
                                  const vector<pair<double, double> > &vals,
                                  const size_t start, const size_t end,
                                  const double lp_sf, const double lp_sb,
                                  const vector<prob_mat> &lm,
                                  const double lp_ft, const double lp_bt,
                                  const betabin &fg_distro,
                                  const betabin &bg_distro,
                                  vector<pair<double, double> > &b) const {
  b[end - 1].first = lp_ft;
  b[end - 1].second = lp_bt;
  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    const double fg_a = fg_distro(vals[k]) + b[k].first;
    const double bg_a = bg_distro(vals[k]) + b[k].second;
    b[i].first = log_sum_log(fg_a + lm[pos[k]].ff, bg_a + lm[pos[k]].fb);
    b[i].second = log_sum_log(fg_a + lm[pos[k]].bf, bg_a + lm[pos[k]].bb);
  }
  return log_sum_log(b[start].first + fg_distro(vals[start]) + lp_sf,
                     b[start].second + bg_distro(vals[start]) + lp_sb);
}


void
TwoStateCTHMM::estimate_emissions(const vector<pair<double, double> > &f,
                                  const vector<pair<double, double> > &b,
                                  vector<double> &fg_probs,
                                  vector<double> &bg_probs) const {
  for (size_t i = 0; i < b.size(); ++i) {
    const double fg = (f[i].first + b[i].first);
    const double bg = (f[i].second + b[i].second);
    const double denom = log_sum_log(fg, bg);
    fg_probs[i] = exp(fg - denom);
    bg_probs[i] = exp(bg - denom);
  }
}



void
TwoStateCTHMM::estimate_transitions(const vector<uint32_t> &pos,
                                    const vector<pair<double, double> > &vals,
                                    const size_t start, const size_t end,
                                    const vector<pair<double, double> > &f,
                                    const vector<pair<double, double> > &b,
                                    const double total,
                                    const betabin &fg_distro,
                                    const betabin &bg_distro,
                                    const vector<prob_mat> &lm,
                                    vector<double> &ff_vals,
                                    vector<double> &fb_vals,
                                    vector<double> &bf_vals,
                                    vector<double> &bb_vals) const {
  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    const double b_first = b[i].first + fg_distro(vals[i]) - total;
    const double b_second = b[i].second + bg_distro(vals[i]) - total;

    const double ff = f[k].first;
    const double bb = f[k].second;

    ff_vals[k] = ff + lm[pos[i]].ff + b_first;
    fb_vals[k] = ff + lm[pos[i]].fb + b_second;

    bf_vals[k] = bb + lm[pos[i]].bf + b_first;
    bb_vals[k] = bb + lm[pos[i]].bb + b_second;
  }
}



static void
count_trans_by_dist(const size_t max_dist, const vector<uint32_t> &pos,
                    const vector<double> &bb, const vector<double> &bf,
                    const vector<double> &fb, const vector<double> &ff,
                    vector<two_by_two> &r) {

  r.clear();
  r.resize(max_dist);

  const size_t n_site_pairs = bb.size();
  for (size_t i = 0; i < n_site_pairs; ++i)
    if (pos[i] < max_dist) {
      r[pos[i]].bb += exp(bb[i]);
      r[pos[i]].bf += exp(bf[i]);
      r[pos[i]].fb += exp(fb[i]);
      r[pos[i]].ff += exp(ff[i]);
    }
}


static void
expectation_J(const double mu0, const double mu1, const double T,
              double &J_00_0, double &J_01_0,
              double &J_10_0, double &J_11_0,
              double &J_00_1, double &J_01_1,
              double &J_10_1, double &J_11_1) {

  const double s = mu0 + mu1;   assert(finite(s));
  const double p = mu0*mu1;     assert(finite(p));
  const double d = mu1 - mu0;   assert(finite(d));
  const double e = exp(-s*T);   assert(finite(e));

  const double C1 = d*(1 - e)/s;                        assert(finite(C1));
  J_00_0 = p*(T*(mu1 - mu0*e) - C1)/(s*(mu1 + mu0*e));  assert(finite(J_00_0));
  J_00_1 = J_00_0;                                      assert(finite(J_00_1));

  J_11_0 = p*(T*(mu0 - mu1*e) + C1)/(s*(mu0 + mu1*e));  assert(finite(J_11_0));
  J_11_1 = J_11_0;                                      assert(finite(J_11_1));

  const double C2 = p*T*(1 + e)/(s*(1 - e));
  const double C3 = (mu0*mu0 + mu1*mu1)/(s*s);
  const double C4 = (2*p)/(s*s);

  J_01_0 = C2 + C3;       assert(finite(J_01_0));
  J_10_1 = J_01_0;        assert(finite(J_10_1));

  J_01_1 = C2 - C4;       assert(finite(J_01_1));
  J_10_0 = J_01_1;        assert(finite(J_10_0));
}


static void
expectation_D(const double mu0, const double mu1, const double T,
              double &D_00_0, double &D_01_0,
              double &D_10_0, double &D_11_0,
              double &D_00_1, double &D_01_1,
              double &D_10_1, double &D_11_1) {

  const double mu00 = mu0*mu0;
  const double mu11 = mu1*mu1;
  const double s = mu0 + mu1;
  const double p = mu0*mu1;
  const double e = exp(-s*T);

  const double C1 = 2*p*(1 - e)/s;
  D_00_0 = ((mu11 + mu00*e)*T + C1)/(s*(mu1 + mu0*e));
  D_00_1 = T - D_00_0;

  D_11_1 = ((mu00 + mu11*e)*T + C1)/(s*(mu0 + mu1*e));
  D_11_0 = T - D_11_1;

  const double C3 = (p - mu00)*(1 - e)/s;
  D_01_1 = ((mu00 - p*e)*T + C3)/(s*(mu0 - mu0*e));
  D_01_0 = T - D_01_1;

  const double C4 = (p - mu11)*(1 - e)/s;
  D_10_0 = ((mu11 - p*e)*T + C4)/(s*(mu1 - mu1*e));
  D_10_1 = T - D_10_0;
}


static void
estimate_rates(const vector<two_by_two> &c, double &mu0, double &mu1) {

  double J0 = 0.0, J1 = 0.0;
  double D0 = 0.0, D1 = 0.0;
  for (size_t i = 2; i < c.size(); ++i) {
    if (c[i].total() > numeric_limits<double>::min()) {
      double J_00_0 = 0.0, J_01_0 = 0.0, J_10_0 = 0.0, J_11_0 = 0.0;
      double J_00_1 = 0.0, J_01_1 = 0.0, J_10_1 = 0.0, J_11_1 = 0.0;
      expectation_J(mu0, mu1, i,
                    J_00_0, J_01_0, J_10_0, J_11_0,
                    J_00_1, J_01_1, J_10_1, J_11_1);

      assert(finite(J_00_0)); assert(finite(J_01_0));
      assert(finite(J_10_0)); assert(finite(J_11_0));

      assert(finite(J_00_1)); assert(finite(J_01_1));
      assert(finite(J_10_1)); assert(finite(J_11_1));

      J0 += (c[i].bb*J_00_0 + c[i].bf*J_01_0 + c[i].fb*J_10_0 + c[i].ff*J_11_0);
      J1 += (c[i].bb*J_00_1 + c[i].bf*J_01_1 + c[i].fb*J_10_1 + c[i].ff*J_11_1);

      double D_00_0 = 0.0, D_01_0 = 0.0, D_10_0 = 0.0, D_11_0 = 0.0;
      double D_00_1 = 0.0, D_01_1 = 0.0, D_10_1 = 0.0, D_11_1 = 0.0;
      expectation_D(mu0, mu1, i,
                    D_00_0, D_01_0, D_10_0, D_11_0,
                    D_00_1, D_01_1, D_10_1, D_11_1);

      assert(finite(D_00_0)); assert(finite(D_01_0));
      assert(finite(D_10_0)); assert(finite(D_11_0));

      assert(finite(D_00_1)); assert(finite(D_01_1));
      assert(finite(D_10_1)); assert(finite(D_11_1));

      D0 += (c[i].bb*D_00_0 + c[i].bf*D_01_0 + c[i].fb*D_10_0 + c[i].ff*D_11_0);
      D1 += (c[i].bb*D_00_1 + c[i].bf*D_01_1 + c[i].fb*D_10_1 + c[i].ff*D_11_1);
    }
  }
  mu0 = J0/D0;
  mu1 = J1/D1;
}


double
TwoStateCTHMM::single_iteration(const vector<uint32_t> &pos,
                                const vector<pair<double, double> > &values,
                                const vector<double> &vals_a,
                                const vector<double> &vals_b,
                                const vector<size_t> &reset_points,
                                vector<pair<double, double> > &forward,
                                vector<pair<double, double> > &backward,
                                double &p_sf, double &p_sb,
                                double &mu0, double &mu1,
                                double &p_ft, double &p_bt,
                                betabin &fg_distro,
                                betabin &bg_distro) const {

  vector<double> log_fg_expected;
  vector<double> log_bg_expected;

  double total_score = 0;

  const double lp_sb = log(p_sb);
  const double lp_sf = log(p_sf);
  const double lp_bt = log(p_bt);
  const double lp_ft = log(p_ft);

  vector<prob_mat> lm(desert_size);
  for (size_t i = 0; i < lm.size(); ++i) {
    lm[i] = prob_mat(mu0, mu1, i);
    lm[i].make_logs();
  }

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ft) && isfinite(lp_bt));

  // for estimating transitions
  vector<double> ff_vals(values.size(), 0);
  vector<double> fb_vals(values.size(), 0);
  vector<double> bf_vals(values.size(), 0);
  vector<double> bb_vals(values.size(), 0);

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(pos, values,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lm,
                                           lp_ft, lp_bt,
                                           fg_distro, bg_distro, forward);
    const double backward_score =
      backward_algorithm(pos, values,
                         reset_points[i], reset_points[i + 1],
                         lp_sf, lp_sb,
                         lm,
                         lp_ft, lp_bt,
                         fg_distro, bg_distro, backward);

    estimate_transitions(pos, values,
                         reset_points[i], reset_points[i + 1],
                         forward, backward,
                         score,
                         fg_distro, bg_distro,
                         lm,
                         ff_vals, fb_vals,
                         bf_vals, bb_vals);

    total_score += score;
  }

  // Subtracting 1 from the limit of the summation because the final
  // term has no meaning since there is no transition to be counted
  // from the final observation (they all must go to terminal state)

  vector<two_by_two> trans_counts(desert_size);
  count_trans_by_dist(desert_size, pos,
                      bb_vals, bf_vals, fb_vals, ff_vals, trans_counts);

  cerr << "trans counts:\n";
  for (size_t i = 0; i < 30; ++i) {
    cerr << i << '\t' << trans_counts[i] << endl;
  }
  cerr << endl;

  double prev_mu0 = 0.0, prev_mu1 = 0.0;

  for (size_t i = 0; i < 100 &&
         (abs(prev_mu0 - mu0) > numeric_limits<double>::min() ||
          abs(prev_mu1 - mu1) > numeric_limits<double>::min()); ++i) {
    prev_mu0 = mu0;
    prev_mu1 = mu1;
    estimate_rates(trans_counts, mu0, mu1);
  }

  const double denom = (1.0/mu0 + 1.0/mu1);
  p_sf = 1.0/(mu0*denom);
  p_sb = 1.0/(mu1*denom);

  p_ft = p_sf;
  p_bt = p_sb;

  // for estimating emissions
  vector<double> fg_probs(values.size());
  vector<double> bg_probs(values.size());
  estimate_emissions(forward, backward, fg_probs, bg_probs);

  fg_distro.fit(vals_a, vals_b, fg_probs);
  bg_distro.fit(vals_a, vals_b, bg_probs);

  return total_score;
}



double
TwoStateCTHMM::BaumWelchTraining(const vector<uint32_t> &pos,
                                 const std::vector<pair<double, double> > &values,
                                 const std::vector<size_t> &reset_points,
                                 vector<double> &start_trans,
                                 double &mu0, double &mu1,
                                 vector<double> &end_trans,
                                 double &fg_alpha, double &fg_beta,
                                 double &bg_alpha, double &bg_beta) const {

  betabin fg_distro(fg_alpha, fg_beta);
  betabin bg_distro(bg_alpha, bg_beta);

  assert(start_trans.size() >= 2 && end_trans.size() >= 2);

  const double score = BaumWelchTraining(pos, values, reset_points,
                                         start_trans[0], start_trans[1],
                                         mu0, mu1,
                                         end_trans[0], end_trans[1],
                                         fg_distro, bg_distro);

  fg_alpha = fg_distro.alpha;
  fg_beta = fg_distro.beta;

  bg_alpha = bg_distro.alpha;
  bg_beta = bg_distro.beta;

  return score;
}



double
TwoStateCTHMM::BaumWelchTraining(const vector<uint32_t> &pos,
                                 const vector<pair<double, double> > &values,
                                 const vector<size_t> &reset_points,
                                 double &p_sf, double &p_sb,
                                 double &mu0, double &mu1,
                                 double &p_ft, double &p_bt,
                                 betabin &fg_distro,
                                 betabin &bg_distro) const {

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  if (VERBOSE)
    cerr << setw(5)  << "ITR"
         << setw(10) << "F size"
         << setw(10) << "B size"
         << setw(18) << "F PARAMS"
         << setw(18) << "B PARAMS"
         << setw(14) << "DELTA"
         << endl;

  double prev_total = -std::numeric_limits<double>::max();

  vector<double> vals_a(values.size()), vals_b(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    vals_a[i] = // values[i].first;
      log(min(max(values[i].first/(values[i].first + values[i].second), 1e-2),
              1.0 - 1e-2));
    vals_b[i] = // values[i].second;
      log(1 - min(max(values[i].first/(values[i].first + values[i].second), 1e-2),
                  1.0 - 1e-2));
  }

  for (size_t i = 0; i < max_iterations; ++i) {
    double p_sf_est = p_sf;
    double p_sb_est = p_sb;
    double mu0_est = mu0;
    double mu1_est = mu1;
    double p_ft_est = p_ft;
    double p_bt_est = p_bt;

    double total = single_iteration(pos, values,
                                    vals_a, vals_b,
                                    reset_points,
                                    forward, backward,
                                    p_sf_est, p_sb_est,
                                    mu0_est, mu1_est,
                                    p_ft_est, p_bt_est,
                                    fg_distro, bg_distro);

    if (VERBOSE) {
      cerr << setw(5) << i + 1
           << setw(10) << 1.0/mu0_est
           << setw(10) << 1.0/mu1_est
           << setw(18) << fg_distro.tostring()
           << setw(18) << bg_distro.tostring()
           << setw(14) << total
           << setw(14) << prev_total
           << setw(14) << (total - prev_total)/std::fabs(total)
           << endl;
    }
    if ((total - prev_total) < tolerance) {
      if (VERBOSE)
        cerr << "CONVERGED" << endl << endl;
      break;
    }

    p_sf = p_sf_est;
    p_sb = p_sb_est;
    mu0 = mu0_est;
    mu1 = mu1_est;
    p_ft = p_ft_est;
    p_bt = p_bt_est;

    prev_total = total;
  }
  return prev_total;
}


double
TwoStateCTHMM::PosteriorDecoding(const vector<uint32_t> &pos,
                                 const vector<pair<double, double> > &values,
                                 const vector<size_t> &reset_points,
                                 const vector<double> &start_trans,
                                 const double mu0, const double mu1,
                                 const vector<double> &end_trans,
                                 const double fg_alpha, const double fg_beta,
                                 const double bg_alpha, const double bg_beta,
                                 vector<bool> &classes,
                                 vector<double> &llr_scores) const {


  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);

  return PosteriorDecoding(pos, values, reset_points,
                           start_trans[0], start_trans[1],
                           mu0, mu1,
                           end_trans[0], end_trans[1],
                           fg_distro, bg_distro, classes, llr_scores);
}


double
TwoStateCTHMM::PosteriorDecoding(const vector<uint32_t> &pos,
                                 const vector<pair<double, double> > &values,
                                 const vector<size_t> &reset_points,
                                 const double p_sf, const double p_sb,
                                 const double mu0, const double mu1,
                                 const double p_ft, const double p_bt,
                                 const betabin &fg_distro,
                                 const betabin &bg_distro,
                                 vector<bool> &classes,
                                 vector<double> &llr_scores) const {

  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ft = log(p_ft);
  const double lp_bt = log(p_bt);

  vector<prob_mat> lm(desert_size);
  for (size_t i = 0; i < lm.size(); ++i) {
    lm[i] = prob_mat(mu0, mu1, i);
    lm[i].make_logs();
  }

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ft) && isfinite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(pos, values,
                                           reset_points[i], reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lm,
                                           lp_ft, lp_bt,
                                           fg_distro, bg_distro,
                                           forward);

    const double backward_score =
      backward_algorithm(pos,
                         values,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sb,
                         lm,
                         lp_ft, lp_bt,
                         fg_distro, bg_distro,
                         backward);

    total_score += score;
  }

  classes.resize(values.size());

  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    classes[i] = static_cast<bool>(fg_state > bg_state);
    llr_scores[i] = exp(fg_state - log_sum_log(fg_state, bg_state));
  }

  return total_score;
}
