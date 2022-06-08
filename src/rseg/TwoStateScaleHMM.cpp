/*
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "TwoStateScaleHMM.hpp"

#include <iomanip>
#include <cmath>
#include <limits>
#include <memory>

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::auto_ptr;
using std::isfinite;


inline double
TwoStateScaleHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {
    return q;
  }
  else if (q == 0) {
    return p;
  }
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


double
TwoStateScaleHMM::forward_algorithm(const vector<double> &vals,
                               const vector<double> &scales,
                               const size_t start, const size_t end,
                               const double lp_sf, const double lp_sb,
                               const double lp_ff, const double lp_fb,
                               const double lp_ft,
                               const double lp_bf, const double lp_bb,
                               const double lp_bt,
                               const Distro &fg_distro,
                               const Distro &bg_distro,
                               vector<pair<double, double> > &f) const {

    f[start].first = fg_distro.log_likelihood(vals[start], scales[start]) + lp_sf;
    f[start].second = bg_distro.log_likelihood(vals[start], scales[start]) + lp_sb;

  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    f[i].first = (fg_distro.log_likelihood(vals[i], scales[i]) +
                  log_sum_log(f[k].first + lp_ff, f[k].second + lp_bf));
    f[i].second = (bg_distro.log_likelihood(vals[i], scales[i]) +
                   log_sum_log(f[k].first + lp_fb, f[k].second + lp_bb));
  }
  return log_sum_log(f[end - 1].first + lp_ft, f[end - 1].second + lp_bt);
}


double
TwoStateScaleHMM::backward_algorithm(const vector<double> &vals,
                                const vector<double> &scales,
                                const size_t start, const size_t end,
                                const double lp_sf, const double lp_sb,
                                const double lp_ff, const double lp_fb,
                                const double lp_ft,
                                const double lp_bf, const double lp_bb,
                                const double lp_bt,
                                const Distro &fg_distro,
                                const Distro &bg_distro,
                                vector<pair<double, double> > &b) const {

  b[end - 1].first = lp_ft;
  b[end - 1].second = lp_bt;

  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    const double fg_a = fg_distro.log_likelihood(vals[k], scales[k]) + b[k].first;
    const double bg_a = bg_distro.log_likelihood(vals[k], scales[k]) + b[k].second;
    b[i].first = log_sum_log(fg_a + lp_ff, bg_a + lp_fb);
    b[i].second = log_sum_log(fg_a + lp_bf, bg_a + lp_bb);
  }
  return log_sum_log(b[start].first +
                     fg_distro.log_likelihood(vals[start], scales[start]) +
                     lp_sf,
                     b[start].second +
                     bg_distro.log_likelihood(vals[start], scales[start]) +
                     lp_sb);
}


double
TwoStateScaleHMM::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
  const vector<double>::const_iterator x =
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(isfinite(sum));
#endif
    }
  }
  return max_val + log(sum);
}


void
TwoStateScaleHMM::estimate_emissions(const vector<pair<double, double> > &f,
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
TwoStateScaleHMM::estimate_transitions(const vector<double> &vals,
                                  const vector<double> &scales,
                                  const size_t start, const size_t end,
                                  const vector<pair<double, double> > &f,
                                  const vector<pair<double, double> > &b,
                                  const double total,
                                  const Distro &fg_distro,
                                  const Distro &bg_distro,
                                  const double lp_ff, const double lp_fb,
                                  const double lp_bf, const double lp_bb,
                                  const double lp_ft, const double lp_bt,
                                  vector<double> &ff_vals,
                                  vector<double> &fb_vals,
                                  vector<double> &bf_vals,
                                  vector<double> &bb_vals) const {

    for (size_t i = start + 1; i < end; ++i)
    {
        const size_t k = i - 1;

        const double lp_fg = fg_distro.log_likelihood(vals[i], scales[i]) - total;
        const double lp_bg = bg_distro.log_likelihood(vals[i], scales[i]) - total;

        ff_vals[k] = f[k].first + lp_ff + lp_fg + b[i].first;
        fb_vals[k] = f[k].first + lp_fb + lp_bg + b[i].second;

        bf_vals[k] = f[k].second + lp_bf + lp_fg + b[i].first;
        bb_vals[k] = f[k].second + lp_bb + lp_bg + b[i].second;
    }
}

double
TwoStateScaleHMM::single_iteration(const vector<double> &values,
                              const vector<double> &scales,
                              const vector<size_t> &reset_points,
                              vector<pair<double, double> > &forward,
                              vector<pair<double, double> > &backward,
                              double &p_sf, double &p_sb,
                              double &p_ff, double &p_fb, double &p_ft,
                              double &p_bf, double &p_bb, double &p_bt,
                              Distro &fg_distro,
                              Distro &bg_distro) const {

  vector<double> log_fg_expected;
  vector<double> log_bg_expected;

  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ff) && isfinite(lp_fb) && isfinite(lp_ft) &&
         isfinite(lp_bf) && isfinite(lp_bb) && isfinite(lp_bt));

  // for estimating transitions
  vector<double> ff_vals(values.size(), 0);
  vector<double> fb_vals(values.size(), 0);
  vector<double> bf_vals(values.size(), 0);
  vector<double> bb_vals(values.size(), 0);

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lp_ff, lp_fb, lp_ft,
                                           lp_bf, lp_bb, lp_bt,
                                           fg_distro, bg_distro, forward);
    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sb,
                         lp_ff, lp_fb, lp_ft,
                         lp_bf, lp_bb, lp_bt,
                         fg_distro, bg_distro, backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    estimate_transitions(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         forward, backward,
                         score,
                         fg_distro, bg_distro,
                         lp_ff, lp_fb, lp_bf,
                         lp_bb, lp_ft, lp_bt,
                         ff_vals, fb_vals,
                         bf_vals, bb_vals);

    total_score += score;
  }

  // Subtracting 1 from the limit of the summation because the final
  // term has no meaning since there is no transition to be counted
  // from the final observation (they all must go to terminal state)
  // SQ: note ff_vals[reset_points[i]] is always euqal to 0 (exp() ==
  // 1), i.e. the start of each region does not contribute to
  // transition emission estimate
  const double p_ff_new_estimate = exp(log_sum_log_vec(ff_vals, values.size()))
      - reset_points.size() + 1;
  const double p_fb_new_estimate = exp(log_sum_log_vec(fb_vals, values.size()))
      - reset_points.size() + 1;
  const double p_bf_new_estimate = exp(log_sum_log_vec(bf_vals, values.size()))
      - reset_points.size() + 1;
  const double p_bb_new_estimate = exp(log_sum_log_vec(bb_vals, values.size()))
      - reset_points.size() + 1;
//   const double p_ff_new_estimate = exp(log_sum_log_vec(ff_vals, values.size() - 1));
//   const double p_fb_new_estimate = exp(log_sum_log_vec(fb_vals, values.size() - 1));
//   const double p_bf_new_estimate = exp(log_sum_log_vec(bf_vals, values.size() - 1));
//   const double p_bb_new_estimate = exp(log_sum_log_vec(bb_vals, values.size() - 1));


  double denom = (p_ff_new_estimate + p_fb_new_estimate);
  p_ff = p_ff_new_estimate/denom - p_ft/2.0;
  p_fb = p_fb_new_estimate/denom - p_ft/2.0;

  if (p_ff < MIN_PROB) {
    if (DEBUG)
      cerr << "p_ff < MIN_PROB" << endl;
    p_ff = MIN_PROB;
  }

  if (p_fb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_fb < MIN_PROB" << endl;
    p_fb = MIN_PROB;
  }

  denom = (p_bf_new_estimate + p_bb_new_estimate);
  p_bf = p_bf_new_estimate/denom - p_bt/2.0;
  p_bb = p_bb_new_estimate/denom - p_bt/2.0;

  if (p_bf < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bf < MIN_PROB" << endl;
    p_bf = MIN_PROB;
  }

  if (p_bb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bb < MIN_PROB" << endl;
    p_bb = MIN_PROB;
  }

  // for estimating emissions
  vector<double> fg_probs(values.size());
  vector<double> bg_probs(values.size());
  estimate_emissions(forward, backward, fg_probs, bg_probs);

  fg_distro.estimate_params_ml(values, scales, fg_probs);
  bg_distro.estimate_params_ml(values, scales, bg_probs);

  return total_score;
}


double
TwoStateScaleHMM::BaumWelchTraining(const std::vector<double> &values,
                               const std::vector<double> &scales,
                               const std::vector<size_t> &reset_points,
                               vector<double> &start_trans,
                               vector<vector<double> > &trans,
                               vector<double> &end_trans,
                               Distro &fg_distro,
                               Distro &bg_distro) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return BaumWelchTraining(values, scales, reset_points,
                           start_trans[0], start_trans[1],
                           trans[0][0], trans[0][1], end_trans[0],
                           trans[1][0], trans[1][1], end_trans[1],
                           fg_distro, bg_distro);
}

double
TwoStateScaleHMM::BaumWelchTraining(const vector<double> &values,
                               const vector<double> &scales,
                               const vector<size_t> &reset_points,
                               double &p_sf, double &p_sb,
                               double &p_ff, double &p_fb, double &p_ft,
                               double &p_bf, double &p_bb, double &p_bt,
                               Distro &fg_distro,
                               Distro &bg_distro) const {

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  if (VERBOSE)
    cout << setw(5)  << "ITR"
         << setw(10) << "F size"
         << setw(10) << "B size"
         << setw(18) << "F PARAMS"
         << setw(18) << "B PARAMS"
         << setw(14) << "DELTA"
         << endl;

  double prev_total = -std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_iterations; ++i) {

    double p_sf_est = p_sf;
    double p_sb_est = p_sb;
    double p_ff_est = p_ff;
    double p_fb_est = p_fb;
    double p_bf_est = p_bf;
    double p_bb_est = p_bb;
    double p_ft_est = p_ft;
    double p_bt_est = p_bt;

    double total = single_iteration(values, scales,
                                    reset_points,
                                    forward, backward,
                                    p_sf_est, p_sb_est,
                                    p_ff_est, p_fb_est, p_ft_est,
                                    p_bf_est, p_bb_est, p_bt_est,
                                    fg_distro, bg_distro);

    if ((prev_total - total)/prev_total < tolerance) {
      if (VERBOSE)
        cout << "CONVERGED" << endl << endl;
      break;
    }

    if (VERBOSE) {
      cout << setw(5) << i + 1
           << setw(10) << 1/p_fb_est
           << setw(10) << 1/p_bf_est
           << setw(18) << fg_distro.tostring()
           << setw(18) << bg_distro.tostring()
           << setw(14) << (prev_total - total)/prev_total
           << endl;
    }

    p_sf = p_sf_est;
    p_sb = p_sb_est;
    p_ff = p_ff_est;
    p_fb = p_fb_est;
    p_bf = p_bf_est;
    p_bb = p_bb_est;
    p_ft = p_ft_est;
    p_bt = p_bt_est;

    prev_total = total;
  }
  return prev_total;
}

void
TwoStateScaleHMM::PosteriorScores(const vector<double> &values,
                             const vector<double> &scales,
                             const vector<size_t> &reset_points,
                             const vector<double> &start_trans,
                             const vector<vector<double> > &trans,
                             const vector<double> &end_trans,
                             const Distro &fg_distro,
                             const Distro &bg_distro,
                             const vector<bool> &classes,
                             vector<double> &llr_scores) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return PosteriorScores(values, scales, reset_points,
                         start_trans[0], start_trans[1],
                         trans[0][0], trans[0][1], end_trans[0],
                         trans[1][0], trans[1][1], end_trans[1],
                         fg_distro, bg_distro, classes, llr_scores);
}


void
TwoStateScaleHMM::PosteriorScores(const vector<double> &values,
                             const vector<double> &scales,
                             const vector<size_t> &reset_points,
                             double p_sf, double p_sb,
                             double p_ff, double p_fb, double p_ft,
                             double p_bf, double p_bb, double p_bt,
                             const Distro &fg_distro,
                             const Distro &bg_distro,
                             const vector<bool> &classes,
                             vector<double> &llr_scores) const {

  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ff) && isfinite(lp_fb) && isfinite(lp_ft) &&
         isfinite(lp_bf) && isfinite(lp_bb) && isfinite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lp_ff, lp_fb, lp_ft,
                                           lp_bf, lp_bb, lp_bt,
                                           fg_distro, bg_distro, forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sb,
                         lp_ff, lp_fb, lp_ft,
                         lp_bf, lp_bb, lp_bt,
                         fg_distro, bg_distro, backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    const double denom = log_sum_log(fg_state, bg_state);
    llr_scores[i] =
        classes[i] ? exp(fg_state - denom) : exp(bg_state - denom);
  }
}


void
TwoStateScaleHMM::PosteriorScores(const vector<double> &values,
                             const vector<double> &scales,
                             const vector<size_t> &reset_points,
                             const vector<double> &start_trans,
                             const vector<vector<double> > &trans,
                             const vector<double> &end_trans,
                             const Distro &fg_distro,
                             const Distro &bg_distro,
                             const bool fg_class,
                             vector<double> &llr_scores) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return PosteriorScores(values, scales, reset_points,
                         start_trans[0], start_trans[1],
                         trans[0][0], trans[0][1], end_trans[0],
                         trans[1][0], trans[1][1], end_trans[1],
                         fg_distro, bg_distro, fg_class, llr_scores);
}


void
TwoStateScaleHMM::PosteriorScores(const vector<double> &values,
                             const vector<double> &scales,
                             const vector<size_t> &reset_points,
                             double p_sf, double p_sb,
                             double p_ff, double p_fb, double p_ft,
                             double p_bf, double p_bb, double p_bt,
                             const Distro &fg_distro,
                             const Distro &bg_distro,
                             const bool fg_class,
                             vector<double> &llr_scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ff) && isfinite(lp_fb) && isfinite(lp_ft) &&
         isfinite(lp_bf) && isfinite(lp_bb) && isfinite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lp_ff, lp_fb, lp_ft,
                                           lp_bf, lp_bb, lp_bt,
                                           fg_distro, bg_distro, forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sb,
                         lp_ff, lp_fb, lp_ft,
                         lp_bf, lp_bb, lp_bt,
                         fg_distro, bg_distro, backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    const double denom = log_sum_log(fg_state, bg_state);
    llr_scores[i] =
        fg_class ? exp(fg_state - denom) : exp(bg_state - denom);
  }
}


void
TwoStateScaleHMM::TransitionPosteriors(const vector<double> &values,
                                  const vector<double> &scales,
                                  const vector<size_t> &reset_points,
                                  const vector<double> &start_trans,
                                  const vector<vector<double> > &trans,
                                  const vector<double> &end_trans,
                                  const Distro &fg_distro,
                                  const Distro &bg_distro,
                                  const size_t transition,
                                  vector<double> &llr_scores) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return TransitionPosteriors(values, scales, reset_points,
                              start_trans[0], start_trans[1],
                              trans[0][0], trans[0][1], end_trans[0],
                              trans[1][0], trans[1][1], end_trans[1],
                              fg_distro, bg_distro, transition, llr_scores);
}


void
TwoStateScaleHMM::TransitionPosteriors(const vector<double> &values,
                                  const std::vector<double> &scales,
                                  const vector<size_t> &reset_points,
                                  double p_sf, double p_sb,
                                  double p_ff, double p_fb, double p_ft,
                                  double p_bf, double p_bb, double p_bt,
                                  const Distro &fg_distro,
                                  const Distro &bg_distro,
                                  const size_t transition,
                                  vector<double> &scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ff) && isfinite(lp_fb) && isfinite(lp_ft) &&
         isfinite(lp_bf) && isfinite(lp_bb) && isfinite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lp_ff, lp_fb, lp_ft,
                                           lp_bf, lp_bb, lp_bt,
                                           fg_distro, bg_distro, forward);

    const double backward_score =
        backward_algorithm(
            values, scales,
            reset_points[i], reset_points[i + 1],
            lp_sf, lp_sb, lp_ff, lp_fb, lp_ft,
            lp_bf, lp_bb, lp_bt,
            fg_distro, bg_distro, backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  scores.resize(values.size());
  size_t j = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i == reset_points[j]) {
      ++j;
      scores[i] = 0;
    }
    else {
      const double fg_to_fg_state = forward[i - 1].first + lp_ff + // transition
        // emission for value i + 1
          fg_distro.log_likelihood(values[i], scales[i]) + backward[i].first;
      const double fg_to_bg_state = forward[i - 1].first + lp_fb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i].second;
      const double bg_to_fg_state = forward[i - 1].second + lp_bf +
        fg_distro.log_likelihood(values[i], scales[i]) + backward[i].first;
      const double bg_to_bg_state = forward[i - 1].second + lp_bb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i].second;
      const double denom = log_sum_log(log_sum_log(fg_to_fg_state, fg_to_bg_state),
                                       log_sum_log(bg_to_fg_state, bg_to_bg_state));
      double numerator = fg_to_fg_state;
      if (transition == 1)
        numerator = fg_to_bg_state;
      if (transition == 2)
        numerator = bg_to_fg_state;
      if (transition == 3)
        numerator = bg_to_bg_state;
      scores[i] = exp(numerator - denom);
    }
  }
}


double
TwoStateScaleHMM::PosteriorDecoding(const vector<double> &values,
                               const std::vector<double> &scales,
                               const vector<size_t> &reset_points,
                               const vector<double> &start_trans,
                               const vector<vector<double> > &trans,
                               const vector<double> &end_trans,
                               const Distro &fg_distro,
                               const Distro &bg_distro,
                               vector<bool> &classes,
                               vector<double> &llr_scores) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return PosteriorDecoding(values, scales, reset_points,
                           start_trans[0], start_trans[1],
                           trans[0][0], trans[0][1], end_trans[0],
                           trans[1][0], trans[1][1], end_trans[1],
                           fg_distro, bg_distro, classes, llr_scores);
}


double
TwoStateScaleHMM::PosteriorDecoding(const vector<double> &values,
                               const std::vector<double> &scales,
                               const vector<size_t> &reset_points,
                               double p_sf, double p_sb,
                               double p_ff, double p_fb, double p_ft,
                               double p_bf, double p_bb, double p_bt,
                               const Distro &fg_distro,
                               const Distro &bg_distro,
                               vector<bool> &classes,
                               vector<double> &llr_scores) const {

  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) && isfinite(lp_sb) &&
         isfinite(lp_ff) && isfinite(lp_fb) && isfinite(lp_ft) &&
         isfinite(lp_bf) && isfinite(lp_bb) && isfinite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sb,
                                           lp_ff, lp_fb, lp_ft,
                                           lp_bf, lp_bb, lp_bt,
                                           fg_distro, bg_distro, forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sb,
                         lp_ff, lp_fb, lp_ft,
                         lp_bf, lp_bb, lp_bt,
                         fg_distro, bg_distro, backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  classes.resize(values.size());
  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    const double denom = log_sum_log(fg_state, bg_state);
    classes[i] = static_cast<bool>(fg_state > bg_state);
    llr_scores[i] =
        classes[i] ? exp(fg_state - denom) : exp(bg_state - denom);
  }

  return total_score;
}



/*************************************************************
 *
 * Functions for Viterbi training and decoding.
 *
 *************************************************************/


double
TwoStateScaleHMM::ViterbiDecoding(const vector<double> &values,
                             const std::vector<double> &scales,
                             const vector<size_t> &reset_points,
                             const vector<double> &start_trans,
                             const vector<vector<double> > &trans,
                             const vector<double> &end_trans,
                             const Distro &fg_distro,
                             const Distro &bg_distro,
                             vector<bool> &classes) const {

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return ViterbiDecoding(values, scales, reset_points,
                         start_trans[0], start_trans[1],
                         trans[0][0], trans[0][1], end_trans[0],
                         trans[1][0], trans[1][1], end_trans[1],
                         fg_distro, bg_distro, classes);
}



double
TwoStateScaleHMM::ViterbiDecoding(const vector<double> &values,
                             const std::vector<double> &scales,
                             const vector<size_t> &reset_points,
                             double p_sf, double p_sb,
                             double p_ff, double p_fb, double p_ft,
                             double p_bf, double p_bb, double p_bt,
                             const Distro &fg_distro,
                             const Distro &bg_distro,
                             vector<bool> &ml_classes) const {

  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  // ml_classes = vector<bool>(values.size());
  double total = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const size_t start = reset_points[i];
      const size_t lim = reset_points[i + 1] - start;

      vector<pair<double, double> > v(lim, pair<double, double>(0, 0));
      vector<pair<size_t, size_t> > trace(lim, pair<size_t, size_t>(0, 0));

      v.front().first =
          lp_sf + fg_distro.log_likelihood(values[start], scales[start]);
      v.front().second =
          lp_sb + bg_distro.log_likelihood(values[start], scales[start]);

      for (size_t j = 1; j < lim; ++j) {

          const double ff = v[j - 1].first + lp_ff;
          const double bf = v[j - 1].second + lp_bf;
          const double fg_log_emmit =
              fg_distro.log_likelihood(values[start + j], scales[start + j]);
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
          const double bg_log_emmit =
              bg_distro.log_likelihood(values[start + j], scales[start + j]);
          if (fb > bb) {
              v[j].second = bg_log_emmit + fb;
              trace[j].second = 0;
          }
          else {
              v[j].second = bg_log_emmit + bb;
              trace[j].second = 1;
          }
      }
      v.back().first += lp_ft;
      v.back().second += lp_bt;

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

      reverse(inner_ml_classes.begin(), inner_ml_classes.end());
      ml_classes.insert(ml_classes.end(), inner_ml_classes.begin(),
                        inner_ml_classes.end());

      total += max(v.back().first, v.back().second);
  }

  return total;
}
