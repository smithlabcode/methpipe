/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Authors: Andrew D. Smith

  This file is part of methpipe.

  methpipe is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  methpipe is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with methpipe; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "TwoStateHMM_PMD.hpp"

// #pragma omp <rest of pragma>

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
using std::unique_ptr;

inline double
TwoStateHMM::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


double
TwoStateHMM::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
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
TwoStateHMM::estimate_emissions(const vector<pair<double, double> > &f,
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
TwoStateHMM::TransitionPosteriors_rep(
                   const vector<vector<pair<double, double> > > &values,
                   const vector<size_t> &reset_points,
                   const vector<double> &start_trans,
                   const vector<vector<double> > &trans,
                   const vector<double> &end_trans,
                   const vector<double> &fg_alpha,
                   const vector<double> &fg_beta,
                   const vector<double> &bg_alpha,
                   const vector<double> &bg_beta,
                   const vector<bool> &array_status,
                   const size_t transition,
                   vector<double> &llr_scores) const {

  vector<unique_ptr<EmissionDistribution> > fg_distro;
  vector<unique_ptr<EmissionDistribution> > bg_distro;

  size_t NREP = values.size();
  for (size_t i = 0; i < NREP; ++i) {
    if(array_status[i]) {
      fg_distro.emplace_back(new Beta(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new Beta(bg_alpha[i], bg_beta[i]));
    }
    else {
      fg_distro.emplace_back(new BetaBinomial(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new BetaBinomial(bg_alpha[i], bg_beta[i]));
    }
  }

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  return TransitionPosteriors_rep(values, reset_points,
                  start_trans[0], start_trans[1],
                  trans[0][0], trans[0][1], end_trans[0],
                  trans[1][0], trans[1][1], end_trans[1],
                  fg_distro, bg_distro, transition, array_status, llr_scores);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////   For multiple replicates       ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double
TwoStateHMM::forward_algorithm_rep(
            const vector<vector<pair<double, double> > > &vals,
                                    const size_t start, const size_t end,
                                    const double lp_sf, const double lp_sb,
                                    const double lp_ff, const double lp_fb,
                                    const double lp_ft, const double lp_bf,
                                    const double lp_bb, const double lp_bt,
                                    const vector<unique_ptr<EmissionDistribution> > &fg_distro,
                                    const vector<unique_ptr<EmissionDistribution> > &bg_distro,
                                    vector<pair<double, double> > &f,
            const vector<bool> &array_status) const {
  size_t NREP = vals.size();
  f[start].first = lp_sf;
  f[start].second = lp_sb;

  for (size_t r = 0; r < NREP; ++r) {
    if (vals[r][start].first + vals[r][start].second  >= 1) {
      f[start].first += (*fg_distro[r])(vals[r][start]);
      f[start].second += (*bg_distro[r])(vals[r][start]);
    }
  }
  for (size_t i = start + 1; i < end; ++i) {
    // assert(isfinite(fg_distro.log_likelihood(vals[i])));
    const size_t k = i - 1;
    f[i].first = log_sum_log(f[k].first + lp_ff, f[k].second + lp_bf);
    f[i].second = log_sum_log(f[k].first + lp_fb, f[k].second + lp_bb);

    for (size_t r = 0; r < NREP; ++r) {
      if(vals[r][i].first + vals[r][i].second >= 1) {
              f[i].first += (*fg_distro[r])(vals[r][i]);
              f[i].second += (*bg_distro[r])(vals[r][i]);
      }
    }
  }
  return log_sum_log(f[end - 1].first + lp_ft, f[end - 1].second + lp_bt);
}


double
TwoStateHMM::backward_algorithm_rep(
             const vector<vector<pair<double, double> > > &vals,
                                     const size_t start, const size_t end,
                                     const double lp_sf, const double lp_sb,
                                     const double lp_ff, const double lp_fb,
                                     const double lp_ft, const double lp_bf,
                                     const double lp_bb, const double lp_bt,
                                     const vector<unique_ptr<EmissionDistribution> > &fg_distro,
                                     const vector<unique_ptr<EmissionDistribution> > &bg_distro,
                                     vector<pair<double, double> > &b,
             const vector<bool> &array_status) const {

  size_t NREP = vals.size();
  b[end - 1].first = lp_ft;
  b[end - 1].second = lp_bt;

  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    double fg_a = b[k].first;
    double bg_a = b[k].second;

    for (size_t r = 0; r < NREP; ++r) {
      if (vals[r][k].first + vals[r][k].second >= 1) {
        fg_a += (*fg_distro[r])(vals[r][k]);
        bg_a += (*bg_distro[r])(vals[r][k]);
      }
    }
    b[i].first = log_sum_log(fg_a + lp_ff, bg_a + lp_fb);
    b[i].second = log_sum_log(fg_a + lp_bf, bg_a + lp_bb);
  }

  double emission_t1_b = 0;
  double emission_t1_f = 0;

  for (size_t r = 0; r < NREP; ++r) {
    if (vals[r][start].first + vals[r][start].second >= 1) {
      emission_t1_b += (*bg_distro[r])(vals[r][start]);
      emission_t1_f += (*fg_distro[r])(vals[r][start]);
    }
  }
  return log_sum_log(b[start].first + emission_t1_f + lp_sf,
                     b[start].second + emission_t1_b + lp_sb);
}


//ff_vals: ksi_t(1,1), where 1 is the S_1, i.e. posterior prob of transitions
void
TwoStateHMM::estimate_transitions_rep(
               const vector<vector<pair<double, double> > > &vals,
                                       const size_t start, const size_t end,
                                       const vector<pair<double, double> > &f,
                                       const vector<pair<double, double> > &b,
                                       const double total,
                                       const vector<unique_ptr<EmissionDistribution> > &fg_distro,
                                       const vector<unique_ptr<EmissionDistribution> > &bg_distro,
                                       const double lp_ff, const double lp_fb,
                                       const double lp_bf, const double lp_bb,
                                       const double lp_ft, const double lp_bt,
                                       vector<double> &ff_vals,
                                       vector<double> &fb_vals,
                                       vector<double> &bf_vals,
                                       vector<double> &bb_vals,
               const vector<bool> &array_status) const {
  size_t NREP = vals.size();
  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    double b_first = b[i].first - total;
    double b_second = b[i].second - total;

    for (size_t r = 0; r < NREP; ++r) {
      if(vals[r][i].first + vals[r][i].second >= 1) {
               b_first += (*fg_distro[r])(vals[r][i]);
               b_second += (*bg_distro[r])(vals[r][i]);
      }
    }
    const double ff = f[k].first;
    const double bb = f[k].second;

    ff_vals[k] = ff + lp_ff + b_first;
    fb_vals[k] = ff + lp_fb + b_second;

    bf_vals[k] = bb + lp_bf + b_first;
    bb_vals[k] = bb + lp_bb + b_second;
  }
}


double
TwoStateHMM::single_iteration_rep(const vector<vector<pair<double, double> > > &values,
                                   const vector<vector<double> > &vals_a_reps,
                                   const vector<vector<double> > &vals_b_reps,
                                   const vector<size_t> &reset_points,
                                   vector<pair<double, double> > &forward,
                                   vector<pair<double, double> > &backward,
                                   double &p_sf, double &p_sb,
                                   double &p_ff, double &p_fb, double &p_ft,
                                   double &p_bf, double &p_bb, double &p_bt,
                                   vector<unique_ptr<EmissionDistribution> > &fg_distro,
                                   vector<unique_ptr<EmissionDistribution> > &bg_distro,
           const vector<bool> &array_status) const {

  size_t NREP= values.size();
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
  vector<double> ff_vals(values[0].size(), 0);
  vector<double> fb_vals(values[0].size(), 0);
  vector<double> bf_vals(values[0].size(), 0);
  vector<double> bb_vals(values[0].size(), 0);

  // #pragma omp parallel
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score =
      forward_algorithm_rep(values,
                            reset_points[i],
                            reset_points[i + 1],
                            lp_sf, lp_sb,
                            lp_ff, lp_fb, lp_ft,
                            lp_bf, lp_bb, lp_bt,
                            fg_distro, bg_distro,
                            forward, array_status);

    const double backward_score =
      backward_algorithm_rep(values,
                             reset_points[i],
                             reset_points[i + 1],
                             lp_sf, lp_sb,
                             lp_ff, lp_fb, lp_ft,
                             lp_bf, lp_bb, lp_bt,
                             fg_distro, bg_distro,
                             backward, array_status);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10) {
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;
    }

    estimate_transitions_rep(values, reset_points[i], reset_points[i + 1],
                             forward, backward, score, fg_distro, bg_distro,
                             lp_ff, lp_fb, lp_bf, lp_bb, lp_ft, lp_bt,
                             ff_vals, fb_vals, bf_vals, bb_vals, array_status);

    total_score += score;
  }

  // Subtracting 1 from the limit of the summation
  // to eliminate the last term in the last block
  // And subtracting  (#blocks-1) from the sum
  // because the final term in each block has no
  // meaning since there is no transition to be counted
  // from the final observation (they all must go to terminal state)
  size_t NBLOCKS= reset_points.size()-1; //should equal to #deserts+1
  const double p_ff_new_estimate =
    exp(log_sum_log_vec(ff_vals, values[0].size() - 1)) - (NBLOCKS - 1);
  const double p_fb_new_estimate =
    exp(log_sum_log_vec(fb_vals, values[0].size() - 1)) - (NBLOCKS - 1);
  const double p_bf_new_estimate =
    exp(log_sum_log_vec(bf_vals, values[0].size() - 1)) - (NBLOCKS - 1);
  const double p_bb_new_estimate =
    exp(log_sum_log_vec(bb_vals, values[0].size() - 1)) - (NBLOCKS - 1);

  double denom = p_ff_new_estimate + p_fb_new_estimate;
  p_ff = (p_ff_new_estimate)/denom - p_ft/2.0;
  p_fb = (p_fb_new_estimate)/denom - p_ft/2.0;

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

  denom = p_bf_new_estimate + p_bb_new_estimate;
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

  p_sb = (p_bb + p_fb)/2.0;
  p_sf = (p_bf + p_ff)/2.0;

  vector<double> fg_probs(values[0].size(), 0);
  vector<double> bg_probs(values[0].size(), 0);
  estimate_emissions(forward, backward, fg_probs, bg_probs);

  vector<double> vals_a, vals_b;
  vector<double> fg_prob, bg_prob;
  for (size_t r = 0; r < NREP; ++r) {
    //individual replicate may have 0 coverage at some sites
    //remove these sites before fitting
    vals_a.clear();
    vals_b.clear();
    fg_prob.clear();
    bg_prob.clear();
    for (size_t i = 0; i < values[0].size(); ++i) {
      if (values[r][i].first + values[r][i].second >= 1) {
               vals_a.push_back(vals_a_reps[r][i]);
               vals_b.push_back(vals_b_reps[r][i]);
               fg_prob.push_back(fg_probs[i]); // use the common posterior prob
               bg_prob.push_back(bg_probs[i]);
      }
    }
    fg_distro[r]->fit(vals_a, vals_b, fg_prob);
    bg_distro[r]->fit(vals_a, vals_b, bg_prob);
  }
  return total_score;
}

double
TwoStateHMM::BaumWelchTraining_rep(
                  const vector<vector<pair<double, double> > > &values,
                  const vector<size_t> &reset_points,
                  vector<double> &start_trans,
                  vector<vector<double> > &trans,
                  vector<double> &end_trans,
                  vector<double> &fg_alpha, vector<double> &fg_beta,
                  vector<double> &bg_alpha, vector<double> &bg_beta,
                  const vector<bool> &array_status) const {

  vector<unique_ptr<EmissionDistribution> > fg_distro;
  vector<unique_ptr<EmissionDistribution> > bg_distro;

  size_t NREP = values.size();
  for (size_t i = 0; i < NREP; ++i) {
    if(array_status[i]) {
      fg_distro.emplace_back(new Beta(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new Beta(bg_alpha[i], bg_beta[i]));
    }
    else {
      fg_distro.emplace_back(new BetaBinomial(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new BetaBinomial(bg_alpha[i], bg_beta[i]));
    }
  }

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  const double score = BaumWelchTraining_rep(values, reset_points,
                                             start_trans[0], start_trans[1],
                                             trans[0][0], trans[0][1],
                                             end_trans[0], trans[1][0],
                                             trans[1][1], end_trans[1],
                                             fg_distro, bg_distro,
                                             array_status);
  for (size_t r = 0; r < NREP; ++r) {
    fg_alpha[r] = fg_distro[r]->getalpha();
    fg_beta[r] = fg_distro[r]->getbeta();
    bg_alpha[r] = bg_distro[r]->getalpha();
    bg_beta[r] = bg_distro[r]->getbeta();
  }

  return score;
}


double
TwoStateHMM::BaumWelchTraining_rep(
                const vector<vector<pair<double, double> > > &values,
                const vector<size_t> &reset_points,
                double &p_sf, double &p_sb,
                double &p_ff, double &p_fb, double &p_ft,
                double &p_bf, double &p_bb, double &p_bt,
                vector<unique_ptr<EmissionDistribution> > &fg_distro,
                vector<unique_ptr<EmissionDistribution> > &bg_distro,
                const vector<bool> &array_status) const {

  size_t NREP = values.size();
  vector<pair<double, double> > forward(values[0].size(),
                                        pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values[0].size(),
                                         pair<double, double>(0, 0));

  if (VERBOSE) {
    cerr << "MAX_ITER=" << max_iterations << "\tTOLERANCE="
         << tolerance << endl;
    cerr << setw(5)  << "ITR"
         << setw(10) << "F size"
         << setw(10) << "B size"
         << setw(14) << "DELTA"
         << endl;
  }

  double prev_total = -std::numeric_limits<double>::max();

  vector<vector<double> > vals_a_reps(NREP,
                                      vector<double>(values[0].size(), 0));
  vector<vector<double> > vals_b_reps(NREP,
                                      vector<double>(values[0].size(), 0));
  for (size_t r = 0; r < NREP; ++r) {
    if(array_status[r]) {
      for (size_t i = 0; i < values[0].size(); ++i) {
        if(values[r][i].second > 0) {
          vals_a_reps[r][i] =
            log(std::min(std::max(values[r][i].first/values[r][i].second,
                1e-2), 1.0 - 1e-2));
          vals_b_reps[r][i] =
            log(1 - std::min(std::max(values[r][i].first/values[r][i].second,
                1e-2), 1.0 - 1e-2));
        }
      }
    }
    else {
      for (size_t i = 0; i < values[0].size(); ++i) {
        if (values[r][i].first + values[r][i].second >= 1) {
              vals_a_reps[r][i] =
                log(std::min(std::max(values[r][i].first/(values[r][i].first +
                                  values[r][i].second), 1e-2), 1.0 - 1e-2));
              vals_b_reps[r][i] =
                log(1 - std::min(std::max(values[r][i].first/(values[r][i].first +
                                      values[r][i].second), 1e-2), 1.0 - 1e-2));
        }
      }
    }
  }

  for (size_t i = 0; i < max_iterations; ++i) {

    double p_sf_est = p_sf;
    double p_sb_est = p_sb;
    double p_ff_est = p_ff;
    double p_fb_est = p_fb;
    double p_bf_est = p_bf;
    double p_bb_est = p_bb;
    double p_ft_est = p_ft;
    double p_bt_est = p_bt;

    double total = single_iteration_rep(values, vals_a_reps, vals_b_reps,
                                        reset_points,
                                        forward, backward,
                                        p_sf_est, p_sb_est,
                                        p_ff_est, p_fb_est, p_ft_est,
                                        p_bf_est, p_bb_est, p_bt_est,
                                        fg_distro, bg_distro,array_status);

    if (DEBUG) {
      cerr << "S_F\t" << p_sf_est << endl
           << "S_B\t" << p_sb_est << endl
           << "F_F\t" << p_ff_est << endl
           << "F_B\t" << p_fb_est << endl
           << "B_F\t" << p_bf_est << endl
           << "B_B\t" << p_bb_est << endl
           << "F_E\t" << p_ft_est << endl
           << "B_E\t" << p_bt_est << endl
           << endl;
      for (size_t r = 0; r < NREP; ++r)
              cerr << "Emission parameters for Rep" << r+1
                   << setw(14) << fg_distro[r]->getalpha() << setw(14)
             << fg_distro[r]->getbeta() << setw(14) << bg_distro[r]->getalpha()
             << setw(14) << bg_distro[r]->getbeta() << endl;
    }

    if (VERBOSE) {
      cerr << setw(5) << i + 1
           << setw(10) << 1/p_fb_est
           << setw(10) << 1/p_bf_est
           << setw(14) << total
           << setw(14) << prev_total
           << setw(14) << total - prev_total
           << setw(14) << (total - prev_total)/std::fabs(total)
           << endl;
    }
    if (std::abs(total - prev_total) < tolerance) {
      if (VERBOSE)
              cerr << "CONVERGED" << "\t" << std::abs(total - prev_total)
             << "\t" << tolerance  << endl << endl;
      break;
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
TwoStateHMM::PosteriorScores_rep(
          const vector<vector<pair<double, double> > > &values,
                                  const vector<size_t> &reset_points,
                                  const vector<double> &start_trans,
                                  const vector<vector<double> > &trans,
                                  const vector<double> &end_trans,
                                  const vector<double> fg_alpha, const vector<double> fg_beta,
                                  const vector<double> bg_alpha, const vector<double> bg_beta,
                                  const vector<bool> &classes,
                                  vector<double> &llr_scores,
          const vector<bool> &array_status) const {

  vector<unique_ptr<EmissionDistribution> > fg_distro;
  vector<unique_ptr<EmissionDistribution> > bg_distro;

  size_t NREP = values.size();
  for (size_t i = 0; i < NREP; ++i) {
    if(array_status[i]) {
      fg_distro.emplace_back(new Beta(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new Beta(bg_alpha[i], bg_beta[i]));
    }
    else {
      fg_distro.emplace_back(new BetaBinomial(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new BetaBinomial(bg_alpha[i], bg_beta[i]));
    }
  }

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return PosteriorScores_rep(values, reset_points,
                     start_trans[0], start_trans[1],
                     trans[0][0], trans[0][1], end_trans[0],
                     trans[1][0], trans[1][1], end_trans[1],
                     fg_distro, bg_distro, classes, llr_scores, array_status);
}


void
TwoStateHMM::PosteriorScores_rep(
          const vector<vector<pair<double, double> > > &values,
                                  const vector<size_t> &reset_points,
                                  double p_sf, double p_sb,
                                  double p_ff, double p_fb, double p_ft,
                                  double p_bf, double p_bb, double p_bt,
                                  const vector<unique_ptr<EmissionDistribution> > &fg_distro,
                                  const vector<unique_ptr<EmissionDistribution> > &bg_distro,
                                  const vector<bool> &classes,
                                  vector<double> &llr_scores,
          const vector<bool> &array_status) const {

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

  vector<pair<double, double> > forward(values[0].size(),
                                        pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values[0].size(),
                                         pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm_rep(values,
                                               reset_points[i],
                                               reset_points[i + 1],
                                               lp_sf, lp_sb,
                                               lp_ff, lp_fb, lp_ft,
                                               lp_bf, lp_bb, lp_bt,
                                               fg_distro, bg_distro,
                                               forward, array_status);

    const double backward_score =
      backward_algorithm_rep(values,
                             reset_points[i],
                             reset_points[i + 1],
                             lp_sf, lp_sb,
                             lp_ff, lp_fb, lp_ft,
                             lp_bf, lp_bb, lp_bt,
                             fg_distro, bg_distro,
                             backward,array_status);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  llr_scores.resize(values[0].size());
  for (size_t i = 0; i < values[0].size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    if (classes[i])
      llr_scores[i] = (fg_state - bg_state);
    else
      llr_scores[i] = (bg_state - fg_state);
  }
}


double
TwoStateHMM::PosteriorDecoding_rep(
            const vector<vector<pair<double, double> > > &values,
                                    const vector<size_t> &reset_points,
                                    const vector<double> &start_trans,
                                    const vector<vector<double> > &trans,
                                    const vector<double> &end_trans,
                                    const vector<double> fg_alpha, const vector<double> fg_beta,
                                    const vector<double> bg_alpha, const vector<double> bg_beta,
                                    vector<bool> &classes,
                                    vector<double> &llr_scores,
            const vector<bool> &array_status) const {

  vector<unique_ptr<EmissionDistribution> > fg_distro;
  vector<unique_ptr<EmissionDistribution> > bg_distro;

  size_t NREP = values.size();
  for (size_t i = 0; i < NREP; ++i) {
    if(array_status[i]) {
      fg_distro.emplace_back(new Beta(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new Beta(bg_alpha[i], bg_beta[i]));
    }
    else {
      fg_distro.emplace_back(new BetaBinomial(fg_alpha[i], fg_beta[i]));
      bg_distro.emplace_back(new BetaBinomial(bg_alpha[i], bg_beta[i]));
    }
  }

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);

  return PosteriorDecoding_rep(values, reset_points,
                               start_trans[0], start_trans[1],
                               trans[0][0], trans[0][1], end_trans[0],
                               trans[1][0], trans[1][1], end_trans[1],
                               fg_distro, bg_distro, classes, llr_scores,array_status);
}


void
TwoStateHMM::TransitionPosteriors_rep(
      const vector<vector<pair<double, double> > > &values,
      const vector<size_t> &reset_points,
      double p_sf, double p_sb,
      double p_ff, double p_fb, double p_ft,
      double p_bf, double p_bb, double p_bt,
      const vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
      const vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
      const size_t transition, const vector<bool> &array_status,
      vector<double> &scores) const {

  double total_score = 0;
  size_t NREP = values.size();
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

  vector<pair<double, double> > forward(values[0].size(),
                                        pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values[0].size(),
                                         pair<double, double>(0, 0));
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm_rep(values,
                       reset_points[i],
                       reset_points[i + 1],
                       lp_sf, lp_sb,
                       lp_ff, lp_fb, lp_ft,
                       lp_bf, lp_bb, lp_bt,
                       fg_distro, bg_distro, forward, array_status);
    const double backward_score =
      backward_algorithm_rep(values, reset_points[i], reset_points[i + 1],
             lp_sf, lp_sb, lp_ff, lp_fb, lp_ft,
             lp_bf, lp_bb, lp_bt,
             fg_distro, bg_distro, backward, array_status);

    if (DEBUG && (fabs(score - backward_score)/
          max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
       << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }
  scores.resize(values[0].size());
  size_t j = 0;
  for (size_t i = 0; i < values[0].size(); ++i) {
    if (i == reset_points[j]) {
      ++j;
      scores[i] = 0;
    }
    else {
      double fg_vals=0;
      double bg_vals=0;
      for(size_t r = 0; r < NREP; ++r) {
        fg_vals+=(*fg_distro[r])(values[r][i]);
        bg_vals+=(*bg_distro[r])(values[r][i]);
      }
      fg_vals /= NREP;
      bg_vals /= NREP;

      const double fg_to_fg_state = forward[i - 1].first + lp_ff +
                                    fg_vals + backward[i].first;
      const double fg_to_bg_state = forward[i - 1].first + lp_fb +
                                    bg_vals + backward[i].second;
      const double bg_to_fg_state = forward[i - 1].second + lp_bf +
                                    fg_vals + backward[i].first;
      const double bg_to_bg_state = forward[i - 1].second + lp_bb +
                                    bg_vals + backward[i].second;
      const double denom = log_sum_log(
                           log_sum_log(fg_to_fg_state, fg_to_bg_state),
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
TwoStateHMM::PosteriorDecoding_rep(
            const vector<vector<pair<double, double> > > &values,
                                    const vector<size_t> &reset_points,
                                    double p_sf, double p_sb,
                                    double p_ff, double p_fb, double p_ft,
                                    double p_bf, double p_bb, double p_bt,
                                    const vector<unique_ptr<EmissionDistribution> > &fg_distro,
                                    const vector<unique_ptr<EmissionDistribution> > &bg_distro,
                                    vector<bool> &classes,
                                    vector<double> &llr_scores,
            const vector<bool> &array_status) const {
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

  vector<pair<double, double> > forward(values[0].size(),
                                        pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values[0].size(),
                                         pair<double, double>(0, 0));
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm_rep(values,
                                               reset_points[i],
                                               reset_points[i + 1],
                                               lp_sf, lp_sb,
                                               lp_ff, lp_fb, lp_ft,
                                               lp_bf, lp_bb, lp_bt,
                                               fg_distro, bg_distro,
                                               forward,array_status);

    const double backward_score =
      backward_algorithm_rep(values,
                             reset_points[i],
                             reset_points[i + 1],
                             lp_sf, lp_sb,
                             lp_ff, lp_fb, lp_ft,
                             lp_bf, lp_bb, lp_bt,
                             fg_distro, bg_distro,
                             backward,array_status);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  classes.resize(values[0].size());

  llr_scores.resize(values[0].size());
  for (size_t i = 0; i < values[0].size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    classes[i] = static_cast<bool>(fg_state > bg_state);
    llr_scores[i] = exp(fg_state - log_sum_log(fg_state, bg_state));
  }

  return total_score;
}

/***********End of functions for multiple replicates**************/
////////////////////////////////////////////////////////////////////////////////
