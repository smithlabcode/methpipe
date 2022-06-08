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

#include "ThreeStateScaleSplitHMM.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>

#include "numerical_utils.hpp"

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


double
ThreeStateScaleSplitHMM::ViterbiDecoding(const vector<double> &values,
                                         const std::vector<double> &scales,
                                    const vector<size_t> &reset_points,
                                    const vector<double> &start_trans,
                                    const vector<vector<double> > &trans,
                                    const vector<double> &end_trans,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    vector<size_t> &ml_classes) const{

    return ViterbiDecoding(values, scales, reset_points,
                         start_trans[0], start_trans[1], start_trans[2],
                         trans[0][0], trans[0][1], trans[0][2],
                         trans[1][0], trans[1][1], trans[1][2],
                         trans[2][0], trans[2][1], trans[2][2],
                         end_trans[0], end_trans[1], end_trans[2],
                         fg_distro, mid_distro, bg_distro, ml_classes);
}


double
ThreeStateScaleSplitHMM::BaumWelchTraining(const vector<double> &values,
                                      const vector<double> &vals_a,
                                      const vector<double> &vals_b,
                      const std::vector<double> &scales,
                                      const vector<size_t> &reset_points,
                                      vector<double> &start_trans,
                                      vector<vector<double> > &trans,
                                      vector<double> &end_trans,
                                      SplitDistro &fg_distro,
                                      SplitDistro &mid_distro,
                                      SplitDistro &bg_distro) const {


    return BaumWelchTraining(values, vals_a, vals_b, scales, reset_points,
                           start_trans[0], start_trans[1], start_trans[2],
                           trans[0][0], trans[0][1], trans[0][2],
                           trans[1][0], trans[1][1], trans[1][2],
                           trans[2][0], trans[2][1], trans[2][2],
                           end_trans[0], end_trans[1], end_trans[2],
                           fg_distro, mid_distro, bg_distro);
}



double
ThreeStateScaleSplitHMM::PosteriorDecoding(const vector<double> &values,
                                           const std::vector<double> &scales,
                                      const vector<size_t> &reset_points,
                                      const vector<double> &start_trans,
                                      const vector<vector<double> > &trans,
                                      const vector<double> &end_trans,
                                      const SplitDistro &fg_distro,
                                      const SplitDistro &mid_distro,
                                      const SplitDistro &bg_distro,
                                      vector<size_t> &classes,
                                      vector<double> &llr_scores) const {

    return PosteriorDecoding(values, scales, reset_points,
                           start_trans[0], start_trans[1], start_trans[2],
                           trans[0][0], trans[0][1], trans[0][2],
                           trans[1][0], trans[1][1], trans[1][2],
                           trans[2][0], trans[2][1], trans[2][2],
                           end_trans[0], end_trans[1], end_trans[2],
                           fg_distro, mid_distro, bg_distro,
                           classes, llr_scores);

}

void
ThreeStateScaleSplitHMM::PosteriorScores(const vector<double> &values,
                                         const std::vector<double> &scales,
                                    const vector<size_t> &reset_points,
                                    const vector<double> &start_trans,
                                    const vector<vector<double> > &trans,
                                    const vector<double> &end_trans,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    const vector<size_t> &classes,
                                    vector<double> &llr_scores) const {

    return PosteriorScores(values, scales, reset_points,
                         start_trans[0], start_trans[1], start_trans[2],
                         trans[0][0], trans[0][1], trans[0][2],
                         trans[1][0], trans[1][1], trans[1][2],
                         trans[2][0], trans[2][1], trans[2][2],
                         end_trans[0], end_trans[1], end_trans[2],
                         fg_distro, mid_distro, bg_distro,
                         classes, llr_scores);
}


void
ThreeStateScaleSplitHMM::PosteriorScores(const vector<double> &values,
                                         const std::vector<double> &scales,
                                    const vector<size_t> &reset_points,
                                    const vector<double> &start_trans,
                                    const vector<vector<double> > &trans,
                                    const vector<double> &end_trans,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    const size_t class_id,
                                    vector<double> &llr_scores) const {

    return PosteriorScores(values, scales, reset_points,
                         start_trans[0], start_trans[1], start_trans[2],
                         trans[0][0], trans[0][1], trans[0][2],
                         trans[1][0], trans[1][1], trans[1][2],
                         trans[2][0], trans[2][1], trans[2][2],
                         end_trans[0], end_trans[1], end_trans[2],
                         fg_distro, mid_distro, bg_distro,
                         class_id, llr_scores);

}

void
ThreeStateScaleSplitHMM::PosteriorScores(const vector<double> &values,
                                         const std::vector<double> &scales,
                                    const vector<size_t> &reset_points,
                                    const vector<double> &start_trans,
                                    const vector<vector<double> > &trans,
                                    const vector<double> &end_trans,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                    vector<double> &fg_scores,
                                    vector<double> &bg_scores) const {

    return PosteriorScores(values, scales, reset_points,
                         start_trans[0], start_trans[1], start_trans[2],
                         trans[0][0], trans[0][1], trans[0][2],
                         trans[1][0], trans[1][1], trans[1][2],
                         trans[2][0], trans[2][1], trans[2][2],
                         end_trans[0], end_trans[1], end_trans[2],
                         fg_distro, mid_distro, bg_distro,
             fg_scores, bg_scores);
}

void
ThreeStateScaleSplitHMM::TransitionPosteriors(const vector<double> &values,
                                              const std::vector<double> &scales,
                                         const vector<size_t> &reset_points,
                                         const vector<double> &start_trans,
                                         const vector<vector<double> > &trans,
                                         const vector<double> &end_trans,
                                         const SplitDistro &fg_distro,
                                         const SplitDistro &mid_distro,
                                         const SplitDistro &bg_distro,
                                         vector<vector<double> > &scores) const {


    return TransitionPosteriors(values, scales, reset_points,
                              start_trans[0], start_trans[1], start_trans[2],
                              trans[0][0], trans[0][1], trans[0][2],
                              trans[1][0], trans[1][1], trans[1][2],
                              trans[2][0], trans[2][1], trans[2][2],
                              end_trans[0], end_trans[1], end_trans[2],
                              fg_distro, mid_distro, bg_distro, scores);
}



static size_t NON_START_END_STATES_COMBINATIONS = 9;

double
ThreeStateScaleSplitHMM::forward_algorithm(const vector<double> &vals,
                                           const std::vector<double> &scales,
                                      const size_t start, const size_t end,
                                      const double lp_sf, const double lp_sm, const double lp_sb,
                                      const double lp_ff, const double lp_fm, const double lp_fb,
                                      const double lp_mf, const double lp_mm, const double lp_mb,
                                      const double lp_bf, const double lp_bm, const double lp_bb,
                                      const double lp_ft, const double lp_bt, const double lp_mt,
                                      const SplitDistro &fg_distro,
                                      const SplitDistro &mid_distro,
                                      const SplitDistro &bg_distro,
                                      vector<vector<double> > &f) const {

    f[start][0] = lp_sf + fg_distro.log_likelihood(vals[start], scales[start]);
    f[start][1] = lp_sm + mid_distro.log_likelihood(vals[start], scales[start]);
    f[start][2] = lp_sb + bg_distro.log_likelihood(vals[start], scales[start]);

    for (size_t i = start + 1; i < end; ++i)
    {
        const size_t k = i - 1;
        f[i][0] = log_sum_log(f[k][0] + lp_ff, f[k][1] + lp_mf, + f[k][2] + lp_bf)
            + fg_distro.log_likelihood(vals[i], scales[i]);

        f[i][1] = log_sum_log(f[k][0] + lp_fm, f[k][1] + lp_mm, + f[k][2] + lp_bm)
            + mid_distro.log_likelihood(vals[i], scales[i]);

        f[i][2] = log_sum_log(f[k][0] + lp_fb, f[k][1] + lp_mb, + f[k][2] + lp_bb)
            + bg_distro.log_likelihood(vals[i], scales[i]);
    }

    return log_sum_log(f[end - 1][0] + lp_ft, f[end - 1][1] + lp_mt,
                       f[end - 1][2] + lp_bt);
}

void
ThreeStateScaleSplitHMM::TransitionPosteriors(const vector<double> &values,
                                              const std::vector<double> &scales,
                                         const vector<size_t> &reset_points,
                                         const vector<double> &start_trans,
                                         const vector<vector<double> > &trans,
                                         const vector<double> &end_trans,
                                         const SplitDistro &fg_distro,
                                         const SplitDistro &mid_distro,
                                         const SplitDistro &bg_distro,
                     vector<vector<vector<double> > > &scores) const {


    return TransitionPosteriors(values, scales, reset_points,
                              start_trans[0], start_trans[1], start_trans[2],
                              trans[0][0], trans[0][1], trans[0][2],
                              trans[1][0], trans[1][1], trans[1][2],
                              trans[2][0], trans[2][1], trans[2][2],
                              end_trans[0], end_trans[1], end_trans[2],
                              fg_distro, mid_distro, bg_distro, scores);
}

void
ThreeStateScaleSplitHMM::TransitionPosteriors(const vector<double> &values,
                                         const vector<double> &scales,
                                         const vector<size_t> &reset_points,
                                         const double p_sf, const double p_sm, const double p_sb,
                                         const double p_ff, const double p_fm, const double p_fb,
                                         const double p_mf, const double p_mm, const double p_mb,
                                         const double p_bf, const double p_bm, const double p_bb,
                                         const double p_ft, const double p_mt, const double p_bt,
                                         const SplitDistro &fg_distro,
                                         const SplitDistro &mid_distro,
                                         const SplitDistro &bg_distro,
                     vector<vector<vector<double> > > &scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  const size_t NUM_OF_STATES = 3;

  scores.resize(NUM_OF_STATES,
                vector<vector<double> >(NUM_OF_STATES,
                                        vector<double>(values.size(), 0)));
  size_t j = 0;
  for (size_t i = 0; i < values.size(); ++i)
  {
      if (i == reset_points[j])
      {
          ++j;
      }
      else
      {
          const double fg_to_fg_state = forward[i - 1][0] + lp_ff +
              fg_distro.log_likelihood(values[i], scales[i]) + backward[i][0];
          const double fg_to_mid_state = forward[i - 1][0] + lp_fm +
              mid_distro.log_likelihood(values[i], scales[i]) + backward[i][1];
          const double fg_to_bg_state = forward[i - 1][0] + lp_fb +
              bg_distro.log_likelihood(values[i], scales[i]) + backward[i][2];

          const double mid_to_fg_state = forward[i - 1][1] + lp_mf +
              fg_distro.log_likelihood(values[i], scales[i]) + backward[i][0];
          const double mid_to_mid_state = forward[i - 1][1] + lp_mm +
              mid_distro.log_likelihood(values[i], scales[i]) + backward[i][1];
          const double mid_to_bg_state = forward[i - 1][1] + lp_mb +
              bg_distro.log_likelihood(values[i], scales[i]) + backward[i][2];

          const double bg_to_fg_state = forward[i - 1][2] + lp_bf +
              fg_distro.log_likelihood(values[i], scales[i]) + backward[i][0];
          const double bg_to_mid_state = forward[i - 1][2] + lp_bm +
              mid_distro.log_likelihood(values[i], scales[i]) + backward[i][1];
          const double bg_to_bg_state = forward[i - 1][2] + lp_bb +
              bg_distro.log_likelihood(values[i], scales[i]) + backward[i][2];

          const double denom =
              log_sum_log(log_sum_log(fg_to_fg_state, fg_to_mid_state, fg_to_bg_state),
                          log_sum_log(mid_to_fg_state, mid_to_mid_state, mid_to_bg_state),
                          log_sum_log(bg_to_fg_state, bg_to_mid_state, bg_to_bg_state));

          scores[0][0][i] = exp(fg_to_fg_state - denom);
          scores[0][1][i] = exp(fg_to_mid_state - denom);
          scores[0][2][i] = exp(fg_to_bg_state - denom);
          scores[1][0][i] = exp(mid_to_fg_state - denom);
          scores[1][1][i] = exp(mid_to_mid_state - denom);
          scores[1][2][i] = exp(mid_to_bg_state - denom);
          scores[2][0][i] = exp(bg_to_fg_state - denom);
          scores[2][1][i] = exp(bg_to_mid_state - denom);
          scores[2][2][i] = exp(bg_to_bg_state - denom);
    }
  }
}




double
ThreeStateScaleSplitHMM::backward_algorithm(const vector<double> &vals,
                                            const std::vector<double> &scales,
                                       const size_t start, const size_t end,
                                       const double lp_sf, const double lp_sm, const double lp_sb,
                                       const double lp_ff, const double lp_fm, const double lp_fb,
                                       const double lp_mf, const double lp_mm, const double lp_mb,
                                       const double lp_bf, const double lp_bm, const double lp_bb,
                                       const double lp_ft, const double lp_bt, const double lp_mt,
                                       const SplitDistro &fg_distro,
                                       const SplitDistro &mid_distro,
                                       const SplitDistro &bg_distro,
                                       vector<vector<double> > &b) const {
  b[end - 1][0] = lp_ft;
  b[end - 1][1] = lp_mt;
  b[end - 1][2] = lp_bt;

  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    const double fg_a = fg_distro.log_likelihood(vals[k], scales[k]) + b[k][0];
    const double mid_a = mid_distro.log_likelihood(vals[k], scales[k]) + b[k][1];
    const double bg_a = bg_distro.log_likelihood(vals[k], scales[k]) + b[k][2];

    b[i][0] = log_sum_log(fg_a + lp_ff,
                          mid_a + lp_fm,
                          bg_a + lp_fb);

    b[i][1] = log_sum_log(fg_a + lp_mf,
                          mid_a + lp_mm,
                          bg_a + lp_mb);

    b[i][2] = log_sum_log(fg_a + lp_bf,
                          mid_a + lp_bm,
                          bg_a + lp_bb);

  }
  return log_sum_log(b[start][0] + fg_distro.log_likelihood(vals[start], scales[start]) + lp_sf,
                     b[start][1] + mid_distro.log_likelihood(vals[start], scales[start]) + lp_sm,
                     b[start][2] + bg_distro.log_likelihood(vals[start], scales[start]) + lp_sb);
}

void
ThreeStateScaleSplitHMM::estimate_emissions(const vector<vector<double> > &f,
                                       const vector<vector<double> > &b,
                                       vector<double> &fg_probs,
                                       vector<double> &mid_probs,
                                       vector<double> &bg_probs) const {
  for (size_t i = 0; i < b.size(); ++i) {
    const double fg = (f[i][0] + b[i][0]);
    const double mid = (f[i][1] + b[i][1]);
    const double bg = (f[i][2] + b[i][2]);
    const double denom = log_sum_log(fg, mid, bg);
    fg_probs[i] = exp(fg - denom);
    mid_probs[i] = exp(mid - denom);
    bg_probs[i] = exp(bg - denom);
  }
}


void
ThreeStateScaleSplitHMM::estimate_transitions(const vector<double> &vals,
                                              const std::vector<double> &scales,
                                         const size_t start, const size_t end,
                                         const vector<vector<double> > &f,
                                         const vector<vector<double> > &b,
                                         const double total,

                                         const double lp_sf, const double lp_sm, const double lp_sb,
                                         const double lp_ff, const double lp_fm, const double lp_fb,
                                         const double lp_mf, const double lp_mm, const double lp_mb,
                                         const double lp_bf, const double lp_bm, const double lp_bb,
                                         const double lp_ft, const double lp_mt, const double lp_bt,

                                         const SplitDistro &fg_distro,
                                         const SplitDistro &mid_distro,
                                         const SplitDistro &bg_distro,

                                         vector<double> &ff_vals,
                                         vector<double> &fm_vals,
                                         vector<double> &fb_vals,

                                         vector<double> &mf_vals,
                                         vector<double> &mm_vals,
                                         vector<double> &mb_vals,

                                         vector<double> &bf_vals,
                                         vector<double> &bm_vals,
                                         vector<double> &bb_vals) const {

    for (size_t i = start + 1; i < end; ++i)
    {
        const size_t k = i - 1;

        const double lp_fg = fg_distro.log_likelihood(vals[i], scales[i]) - total;
        const double lp_mg = mid_distro.log_likelihood(vals[i], scales[i]) - total;
        const double lp_bg = bg_distro.log_likelihood(vals[i], scales[i]) - total;

        ff_vals[k] = f[k][0] + lp_ff + lp_fg + b[i][0];
        fm_vals[k] = f[k][0] + lp_fm + lp_mg + b[i][1];
        fb_vals[k] = f[k][0] + lp_fb + lp_bg + b[i][2];

        mf_vals[k] = f[k][1] + lp_mf + lp_fg + b[i][0];
        mm_vals[k] = f[k][1] + lp_mm + lp_mg + b[i][1];
        mb_vals[k] = f[k][1] + lp_mb + lp_bg + b[i][2];

        bf_vals[k] = f[k][2] + lp_bf + lp_fg + b[i][0];
        bm_vals[k] = f[k][2] + lp_bm + lp_mg + b[i][1];
        bb_vals[k] = f[k][2] + lp_bb + lp_bg + b[i][2];
    }
}



double
ThreeStateScaleSplitHMM::single_iteration(const vector<double> &values,
                                     const vector<double> &vals_a,
                                     const vector<double> &vals_b,
                     const std::vector<double> &scales,
                                     const vector<size_t> &reset_points,

                                     vector<vector<double> > &forward,
                                     vector<vector<double> > &backward,

                                     double &p_sf, double &p_sm, double &p_sb,
                                     double &p_ff, double &p_fm, double &p_fb,
                                     double &p_mf, double &p_mm, double &p_mb,
                                     double &p_bf, double &p_bm, double &p_bb,
                                     double &p_ft, double &p_mt, double &p_bt,

                                     SplitDistro &fg_distro,
                                     SplitDistro &mid_distro,
                                     SplitDistro &bg_distro) const {

  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  // for estimating transitions
  vector<double> ff_vals(values.size(), 0);
  vector<double> fm_vals(values.size(), 0);
  vector<double> fb_vals(values.size(), 0);

  vector<double> mf_vals(values.size(), 0);
  vector<double> mm_vals(values.size(), 0);
  vector<double> mb_vals(values.size(), 0);

  vector<double> bf_vals(values.size(), 0);
  vector<double> bm_vals(values.size(), 0);
  vector<double> bb_vals(values.size(), 0);

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);
    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

    if (DEBUG && fabs((score - backward_score)/
        max(score, backward_score)) > 1e-20)
    {
      cerr << "fabs(score - backward_score)/max(score, backward_score) > 1e-20" << endl;
      cerr << std::setprecision(12) << score << "\t" << std::setprecision(12) << backward_score << endl;

    }

    estimate_transitions(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         forward, backward,
                         score,
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,

                         fg_distro,
                         mid_distro,
                         bg_distro,

                         ff_vals,
                         fm_vals,
                         fb_vals,

                         mf_vals,
                         mm_vals,
                         mb_vals,

                         bf_vals,
                         bm_vals,
                         bb_vals);

    total_score += score;
  }

  // Subtracting 1 from the limit of the summation because the final
  // term has no meaning since there is no transition to be counted
  // from the final observation (they all must go to terminal state)
  //  SQ: ff_vals[reset_points[i] - 1] is always equal to 0, and should not be considered
  // i.e. the start of each region does not contribute to transition emission estimate
  const double p_ff_new_estimate = exp(log_sum_log_vec(ff_vals, values.size()))
      - reset_points.size() + 1;
  const double p_fm_new_estimate = exp(log_sum_log_vec(fm_vals, values.size()))
      - reset_points.size() + 1;
  const double p_fb_new_estimate = exp(log_sum_log_vec(fb_vals, values.size()))
      - reset_points.size() + 1;

  const double p_mf_new_estimate = exp(log_sum_log_vec(mf_vals, values.size()))
      - reset_points.size() + 1;
  const double p_mm_new_estimate = exp(log_sum_log_vec(mm_vals, values.size()))
      - reset_points.size() + 1;
  const double p_mb_new_estimate = exp(log_sum_log_vec(mb_vals, values.size()))
      - reset_points.size() + 1;

  const double p_bf_new_estimate = exp(log_sum_log_vec(bf_vals, values.size()))
      - reset_points.size() + 1;
  const double p_bm_new_estimate = exp(log_sum_log_vec(bm_vals, values.size()))
      - reset_points.size() + 1;
  const double p_bb_new_estimate = exp(log_sum_log_vec(bb_vals, values.size()))
      - reset_points.size() + 1;


  double denom = (p_ff_new_estimate + p_fm_new_estimate + p_fb_new_estimate);
  p_ff = p_ff_new_estimate/denom;
  p_fm = p_fm_new_estimate / denom;
  p_fb = p_fb_new_estimate / denom;


  if (p_ff < MIN_PROB) {
    if (DEBUG)
      cerr << "p_ff < MIN_PROB" << endl;
    p_ff = MIN_PROB;
  }
  if (p_fm < MIN_PROB) {
    if (DEBUG)
      cerr << "p_fm < MIN_PROB" << endl;
    p_fm = MIN_PROB;
  }
  if (p_fb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_fb < MIN_PROB" << endl;
    p_fb = MIN_PROB;
  }

  denom = (p_mf_new_estimate + p_mm_new_estimate + p_mb_new_estimate);
  p_mm = p_mm_new_estimate/denom;
  p_mf = p_mf_new_estimate / denom;
  p_mb = p_mb_new_estimate / denom;

  if (p_mf < MIN_PROB) {
    if (DEBUG)
      cerr << "p_mf < MIN_PROB" << endl;
    p_mf = MIN_PROB;
  }
  if (p_mm < MIN_PROB) {
    if (DEBUG)
      cerr << "p_mm < MIN_PROB" << endl;
    p_mm = MIN_PROB;
  }
  if (p_mb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_mb < MIN_PROB" << endl;
    p_mb = MIN_PROB;
  }


  denom = (p_bf_new_estimate + p_bm_new_estimate + p_bb_new_estimate);
  p_bb = p_bb_new_estimate/denom;
  p_bf = p_bf_new_estimate / denom;
  p_bm = p_bm_new_estimate / denom;


  if (p_bf < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bf < MIN_PROB" << endl;
    p_bf = MIN_PROB;
  }
  if (p_bm < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bm < MIN_PROB" << endl;
    p_bm = MIN_PROB;
  }
  if (p_bb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bb < MIN_PROB" << endl;
    p_bb = MIN_PROB;
  }

  // for estimating emissions
  vector<double> fg_probs(values.size());
  vector<double> mid_probs(values.size());
  vector<double> bg_probs(values.size());
  estimate_emissions(forward, backward, fg_probs, mid_probs, bg_probs);

  fg_distro.estimate_params_ml(vals_a, vals_b, scales, fg_probs);
  mid_distro.estimate_params_ml(vals_a, vals_b, scales, mid_probs);
  bg_distro.estimate_params_ml(vals_a, vals_b, scales, bg_probs);

  // make sure the mean satisfies fg < mid < bg
  if (fg_distro.get_mean() < mid_distro.get_mean())
  {
      SplitDistro tmp_distro(fg_distro);
      fg_distro = mid_distro;
      mid_distro = tmp_distro;

      std::swap(p_sf, p_sm);
      std::swap(p_ff, p_mm);
      std::swap(p_fm, p_mf);
      std::swap(p_fb, p_mb);
      std::swap(p_bf, p_bm);
      std::swap(p_ft, p_mt);
  }

  if (fg_distro.get_mean() < bg_distro.get_mean())
  {
      SplitDistro tmp_distro(fg_distro);
      fg_distro = bg_distro;
      bg_distro = tmp_distro;

      std::swap(p_sf, p_sb);
      std::swap(p_ff, p_bb);
      std::swap(p_fm, p_bm);
      std::swap(p_fb, p_bf);
      std::swap(p_mf, p_mb);
      std::swap(p_ft, p_bt);
  }

  if (mid_distro.get_mean() < bg_distro.get_mean())
  {
      SplitDistro tmp_distro(mid_distro);
      mid_distro = bg_distro;
      bg_distro = tmp_distro;

      std::swap(p_sm, p_sb);
      std::swap(p_mf, p_bf);
      std::swap(p_mm, p_bb);
      std::swap(p_mb, p_bm);
      std::swap(p_fm, p_fb);
      std::swap(p_mt, p_bt);
  }

  return total_score;
}



double
ThreeStateScaleSplitHMM::BaumWelchTraining(const vector<double> &values,
                                      const vector<double> &vals_a,
                                      const vector<double> &vals_b,
                      const std::vector<double> &scales,
                                      const vector<size_t> &reset_points,
                                      double &p_sf, double &p_sm, double &p_sb,
                                      double &p_ff, double &p_fm, double &p_fb,
                                      double &p_mf, double &p_mm, double &p_mb,
                                      double &p_bf, double &p_bm, double &p_bb,
                                      double &p_ft, double &p_mt, double &p_bt,
                                      SplitDistro &fg_distro,
                                      SplitDistro &mid_distro,
                                      SplitDistro &bg_distro) const {

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  if (VERBOSE)
    cout << setw(4) << "ITR"
         << setw(8) << "F BINS"
         << setw(8) << "M BINS"
         << setw(8) << "B BINS"
         << setw(32) << "F PARAMS"
         << setw(32) << "M PARAMS"
         << setw(32) << "B PARAMS"
         << setw(14) << "DELTA"
         << endl;

  double prev_total = -std::numeric_limits<double>::max();

  for (size_t i = 0; i < max_iterations; ++i) {


      double total = single_iteration(values, vals_a, vals_b, scales, reset_points,
                                    forward, backward,

                                    p_sf, p_sm, p_sb,
                                    p_ff, p_fm, p_fb,
                                    p_mf, p_mm, p_mb,
                                    p_bf, p_bm, p_bb,
                                    p_ft, p_mt, p_bt,

                                    fg_distro, mid_distro, bg_distro);

    if (std::fabs((prev_total - total)/prev_total) < tolerance) {
      if (VERBOSE)
        cout << "CONVERGED" << endl << endl;
      break;
    }

    if (total < prev_total) {
        if (VERBOSE)
            cout << "Oscillaion point." << endl;
        break;
    }

    if (VERBOSE) {
      cout << setw(4) << i + 1
           << setw(8) << 1.0/(p_fb + p_fm)
           << setw(8) << 1.0/(p_mf + p_mb)
           << setw(8) << 1.0/(p_bf + p_bm)
           << setw(32) << fg_distro
           << setw(32) << mid_distro
           << setw(32) << bg_distro
           << setw(14) << (prev_total - total)/prev_total
           << "\t" << total
           << endl;
    }

    prev_total = total;
  }
  return prev_total;
}


void
ThreeStateScaleSplitHMM::PosteriorScores(const vector<double> &values,
                                const std::vector<double> &scales,
                    const vector<size_t> &reset_points,
                                    const double p_sf, const double p_sm, const double p_sb,
                                    const double p_ff, const double p_fm, const double p_fb,
                                    const double p_mf, const double p_mm, const double p_mb,
                                    const double p_bf, const double p_bm, const double p_bb,
                                    const double p_ft, const double p_mt, const double p_bt,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    const vector<size_t> &classes,
                                    vector<double> &llr_scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i][0] + backward[i][0];
    const double mid_state = forward[i][1] + backward[i][1];
    const double bg_state = forward[i][2] + backward[i][2];
    const double denom = log_sum_log(fg_state, mid_state, bg_state);

    if (classes[i] == 0)
        llr_scores[i] = exp(fg_state - denom);
    else if (classes[i] == 1)
        llr_scores[i] = exp(mid_state - denom);
    else
        llr_scores[i] = exp(bg_state - denom);
  }
}


void
ThreeStateScaleSplitHMM::PosteriorScores(const vector<double> &values,
                                         const std::vector<double> &scales,
                                    const vector<size_t> &reset_points,
                                    const double p_sf, const double p_sm, const double p_sb,
                                    const double p_ff, const double p_fm, const double p_fb,
                                    const double p_mf, const double p_mm, const double p_mb,
                                    const double p_bf, const double p_bm, const double p_bb,
                                    const double p_ft, const double p_mt, const double p_bt,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    const size_t class_id,
                                    vector<double> &llr_scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const double score = forward_algorithm(values,  scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i][0] + backward[i][0];
    const double mid_state = forward[i][1] + backward[i][1];
    const double bg_state = forward[i][2] + backward[i][2];
    const double denom = log_sum_log(fg_state, mid_state, bg_state);

    if (class_id == 0)
        llr_scores[i] = exp(fg_state - denom);
    else if (class_id == 1)
        llr_scores[i] = exp(mid_state - denom);
    else
        llr_scores[i] = exp(bg_state - denom);
  }
}


void
ThreeStateScaleSplitHMM::PosteriorScores(const vector<double> &values,
                            const std::vector<double> &scales,
                    const vector<size_t> &reset_points,
                                    const double p_sf, const double p_sm, const double p_sb,
                                    const double p_ff, const double p_fm, const double p_fb,
                                    const double p_mf, const double p_mm, const double p_mb,
                                    const double p_bf, const double p_bm, const double p_bb,
                                    const double p_ft, const double p_mt, const double p_bt,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    vector<double> &fg_scores,
                                    vector<double> &bg_scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  fg_scores.resize(values.size());
  bg_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i][0] + backward[i][0];
    const double mid_state = forward[i][1] + backward[i][1];
    const double bg_state = forward[i][2] + backward[i][2];
    const double denom = log_sum_log(fg_state, mid_state, bg_state);

    fg_scores[i] = exp(fg_state - denom);
    bg_scores[i] = exp(bg_state - denom);
  }
}


void
ThreeStateScaleSplitHMM::TransitionPosteriors(const vector<double> &values,
                                              const std::vector<double> &scales,
                                         const vector<size_t> &reset_points,
                                         const double p_sf, const double p_sm, const double p_sb,
                                         const double p_ff, const double p_fm, const double p_fb,
                                         const double p_mf, const double p_mm, const double p_mb,
                                         const double p_bf, const double p_bm, const double p_bb,
                                         const double p_ft, const double p_mt, const double p_bt,
                                         const SplitDistro &fg_distro,
                                         const SplitDistro &mid_distro,
                                         const SplitDistro &bg_distro,
                                         vector<vector<double> > &scores) const {


  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

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
      scores[i] = vector<double>(NON_START_END_STATES_COMBINATIONS, 0);
    }
    else {
      const double fg_to_fg_state = forward[i - 1][0] + lp_ff +
          fg_distro.log_likelihood(values[i], scales[i]) + backward[i][0];
      const double fg_to_mid_state = forward[i - 1][0] + lp_fm +
        mid_distro.log_likelihood(values[i], scales[i]) + backward[i][1];
      const double fg_to_bg_state = forward[i - 1][0] + lp_fb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i][2];

      const double mid_to_fg_state = forward[i - 1][1] + lp_mf +
        fg_distro.log_likelihood(values[i], scales[i]) + backward[i][0];
      const double mid_to_mid_state = forward[i - 1][1] + lp_mm +
        mid_distro.log_likelihood(values[i], scales[i]) + backward[i][1];
      const double mid_to_bg_state = forward[i - 1][1] + lp_mb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i][2];

      const double bg_to_fg_state = forward[i - 1][2] + lp_bf +
        fg_distro.log_likelihood(values[i], scales[i]) + backward[i][0];
      const double bg_to_mid_state = forward[i - 1][2] + lp_bm +
        mid_distro.log_likelihood(values[i], scales[i]) + backward[i][1];
      const double bg_to_bg_state = forward[i - 1][2] + lp_bb +
        bg_distro.log_likelihood(values[i], scales[i]) + backward[i][2];

      const double denom =
        log_sum_log(log_sum_log(log_sum_log(log_sum_log(fg_to_fg_state, fg_to_mid_state),
                                            fg_to_bg_state),
                                log_sum_log(log_sum_log(mid_to_fg_state, mid_to_mid_state),
                                            mid_to_bg_state)),
                    log_sum_log(log_sum_log(bg_to_fg_state, bg_to_mid_state),
                                bg_to_bg_state));

      scores[i].resize(NON_START_END_STATES_COMBINATIONS);
      scores[i][0] = exp(fg_to_fg_state - denom);
      scores[i][1] = exp(fg_to_mid_state - denom);
      scores[i][2] = exp(fg_to_bg_state - denom);
      scores[i][3] = exp(mid_to_fg_state - denom);
      scores[i][4] = exp(mid_to_mid_state - denom);
      scores[i][5] = exp(mid_to_bg_state - denom);
      scores[i][6] = exp(bg_to_fg_state - denom);
      scores[i][7] = exp(bg_to_mid_state - denom);
      scores[i][8] = exp(bg_to_bg_state - denom);

    }
  }
}




double
ThreeStateScaleSplitHMM::PosteriorDecoding(const vector<double> &values,
                                           const std::vector<double> &scales,
                                      const vector<size_t> &reset_points,
                                      const double p_sf, const double p_sm, const double p_sb,
                                      const double p_ff, const double p_fm, const double p_fb,
                                      const double p_mf, const double p_mm, const double p_mb,
                                      const double p_bf, const double p_bm, const double p_bb,
                                      const double p_ft, const double p_mt, const double p_bt,
                                      const SplitDistro &fg_distro,
                                      const SplitDistro &mid_distro,
                                      const SplitDistro &bg_distro,
                                      vector<size_t> &classes,
                                      vector<double> &llr_scores) const {
  double total_score = 0;

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  assert(isfinite(lp_sf) &&
         isfinite(lp_sm) &&
         isfinite(lp_sb) &&
         isfinite(lp_ff) &&
         isfinite(lp_fm) &&
         isfinite(lp_fb) &&
         isfinite(lp_ft) &&
         isfinite(lp_mf) &&
         isfinite(lp_mm) &&
         isfinite(lp_mb) &&
         isfinite(lp_mt) &&
         isfinite(lp_bf) &&
         isfinite(lp_bm) &&
         isfinite(lp_bb) &&
         isfinite(lp_bt));

  vector<vector<double> > forward(values.size(), vector<double>(3, 0));
  vector<vector<double> > backward(values.size(), vector<double>(3, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const double score = forward_algorithm(values, scales,
                                           reset_points[i],
                                           reset_points[i + 1],
                                           lp_sf, lp_sm, lp_sb,
                                           lp_ff, lp_fm, lp_fb,
                                           lp_mf, lp_mm, lp_mb,
                                           lp_bf, lp_bm, lp_bb,
                                           lp_ft, lp_mt, lp_bt,
                                           fg_distro,
                                           mid_distro,
                                           bg_distro,
                                           forward);

    const double backward_score =
        backward_algorithm(values, scales,
                         reset_points[i],
                         reset_points[i + 1],
                         lp_sf, lp_sm, lp_sb,
                         lp_ff, lp_fm, lp_fb,
                         lp_mf, lp_mm, lp_mb,
                         lp_bf, lp_bm, lp_bb,
                         lp_ft, lp_mt, lp_bt,
                         fg_distro,
                         mid_distro,
                         bg_distro,
                         backward);

    if (DEBUG && (fabs(score - backward_score)/
                  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
           << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }

  classes.resize(values.size());
  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i][0] + backward[i][0];
    const double mid_state = forward[i][1] + backward[i][1];
    const double bg_state = forward[i][2] + backward[i][2];
    const double denom = log_sum_log(fg_state, mid_state, bg_state);

    if (fg_state > mid_state && fg_state > bg_state)
      classes[i] = 0;
    else if (mid_state > fg_state && mid_state > bg_state)
      classes[i] = 1;
    else
      classes[i] = 2;

    if (classes[i] == 0)
        llr_scores[i] = exp(fg_state - denom);
    else if (classes[i] == 1)
        llr_scores[i] = exp(mid_state - denom);
    else
        llr_scores[i] = exp(bg_state - denom);
  }

  return total_score;
}



/*************************************************************
 *
 * Functions for Viterbi training and decoding.
 *
 *************************************************************/

double
ThreeStateScaleSplitHMM::ViterbiDecoding(const vector<double> &values,
                                         const std::vector<double> &scales,
                                    const vector<size_t> &reset_points,
                                    const double p_sf, const double p_sm, const double p_sb,
                                    const double p_ff, const double p_fm, const double p_fb,
                                    const double p_mf, const double p_mm, const double p_mb,
                                    const double p_bf, const double p_bm, const double p_bb,
                                    const double p_ft, const double p_bt, const double p_mt,
                                    const SplitDistro &fg_distro,
                                    const SplitDistro &mid_distro,
                                    const SplitDistro &bg_distro,
                                    vector<size_t> &ml_classes) const {

  const double lp_sf = log(p_sf);
  const double lp_sm = log(p_sm);
  const double lp_sb = log(p_sb);

  const double lp_ff = log(p_ff);
  const double lp_fm = log(p_fm);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);

  const double lp_mf = log(p_mf);
  const double lp_mm = log(p_mm);
  const double lp_mb = log(p_mb);
  const double lp_mt = log(p_mt);

  const double lp_bf = log(p_bf);
  const double lp_bm = log(p_bm);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);

  // ml_classes = vector<bool>(values.size());
  double total = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

      const size_t start = reset_points[i];
    const size_t lim = reset_points[i + 1] - reset_points[i];

    vector<vector<double> > v(lim, vector<double>(3, 0));
    vector<vector<size_t> > trace(lim, vector<size_t>(3, 0));

    v.front()[0] = fg_distro.log_likelihood(values[start], scales[start]) + lp_sf;
    v.front()[1] = mid_distro.log_likelihood(values[start], scales[start]) + lp_sm;
    v.front()[2] = bg_distro.log_likelihood(values[start], scales[start]) + lp_sb;

    for (size_t j = 1; j < lim; ++j) {

      const double ff = v[j - 1][0] + lp_ff;
      const double mf = v[j - 1][1] + lp_mf;
      const double bf = v[j - 1][2] + lp_bf;

      const double fg_log_emmit =
          fg_distro.log_likelihood(values[start + j], scales[start + j]);
      if (ff >= bf && ff >= mf) {
        v[j][0] = fg_log_emmit + ff;
        trace[j][0] = 0;
      }
      else if (mf >= ff && mf >= bf) {
        v[j][0] = fg_log_emmit + mf;
        trace[j][0] = 1;
      }
      else {
        v[j][0] = fg_log_emmit + bf;
        trace[j][0] = 2;
      }

      const double fm = v[j - 1][0] + lp_fm;
      const double mm = v[j - 1][1] + lp_mm;
      const double bm = v[j - 1][2] + lp_bm;

      const double mid_log_emmit =
        mid_distro.log_likelihood(values[start + j], scales[start + j]);
      if (fm >= bm && fm >= mm) {
        v[j][1] = mid_log_emmit + fm;
        trace[j][1] = 0;
      }
      else if (mm >= fm && mm >= bm) {
        v[j][1] = mid_log_emmit + mm;
        trace[j][1] = 1;
      }
      else {
        v[j][1] = mid_log_emmit + bm;
        trace[j][1] = 2;
      }

      const double fb = v[j - 1][0] + lp_fb;
      const double mb = v[j - 1][1] + lp_mb;
      const double bb = v[j - 1][2] + lp_bb;

      const double bg_log_emmit =
        bg_distro.log_likelihood(values[start + j], scales[start + j]);
      if (fb >= bb && fb >= mb) {
        v[j][2] = bg_log_emmit + fb;
        trace[j][2] = 0;
      }
      else if (mb >= fb && mb >= bb) {
        v[j][2] = bg_log_emmit + mb;
        trace[j][2] = 1;
      }
      else {
        v[j][2] = bg_log_emmit + bb;
        trace[j][2] = 2;
      }
    }

    v.back()[0] += lp_ft;
    v.back()[1] += lp_mt;
    v.back()[2] += lp_bt;

    vector<size_t> inner_ml_classes;

    // do the traceback
    size_t prev = 0;
    if (v.back()[0] >= v.back()[1] && v.back()[0] >= v.back()[2]) {
      inner_ml_classes.push_back(0);
      prev = trace.back()[0];
    }
    else if (v.back()[1] >= v.back()[0] && v.back()[1] >= v.back()[2]) {
      inner_ml_classes.push_back(1);
      prev = trace.back()[1];
    }
    else {
      inner_ml_classes.push_back(2);
      prev = trace.back()[2];
    }

    for (size_t j = trace.size() - 1; j > 0; --j) {
      const size_t k = j - 1;
      if (prev == 0) {
        inner_ml_classes.push_back(0);
        prev = trace[k][0];
      }
      else if (prev == 1) {
        inner_ml_classes.push_back(1);
        prev = trace[k][1];
      }
      else {
        inner_ml_classes.push_back(2);
        prev = trace[k][2];
      }
    }

    reverse(inner_ml_classes.begin(), inner_ml_classes.end());
    ml_classes.insert(ml_classes.end(), inner_ml_classes.begin(),
                      inner_ml_classes.end());

    total += max(v.back()[0], max(v.back()[1], v.back()[2]));

  }

  return total;
}
