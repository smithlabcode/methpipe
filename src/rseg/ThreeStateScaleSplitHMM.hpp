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

#ifndef THREE_STATE_SCALE_SPLIT_HMM_HPP
#define THREE_STATE_SCALE_SPLIT_HMM_HPP

#include <memory>

#include "smithlab_utils.hpp"
#include "SplitDistro.hpp"

class ThreeStateScaleSplitHMM {
public:

  ThreeStateScaleSplitHMM(const double mp, const double tol,
		     const size_t max_itr, const bool v, bool d = false) :
    MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),
    VERBOSE(v), DEBUG(d) {}

  double
  ViterbiDecoding(const std::vector<double> &values,
                  const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans,
		  const std::vector<std::vector<double> > &trans,
		  const std::vector<double> &end_trans,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
		  std::vector<size_t> &ml_classes) const;

  double
  BaumWelchTraining(const std::vector<double> &values,
		    const std::vector<double> &vals_a,
		    const std::vector<double> &vals_b,
            const std::vector<double> &scales,
		    const std::vector<size_t> &reset_points,
		    std::vector<double> &start_trans,
		    std::vector<std::vector<double> > &trans,
		    std::vector<double> &end_trans,
		    SplitDistro &fg_distro,
		    SplitDistro &mid_distro,
		    SplitDistro &bg_distro) const;

  double
  PosteriorDecoding(const std::vector<double> &values,
                    const std::vector<double> &scales,
		    const std::vector<size_t> &reset_points,
		    const std::vector<double> &start_trans,
		    const std::vector<std::vector<double> > &trans,
		    const std::vector<double> &end_trans,
		    const SplitDistro &fg_distro,
		    const SplitDistro &mid_distro,
		    const SplitDistro &bg_distro,
		    std::vector<size_t> &classes,
		    std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<double> &values,
                  const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans,
		  const std::vector<std::vector<double> > &trans,
		  const std::vector<double> &end_trans,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
		  const std::vector<size_t> &classes,
		  std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<double> &values,
                  const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans,
		  const std::vector<std::vector<double> > &trans,
		  const std::vector<double> &end_trans,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
		  const size_t class_id,
		  std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<double> &values,
          const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans,
		  const std::vector<std::vector<double> > &trans,
		  const std::vector<double> &end_trans,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
          std::vector<double> &fg_scores,
          std::vector<double> &bg_scores) const;

  void
  TransitionPosteriors(const std::vector<double> &values,
                       const std::vector<double> &scales,
		       const std::vector<size_t> &reset_points,
		       const std::vector<double> &start_trans,
		       const std::vector<std::vector<double> > &trans,
		       const std::vector<double> &end_trans,
		       const SplitDistro &fg_distro,
		       const SplitDistro &mid_distro,
		       const SplitDistro &bg_distro,
		       std::vector<std::vector<double> > &scores) const;

  void
  TransitionPosteriors(const std::vector<double> &values,
                       const std::vector<double> &scales,
		       const std::vector<size_t> &reset_points,
		       const std::vector<double> &start_trans,
		       const std::vector<std::vector<double> > &trans,
		       const std::vector<double> &end_trans,
		       const SplitDistro &fg_distro,
		       const SplitDistro &mid_distro,
		       const SplitDistro &bg_distro,
               std::vector<std::vector<std::vector<double> > > &scores) const;

  std::string
  error_log() const;

private:

  double
  ViterbiDecoding(const std::vector<double> &values,
                  const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const double p_sf, const double p_sm, const double p_sb,
		  const double p_ff, const double p_fm, const double p_fb,
		  const double p_mf, const double p_mm, const double p_mb,
		  const double p_bf, const double p_bm, const double p_bb,
		  const double p_ft, const double p_bt, const double p_mt,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
		  std::vector<size_t> &ml_classes) const;

  double
  BaumWelchTraining(const std::vector<double> &values,
		    const std::vector<double> &vals_a,
		    const std::vector<double> &vals_b,
            const std::vector<double> &scales,
		    const std::vector<size_t> &reset_points,
		    double &p_sf, double &p_sm, double &p_sb,
		    double &p_ff, double &p_fm, double &p_fb,
		    double &p_mf, double &p_mm, double &p_mb,
		    double &p_bf, double &p_bm, double &p_bb,
		    double &p_ft, double &p_bt, double &p_mt,
		    SplitDistro &fg_distro,
		    SplitDistro &mid_distro,
		    SplitDistro &bg_distro) const;

  double
  PosteriorDecoding(const std::vector<double> &values,
                    const std::vector<double> &scales,
		    const std::vector<size_t> &reset_points,
		    const double p_sf, const double p_sm, const double p_sb,
		    const double p_ff, const double p_fm, const double p_fb,
		    const double p_mf, const double p_mm, const double p_mb,
		    const double p_bf, const double p_bm, const double p_bb,
		    const double p_ft, const double p_bt, const double p_mt,
		    const SplitDistro &fg_distro,
		    const SplitDistro &mid_distro,
		    const SplitDistro &bg_distro,
		    std::vector<size_t> &classes,
		    std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<double> &values,
                  const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const double p_sf, const double p_sm, const double p_sb,
		  const double p_ff, const double p_fm, const double p_fb,
		  const double p_mf, const double p_mm, const double p_mb,
		  const double p_bf, const double p_bm, const double p_bb,
		  const double p_ft, const double p_bt, const double p_mt,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
		  const std::vector<size_t> &classes,
		  std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<double> &values,
                  const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const double p_sf, const double p_sm, const double p_sb,
		  const double p_ff, const double p_fm, const double p_fb,
		  const double p_mf, const double p_mm, const double p_mb,
		  const double p_bf, const double p_bm, const double p_bb,
		  const double p_ft, const double p_bt, const double p_mt,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
		  const size_t class_id,
		  std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<double> &values,
          const std::vector<double> &scales,
		  const std::vector<size_t> &reset_points,
		  const double p_sf, const double p_sm, const double p_sb,
		  const double p_ff, const double p_fm, const double p_fb,
		  const double p_mf, const double p_mm, const double p_mb,
		  const double p_bf, const double p_bm, const double p_bb,
		  const double p_ft, const double p_bt, const double p_mt,
		  const SplitDistro &fg_distro,
		  const SplitDistro &mid_distro,
		  const SplitDistro &bg_distro,
          std::vector<double> &fg_scores,
		  std::vector<double> &bg_scores) const;

  void
  TransitionPosteriors(const std::vector<double> &values,
                       const std::vector<double> &scales,
		       const std::vector<size_t> &reset_points,
		       const double p_sf, const double p_sm, const double p_sb,
		       const double p_ff, const double p_fm, const double p_fb,
		       const double p_mf, const double p_mm, const double p_mb,
		       const double p_bf, const double p_bm, const double p_bb,
		       const double p_ft, const double p_bt, const double p_mt,
		       const SplitDistro &fg_distro,
		       const SplitDistro &mid_distro,
		       const SplitDistro &bg_distro,
		       std::vector<std::vector<double> > &scores) const;

  void
  TransitionPosteriors(const std::vector<double> &values,
                       const std::vector<double> &scales,
		       const std::vector<size_t> &reset_points,
		       const double p_sf, const double p_sm, const double p_sb,
		       const double p_ff, const double p_fm, const double p_fb,
		       const double p_mf, const double p_mm, const double p_mb,
		       const double p_bf, const double p_bm, const double p_bb,
		       const double p_ft, const double p_bt, const double p_mt,
		       const SplitDistro &fg_distro,
		       const SplitDistro &mid_distro,
		       const SplitDistro &bg_distro,
               std::vector<std::vector<std::vector<double> > > &scores) const;


  double
  single_iteration(const std::vector<double> &values,
		   const std::vector<double> &vals_a,
		   const std::vector<double> &vals_b,
           const std::vector<double> &scales,
           const std::vector<size_t> &reset_points,
		   std::vector<std::vector<double> > &forward,
		   std::vector<std::vector<double> > &backward,
		   double &p_sf, double &p_sm, double &p_sb,
		   double &p_ff, double &p_fm, double &p_fb,
		   double &p_mf, double &p_mm, double &p_mb,
		   double &p_bf, double &p_bm, double &p_bb,
		   double &p_ft, double &p_bt, double &p_mt,
		   SplitDistro &fg_distro,
		   SplitDistro &mid_distro,
		   SplitDistro &bg_distro) const;

  double
  forward_algorithm(const std::vector<double> &vals,
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
		    std::vector<std::vector<double> > &f) const;
  double
  backward_algorithm(const std::vector<double> &vals,
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
		     std::vector<std::vector<double> > &b) const;

  void
  estimate_emissions(const std::vector<std::vector<double> > &f,
		     const std::vector<std::vector<double> > &b,
		     std::vector<double> &fg_probs,
		     std::vector<double> &mid_probs,
		     std::vector<double> &bg_probs) const;

  void
  estimate_transitions(const std::vector<double> &vals,
            const std::vector<double> &scales,
		       const size_t start, const size_t end,
		       const std::vector<std::vector<double> > &f,
		       const std::vector<std::vector<double> > &b,
		       const double total,

		       const double lp_sf, const double lp_sm, const double lp_sb,
		       const double lp_ff, const double lp_fm, const double lp_fb,
		       const double lp_mf, const double lp_mm, const double lp_mb,
		       const double lp_bf, const double lp_bm, const double lp_bb,
		       const double lp_ft, const double lp_bt, const double lp_mt,

		       const SplitDistro &fg_distro,
		       const SplitDistro &mid_distro,
		       const SplitDistro &bg_distro,

		       std::vector<double> &ff_vals,
		       std::vector<double> &fm_vals,
		       std::vector<double> &fb_vals,

		       std::vector<double> &mf_vals,
		       std::vector<double> &mm_vals,
		       std::vector<double> &mb_vals,

		       std::vector<double> &bf_vals,
		       std::vector<double> &bm_vals,
		       std::vector<double> &bb_vals) const;

  double MIN_PROB;
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;

  mutable size_t emission_correction_count;
};

#endif
