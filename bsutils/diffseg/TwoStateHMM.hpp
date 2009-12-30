/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Authors: Andrew D. Smith

  This file is part of rmap.

  rmap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  rmap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rmap; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef TWO_STATE_HMM_HPP
#define TWO_STATE_HMM_HPP

#include "rmap_utils.hpp"
#include <memory>
#include <vector>

using std::vector;

struct geometric;

class TwoStateHMMB {
public:
  
  TwoStateHMMB(const double mp, const double tol,
	       const size_t max_itr, const bool v, bool d = false) : 
    MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),    
    VERBOSE(v), DEBUG(d) {}
  
  double
  ViterbiDecoding(const vector<size_t> &widths,
		  const vector<size_t> &reset_points,
		  const vector<double> &start_trans, 
		  const vector<vector<double> > &trans, 
		  const vector<double> &end_trans, 
		  const double fg_lambda,
		  const double bg_lambda,
		  vector<bool> &ml_classes) const;
  
  
  double
  BaumWelchTraining(const vector<size_t> &widths,
		    const vector<size_t> &reset_points,
		    vector<double> &start_trans,
		    vector<vector<double> > &trans, 
		    vector<double> &end_trans,
		    double &fg_lambda, double &bg_lambda) const;
  
  double
  PosteriorDecoding(const vector<size_t> &widths,
		    const vector<size_t> &reset_points,
		    const vector<double> &start_trans, 
		    const vector<vector<double> > &trans, 
		    const vector<double> &end_trans, 
		    const double fg_lambda,
		    const double bg_lambda,
		    vector<bool> &classes,
		    vector<double> &llr_scores) const;
  
  void
  PosteriorScores(const vector<size_t> &widths,
		  const vector<size_t> &reset_points,
		  const vector<double> &start_trans, 
		  const vector<vector<double> > &trans, 
		  const vector<double> &end_trans, 
		  const double fg_lambda,
		  const double bg_lambda,
		  const vector<bool> &classes,
		  vector<double> &llr_scores) const;
  
  void
  PosteriorScores(const vector<size_t > &widths,
		  const vector<size_t> &reset_points,
		  const vector<double> &start_trans, 
		  const vector<vector<double> > &trans, 
		  const vector<double> &end_trans, 
		  const double fg_lambda,
		  const double bg_lambda,
		  const bool class_id,
		  vector<double> &llr_scores) const;
  
  void
  TransitionPosteriors(const vector<size_t > &widths,
		       const vector<size_t> &reset_points,
		       const vector<double> &start_trans, 
		       const vector<vector<double> > &trans, 
		       const vector<double> &end_trans, 
		       const double fg_lambda,
		       const double bg_lambda,
		       const size_t transition,
		       vector<double> &scores) const;
  
  std::string
  error_log() const;
  
  static const size_t FG_TO_BG_TRANSITION = 1;
  static const size_t BG_TO_FG_TRANSITION = 2;
    
private:

  double
  ViterbiDecoding(const vector<size_t > &widths,
		  const vector<size_t> &reset_points,
		  double p_sf, double p_sb,
		  double p_ff, double p_fb, double p_ft,
		  double p_bf, double p_bb, double p_bt,
		  const geometric &fg_distro,
		  const geometric &bg_distro,
		  vector<bool> &ml_classes) const;

  double
  BaumWelchTraining(const vector<size_t > &widths,
		    const vector<size_t> &reset_points,
		    double &p_sf, double &p_sb,
		    double &p_ff, double &p_fb, double &p_ft,
		    double &p_bf, double &p_bb, double &p_bt,
		    geometric &fg_distro,
		    geometric &bg_distro) const;

  double
  PosteriorDecoding(const vector<size_t > &widths,
		    const vector<size_t> &reset_points,
		    double p_sf, double p_sb,
		    double p_ff, double p_fb, double p_ft,
		    double p_bf, double p_bb, double p_bt,
		    const geometric &fg_distro,
		    const geometric &bg_distro,
		    vector<bool> &classes,
		    vector<double> &llr_scores) const;

  void
  PosteriorScores(const vector<size_t > &widths,
		  const vector<size_t> &reset_points,
		  double p_sf, double p_sb,
		  double p_ff, double p_fb, double p_ft,
		  double p_bf, double p_bb, double p_bt,
		  const geometric &fg_distro,
		  const geometric &bg_distro,
		  const vector<bool> &classes,
		  vector<double> &llr_scores) const;

  void
  PosteriorScores(const vector<size_t > &widths,
		  const vector<size_t> &reset_points,
		  double p_sf, double p_sb,
		  double p_ff, double p_fb, double p_ft,
		  double p_bf, double p_bb, double p_bt,
		  const geometric &fg_distro,
		  const geometric &bg_distro,
		  const bool class_id,
		  vector<double> &llr_scores) const;

  void
  TransitionPosteriors(const vector<size_t > &widths,
		       const vector<size_t> &reset_points,
		       double p_sf, double p_sb,
		       double p_ff, double p_fb, double p_ft,
		       double p_bf, double p_bb, double p_bt,
		       const geometric &fg_distro,
		       const geometric &bg_distro,
		       const size_t transition,
		       vector<double> &scores) const;
  
  double
  single_iteration(const vector<size_t > &vals,
		   const vector<size_t> &reset_points,
		   vector<std::pair<double, double> > &forward,
		   vector<std::pair<double, double> > &backward,
		   double &p_sf, double &p_sb,
		   double &p_ff, double &p_fb, double &p_ft,
		   double &p_bf, double &p_bb, double &p_bt,
		   geometric &fg_distro,
		   geometric &bg_distro) const;

  double 
  forward_algorithm(const vector<size_t > &vals,
		    const size_t start, const size_t end,
		    const double lp_sf, const double lp_sb,
		    const double lp_ff, const double lp_fb, const double lp_ft,
		    const double lp_bf, const double lp_bb, const double lp_bt,
		    const geometric &fg_distro,
		    const geometric &bg_distro,
		    vector<std::pair<double, double> > &f) const;
  double 
  backward_algorithm(const vector<size_t > &vals,
		     const size_t start, const size_t end,
		     const double lp_sf, const double lp_sb,
		     const double lp_ff, const double lp_fb, const double lp_ft,
		     const double lp_bf, const double lp_bb, const double lp_bt,
		     const geometric &fg_distro,
		     const geometric &bg_distro,
		     vector<std::pair<double, double> > &b) const;

  double
  log_sum_log_vec(const vector<double> &vals, size_t limit) const;
  
  void 
  estimate_emissions(const vector<std::pair<double, double> > &f,
		     const vector<std::pair<double, double> > &b,
		     vector<double> &fg_probs,
		     vector<double> &bg_probs) const;
  
  void
  estimate_transitions(const vector<size_t > &vals,
		       const size_t start, const size_t end,
		       const vector<std::pair<double, double> > &f,
		       const vector<std::pair<double, double> > &b,
		       const double total,
		       const geometric &fg_distro,
		       const geometric &bg_distro,
		       const double lp_ff, const double lp_fb,
		       const double lp_bf, const double lp_bb,
		       const double lp_ft, const double lp_bt,
		       vector<double> &ff_vals,
		       vector<double> &fb_vals,
		       vector<double> &bf_vals,
		       vector<double> &bb_vals) const;
  
  double
  log_sum_log(const double p, const double q) const;
  
  double MIN_PROB;
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;
  
  mutable size_t emission_correction_count;
};

#endif
