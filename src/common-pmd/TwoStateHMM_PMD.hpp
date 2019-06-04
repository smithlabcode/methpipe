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

#ifndef TWO_STATE_HMM_PMD_HPP
#define TWO_STATE_HMM_PMD_HPP
#include <memory>
#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include "smithlab_utils.hpp"
#include "EmissionDistribution.hpp"

class TwoStateHMM {
public:

  TwoStateHMM(const double mp, const double tol,
	       const size_t max_itr, const bool v, bool d = false) :
    MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),
    VERBOSE(v), DEBUG(d) {}

  /***************************/
  /* for multiple replicates */
  double
  BaumWelchTraining_rep(
      const std::vector<std::vector<std::pair<double, double> > > &values,
			const std::vector<size_t> &reset_points,
			std::vector<double> &start_trans,
			std::vector<std::vector<double> > &trans,
			std::vector<double> &end_trans,
			std::vector<double> &fg_alpha, std::vector<double> &fg_beta,
			std::vector<double> &bg_alpha, std::vector<double> &bg_beta,
      const std::vector<bool> &array_status) const;

  double
  PosteriorDecoding_rep(
      const std::vector<std::vector<std::pair<double, double> > > &values,
			const std::vector<size_t> &reset_points,
			const std::vector<double> &start_trans,
			const std::vector<std::vector<double> > &trans,
			const std::vector<double> &end_trans,
			const std::vector<double> fg_alpha, const std::vector<double> fg_beta,
			const std::vector<double> bg_alpha, const std::vector<double> bg_beta,
			std::vector<bool> &classes,
			std::vector<double> &llr_scores,
      const std::vector<bool> &array_status) const;

  void
  PosteriorScores_rep(
          const std::vector<std::vector<std::pair<double, double> > > &values,
		      const std::vector<size_t> &reset_points,
		      const std::vector<double> &start_trans,
		      const std::vector<std::vector<double> > &trans,
		      const std::vector<double> &end_trans,
		      const std::vector<double> fg_alpha, const std::vector<double> fg_beta,
		      const std::vector<double> bg_alpha, const std::vector<double> bg_beta,
		      const std::vector<bool> &classes,
		      std::vector<double> &llr_scores,
          const std::vector<bool> &array_status) const;

  void
  TransitionPosteriors_rep(
          const std::vector<std::vector<std::pair<double, double> > > &values,
          const std::vector<size_t> &reset_points,
          const std::vector<double> &start_trans,
          const std::vector<std::vector<double> > &trans,
          const std::vector<double> &end_trans,
          const std::vector<double> &fg_alpha,
          const std::vector<double> &fg_beta,
          const std::vector<double> &bg_alpha,
          const std::vector<double> &bg_beta,
          const std::vector<bool> &array_status,
          const size_t transition,
          std::vector<double> &scores) const;

  /***************************/

  std::string
  error_log() const;

  static const size_t FG_TO_BG_TRANSITION = 1;
  static const size_t BG_TO_FG_TRANSITION = 2;

private:

  double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;

  void
  estimate_emissions(const std::vector<std::pair<double, double> > &f,
		     const std::vector<std::pair<double, double> > &b,
		     std::vector<double> &fg_probs,
		     std::vector<double> &bg_probs) const;

  void
  estimate_transitions(const std::vector<std::pair<double, double> > &vals,
		       const size_t start, const size_t end,
		       const std::vector<std::pair<double, double> > &f,
		       const std::vector<std::pair<double, double> > &b,
		       const double total,
		       const EmissionDistribution &fg_distro,
		       const EmissionDistribution &bg_distro,
		       const double lp_ff, const double lp_fb,
		       const double lp_bf, const double lp_bb,
		       const double lp_ft, const double lp_bt,
		       std::vector<double> &ff_vals,
		       std::vector<double> &fb_vals,
		       std::vector<double> &bf_vals,
		       std::vector<double> &bb_vals) const;


  /***************************/
  /* for multiple replicates */
  double
  forward_algorithm_rep(
      const std::vector<std::vector<std::pair<double, double> > > &vals,
			const size_t start, const size_t end,
			const double lp_sf, const double lp_sb,
			const double lp_ff, const double lp_fb, const double lp_ft,
			const double lp_bf, const double lp_bb, const double lp_bt,
			const std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
			const std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
			std::vector<std::pair<double, double> > &f,
      const std::vector<bool> &array_status) const;

  double
  backward_algorithm_rep(
       const std::vector<std::vector<std::pair<double, double> > > &vals,
			 const size_t start, const size_t end,
			 const double lp_sf, const double lp_sb,
			 const double lp_ff, const double lp_fb, const double lp_ft,
			 const double lp_bf, const double lp_bb, const double lp_bt,
			 const std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
			 const std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
			 std::vector<std::pair<double, double> > &b,
       const std::vector<bool> &array_status) const;

  void
  estimate_transitions_rep(
         const std::vector<std::vector<std::pair<double, double> > > &vals,
			   const size_t start, const size_t end,
			   const std::vector<std::pair<double, double> > &f,
			   const std::vector<std::pair<double, double> > &b,
			   const double total,
			   const std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
			   const std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
			   const double lp_ff, const double lp_fb,
			   const double lp_bf, const double lp_bb,
			   const double lp_ft, const double lp_bt,
			   std::vector<double> &ff_vals,
			   std::vector<double> &fb_vals,
			   std::vector<double> &bf_vals,
			   std::vector<double> &bb_vals,
         const std::vector<bool> &array_status) const;

  double
  single_iteration_rep(
           const std::vector<std::vector<std::pair<double, double> > > &values,
		       const std::vector<std::vector<double> > &vals_a,
		       const std::vector<std::vector<double> > &vals_b,
		       const std::vector<size_t> &reset_points,
		       std::vector<std::pair<double, double> > &forward,
		       std::vector<std::pair<double, double> > &backward,
		       double &p_sf, double &p_sb,
		       double &p_ff, double &p_fb, double &p_ft,
		       double &p_bf, double &p_bb, double &p_bt,
		       std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
		       std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
           const std::vector<bool> &array_status) const;


  double
  BaumWelchTraining_rep(
      const std::vector<std::vector<std::pair<double, double> > > &values,
			const std::vector<size_t> &reset_points,
			double &p_sf, double &p_sb,
			double &p_ff, double &p_fb, double &p_ft,
			double &p_bf, double &p_bb, double &p_bt,
			std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
			std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
      const std::vector<bool> &array_status) const;

  void
  PosteriorScores_rep(
          const std::vector<std::vector<std::pair<double, double> > > &values,
		      const std::vector<size_t> &reset_points,
		      double p_sf, double p_sb,
		      double p_ff, double p_fb, double p_ft,
		      double p_bf, double p_bb, double p_bt,
		      const std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
		      const std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
		      const std::vector<bool> &classes,
		      std::vector<double> &llr_scores,
          const std::vector<bool> &array_status) const;

  double
  PosteriorDecoding_rep(
      const std::vector< std::vector<std::pair<double, double> > > &values,
			const std::vector<size_t> &reset_points,
			double p_sf, double p_sb,
			double p_ff, double p_fb, double p_ft,
			double p_bf, double p_bb, double p_bt,
			const std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
			const std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
			std::vector<bool> &classes,
			std::vector<double> &llr_scores,
      const std::vector<bool> &array_status) const;

  void
  TransitionPosteriors_rep(
      const std::vector< std::vector<std::pair<double, double> > > &values,
      const std::vector<size_t> &reset_points,
      double p_sf, double p_sb,
      double p_ff, double p_fb, double p_ft,
      double p_bf, double p_bb, double p_bt,
      const std::vector<std::unique_ptr<EmissionDistribution> > &fg_distro,
      const std::vector<std::unique_ptr<EmissionDistribution> > &bg_distro,
      const size_t transition, const std::vector<bool> &array_status,
      std::vector<double> &scores) const;



  /***************************/

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
