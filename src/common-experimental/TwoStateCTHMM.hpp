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

#ifndef TWO_STATE_CTHMM_HPP
#define TWO_STATE_CTHMM_HPP

#include "smithlab_utils.hpp"
#include <memory>

struct betabin;
struct prob_mat;

class TwoStateCTHMM {
public:

  TwoStateCTHMM(const double ds, const double mp, const double tol,
                const size_t max_itr, const bool v, bool d = false) :
    desert_size(ds), MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),
    VERBOSE(v), DEBUG(d) {}

  double
  ViterbiDecoding(const std::vector<uint32_t> &pos,
                  const std::vector<std::pair<double, double> > &values,
                  const std::vector<size_t> &reset_points,
                  const std::vector<double> &start_trans,
                  const std::vector<std::vector<double> > &trans,
                  const std::vector<double> &end_trans,
                  const double fg_alpha, const double fg_beta,
                  const double bg_alpha, const double bg_beta,
                  std::vector<bool> &ml_classes) const;


  double
  BaumWelchTraining(const std::vector<uint32_t> &pos,
                    const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    std::vector<double> &start_trans,
                    double &mu0, double &mu1,
                    std::vector<double> &end_trans,
                    double &fg_alpha, double &fg_beta,
                    double &bg_alpha, double &bg_beta) const;

  double
  PosteriorDecoding(const std::vector<uint32_t> &pos,
                    const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    const std::vector<double> &start_trans,
                    const double mu0, const double mu1,
                    const std::vector<double> &end_trans,
                    const double fg_alpha, const double fg_beta,
                    const double bg_alpha, const double bg_beta,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;

  std::string
  error_log() const;

  static const size_t FG_TO_BG_TRANSITION = 1;
  static const size_t BG_TO_FG_TRANSITION = 2;

private:

  double
  BaumWelchTraining(const std::vector<uint32_t> &pos,
                    const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    double &p_sf, double &p_sb,
                    double &mu0, double &mu1,
                    double &p_ft, double &p_bt,
                    betabin &fg_distro, betabin &bg_distro) const;

  double
  PosteriorDecoding(const std::vector<uint32_t> &pos,
                    const std::vector<std::pair<double, double> > &values,
                    const std::vector<size_t> &reset_points,
                    const double p_sf, const double p_sb,
                    const double mu0, const double mu1,
                    const double p_ft, const double p_bt,
                    const betabin &fg_distro,
                    const betabin &bg_distro,
                    std::vector<bool> &classes,
                    std::vector<double> &llr_scores) const;

  double
  single_iteration(const std::vector<uint32_t> &pos,
                   const std::vector<std::pair<double, double> > &values,
                   const std::vector<double> &vals_a,
                   const std::vector<double> &vals_b,
                   const std::vector<size_t> &reset_points,
                   std::vector<std::pair<double, double> > &forward,
                   std::vector<std::pair<double, double> > &backward,
                   double &p_sf, double &p_sb,
                   double &mu0, double &mu1,
                   double &p_ft, double &p_bt,
                   betabin &fg_distro, betabin &bg_distro) const;

  double
  forward_algorithm(const std::vector<uint32_t> &pos,
                  const std::vector<std::pair<double, double> > &vals,
                    const size_t start, const size_t end,
                    const double lp_sf, const double lp_sb,
                    const std::vector<prob_mat> &lm,
                    const double lp_ft, const double lp_bt,
                    const betabin &fg_distro,
                    const betabin &bg_distro,
                    std::vector<std::pair<double, double> > &f) const;
  double
  backward_algorithm(const std::vector<uint32_t> &pos,
                     const std::vector<std::pair<double, double> > &vals,
                     const size_t start, const size_t end,
                     const double lp_sf, const double lp_sb,
                     const std::vector<prob_mat> &lm,
                     const double lp_ft, const double lp_bt,
                     const betabin &fg_distro,
                     const betabin &bg_distro,
                     std::vector<std::pair<double, double> > &b) const;

  double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;

  void
  estimate_emissions(const std::vector<std::pair<double, double> > &f,
                     const std::vector<std::pair<double, double> > &b,
                     std::vector<double> &fg_probs,
                     std::vector<double> &bg_probs) const;

  void
  estimate_transitions(const std::vector<uint32_t> &pos,
                  const std::vector<std::pair<double, double> > &vals,
                       const size_t start, const size_t end,
                       const std::vector<std::pair<double, double> > &f,
                       const std::vector<std::pair<double, double> > &b,
                       const double total,
                       const betabin &fg_distro,
                       const betabin &bg_distro,
                       const std::vector<prob_mat> &lm,
                       std::vector<double> &ff_vals,
                       std::vector<double> &fb_vals,
                       std::vector<double> &bf_vals,
                       std::vector<double> &bb_vals) const;

  double
  log_sum_log(const double p, const double q) const;

  uint32_t desert_size;
  double MIN_PROB;
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;
};

#endif
