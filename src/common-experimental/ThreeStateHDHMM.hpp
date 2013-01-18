/*
  Copyright (C) 2011 University of Southern California
  Authors: Andrew D. Smith, Song Qiang

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

#ifndef THREE_STATE_HDHMM_HPP
#define THREE_STATE_HDHMM_HPP

#include <memory>
#include <utility>
#include <string>
#include <vector>

#include "smithlab_utils.hpp"
#include "Distro.hpp"

enum STATE_LABELS {GAIN, SAME, LOSS};
struct Triplet { double gain, same, loss; };

class ThreeStateHDHMM {
public:
    ThreeStateHDHMM(
        const std::vector<double> &_observations,
        const std::vector<size_t> &_reset_points, 
        const double mp, const double tol,
        const size_t max_itr, const bool v,
        const size_t _MAX_LEN);

    void
    set_parameters(const Distro & _gain_emission,
                   const Distro & _same_emission,
                   const Distro & _loss_emission,
                   const Distro & _gain_duration,
                   const Distro & _same_duration,
                   const Distro & _loss_duration,
                   const std::vector<std::vector<double> > & _trans);
    void
    get_parameters(Distro & _gain_emission,
                   Distro & _same_emission,
                   Distro & _loss_emission,
                   Distro & _gain_duration,
                   Distro & _same_duration,
                   Distro & _loss_duration,
                   std::vector<std::vector<double> > & _trans) const;
    
    double
    BaumWelchTraining();

    double
    PosteriorDecoding();

    void
    get_posterior_scores(std::vector<Triplet> &scores,
                         std::vector<STATE_LABELS>  &classes);
    
private:

    //////////// methods ////////////
    double
    single_iteration();
    double 
    forward_algorithm(const size_t start, const size_t end);
    double 
    backward_algorithm(const size_t start, const size_t end);

    double
    gain_segment_log_likelihood(const size_t start, const size_t end);

    double
    same_segment_log_likelihood(const size_t start, const size_t end);

    double
    loss_segment_log_likelihood(const size_t start, const size_t end);

    void
    estimate_state_posterior(const size_t start, const size_t end);

    void estimate_parameters();

    void update_observation_likelihood();

    ////////   data   ////////
    std::vector<double> observations;
    std::vector<size_t> reset_points;
    std::vector<double> meth_lp, unmeth_lp;
    std::vector<double> gain_log_likelihood, same_log_likelihood, loss_log_likelihood;

    //  HMM internal data 
    Distro gain_emission, same_emission, loss_emission;
    Distro gain_duration, same_duration, loss_duration;
    
    Triplet lp_start, lp_end;
    std::vector<std::vector<double> > trans;
    
    std::vector<Triplet> forward;
    std::vector<Triplet> backward;
    std::vector<double> gain_posteriors, same_posteriors, loss_posteriors;

    // parameters
    double MIN_PROB;
    double tolerance;
    size_t max_iterations;
    bool VERBOSE;
    size_t MAX_LEN;
};
// }

#endif
