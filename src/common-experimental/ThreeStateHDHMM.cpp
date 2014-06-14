/*
  Copyright (C) 2011 University of Southern California
  Authors: Song Qiang

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

// #define DEBUG

#include "ThreeStateHDHMM.hpp"
#include "numerical_utils.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

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


static STATE_LABELS 
get_state(const double gain, const double same, const double loss)
{
    if (gain > same && gain > loss)
        return GAIN;
    else if (loss > gain && loss > same)
        return LOSS;
    else
        return SAME;
}

ThreeStateHDHMM::ThreeStateHDHMM(
    const std::vector<double> &_observations,
    const std::vector<size_t> &_reset_points, 
    const double tol,
    const size_t max_itr, const bool v,
    const size_t _MAX_LEN):
    observations(_observations), reset_points(_reset_points),
    meth_lp(_observations.size()), unmeth_lp(_observations.size()),
    gain_log_likelihood(_observations.size()),
    same_log_likelihood(_observations.size()),
    loss_log_likelihood(_observations.size()),
    trans(3, vector<double>(3, 0)),
    forward(_observations.size()), backward(_observations.size()),
    gain_posteriors(_observations.size()),
    same_posteriors(_observations.size()),
    loss_posteriors(_observations.size()),
    tolerance(tol), max_iterations(max_itr),
    VERBOSE(v), MAX_LEN(_MAX_LEN)
{
    for (size_t i = 0; i < observations.size(); ++i)
    {
        meth_lp[i] = 
            log(std::min(std::max(observations[i], 1e-10), 1.0 - 1e-10));
        unmeth_lp[i] = 
            log(std::min(std::max(1 - observations[i], 1e-10), 1.0 - 1e-10));
        observations[i] =
            std::min(std::max(observations[i], 1e-10), 1.0 - 1e-10);
    }
}

void
ThreeStateHDHMM::set_parameters(
    const Distro & _gain_emission,
    const Distro & _same_emission,
    const Distro & _loss_emission,
    const Distro & _gain_duration,
    const Distro & _same_duration,
    const Distro & _loss_duration,
    const vector<vector<double> > & _trans)
{
#ifdef DEBUG
    cerr << "check enter initialize_parameters: "<< "OK" << endl;
#endif
    
    // initialize emission parameters
    gain_emission = _gain_emission;
    same_emission = _same_emission;
    loss_emission = _loss_emission;
    update_observation_likelihood();

    // initialize duration parameters
    gain_duration = _gain_duration;
    same_duration = _same_duration;
    loss_duration = _loss_duration;

    // initialize transition matrix 
    trans  = _trans;

#ifdef DEBUG
    cerr << "check exit initialize_parameters: "<< "OK" << endl;
#endif
}

void
ThreeStateHDHMM::get_parameters(
    Distro & _gain_emission,
    Distro & _same_emission,
    Distro & _loss_emission,
    Distro & _gain_duration,
    Distro & _same_duration,
    Distro & _loss_duration,
    vector<vector<double> > & _trans) const
{
    _gain_emission = gain_emission;
    _same_emission = same_emission;
    _loss_emission = loss_emission;

    _gain_duration = gain_duration;
    _same_duration = same_duration;
    _loss_duration = loss_duration;

    _trans  = trans;
}


//////////////////////////////////////////////
////// forward and backward algorithms  //////
//////////////////////////////////////////////
void
ThreeStateHDHMM::update_observation_likelihood()
{
    assert(observations.size() == gain_log_likelihood.size());
    assert(observations.size() == same_log_likelihood.size());
    assert(observations.size() == loss_log_likelihood.size());
    
    gain_log_likelihood.front() =
        gain_emission.log_likelihood(observations.front());
    same_log_likelihood.front() =
        same_emission.log_likelihood(observations.front());
    loss_log_likelihood.front() =
        loss_emission.log_likelihood(observations.front());
    
    for (size_t i = 1; i < observations.size(); ++i)
    {
        gain_log_likelihood[i] = gain_log_likelihood[i - 1]
            + gain_emission.log_likelihood(observations[i]);
        same_log_likelihood[i] = same_log_likelihood[i - 1]
            + same_emission.log_likelihood(observations[i]);
        loss_log_likelihood[i] = loss_log_likelihood[i - 1]
            + loss_emission.log_likelihood(observations[i]);

        assert(isfinite(gain_log_likelihood[i]));
        assert(isfinite(same_log_likelihood[i]));
        assert(isfinite(loss_log_likelihood[i]));
    }
}

double
ThreeStateHDHMM::gain_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? gain_log_likelihood[end - 1]
        : gain_log_likelihood[end - 1] - gain_log_likelihood[start - 1];
}

double
ThreeStateHDHMM::same_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? same_log_likelihood[end - 1]
        : same_log_likelihood[end - 1] - same_log_likelihood[start - 1];
}

double
ThreeStateHDHMM::loss_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? loss_log_likelihood[end - 1]
        : loss_log_likelihood[end - 1] - loss_log_likelihood[start - 1];
}

double
ThreeStateHDHMM::forward_algorithm(const size_t start, const size_t end) 
{
#ifdef DEBUG
    cerr << "check enter forward_algorithm: "<< "OK" << endl;
#endif
    assert(start < end);
    
    for (size_t i = start; i < end; ++i)
        forward[i].gain = forward[i].same = forward[i].loss = 0.0;

    forward[start].gain =
        lp_start.gain + gain_segment_log_likelihood(start, start + 1)
        + gain_duration.log_likelihood(1);
    forward[start].same =
        lp_start.same + same_segment_log_likelihood(start, start + 1);
    forward[start].loss =
        lp_start.loss + loss_segment_log_likelihood(start, start + 1)
        + loss_duration.log_likelihood(1);

    for (size_t i = start + 1; i < end; ++i)
    {
        // in segment gaining methylation
        const size_t max_len = min(i - start + 1, MAX_LEN) + 1;
        for (size_t l = 1; l < max_len; ++l)
        {
            const size_t beginning = i - l + 1; // inclusive
            const size_t ending = i + 1;       // exclusive
            
            // transition from no-change segment
            const double gain_seg_llh =  (beginning == start) ?
                lp_start.gain  
                + gain_segment_log_likelihood(beginning, ending)
                + gain_duration.log_likelihood(l)
                :
                forward[beginning - 1].same + log(trans[SAME][GAIN])
                + gain_segment_log_likelihood(beginning, ending)
                + gain_duration.log_likelihood(l);
            
            forward[i].gain = log_sum_log(forward[i].gain, gain_seg_llh);
        }

        // in non-change segment 
        forward[i].same =
            log_sum_log(forward[i-1].gain + log(trans[GAIN][SAME]),
                        forward[i-1].loss + log(trans[LOSS][SAME]),
                        forward[i-1].same + log(trans[SAME][SAME]))
            + same_segment_log_likelihood(i, i + 1);
        
        // in segment losing methylation
        for (size_t l = 1; l < max_len; ++l)
        {
            const size_t beginning = i - l + 1; // inclusive
            const size_t ending = i + 1;       // exclusive
            
            // transition from no-change segment
            const double loss_seg_llh =  (beginning == start) ?
                lp_start.loss  
                + loss_segment_log_likelihood(beginning, ending)
                + loss_duration.log_likelihood(l)
                :
                forward[beginning - 1].same + log(trans[SAME][LOSS])
                + loss_segment_log_likelihood(beginning, ending)
                + loss_duration.log_likelihood(l);
            
            forward[i].loss = log_sum_log(forward[i].loss, loss_seg_llh);
        }
        
        assert(isfinite(forward[i].gain));
        assert(isfinite(forward[i].same));
        assert(isfinite(forward[i].loss));
    }

#ifdef DEBUG
    cerr << "check exit forward_algorithm: "<< "OK" << endl;
#endif

    return log_sum_log(forward[end - 1].gain + lp_end.gain,
                       forward[end - 1].same + lp_end.same,
                       forward[end - 1].loss + lp_end.loss);
}

double
ThreeStateHDHMM::backward_algorithm(const size_t start, const size_t end)
{
#ifdef DEBUG
    cerr << "check enter backward_algorithm: "<< "OK" << endl;
#endif
    const int start_int(start), end_int(end);
    
    for (int i = start_int; i < end_int; ++i)
        backward[i].gain = backward[i].same = backward[i].loss = 0.0;
 
    backward[end - 1].gain = lp_end.gain;
    backward[end - 1].same = lp_end.same;
    backward[end - 1].loss = lp_end.loss;

    
    for (int i = end_int - 2; i >= start_int; --i)
    {
        assert(i >= start_int && i < end_int);
        // Terminate a segment gaining methylation
        backward[i].gain =
            log(trans[GAIN][SAME])
            + same_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].same;

        // Terminate a segment losing methylation
        backward[i].loss =
            log(trans[LOSS][SAME])
            + same_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].same;

        // Remain at a no-change segment
        backward[i].same =
            log(trans[SAME][SAME])
            + same_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].same;

        // Terminate a no-change segment
        const size_t max_len = min(end - i - 1, MAX_LEN) + 1;
        for (size_t l = 1; l < max_len; ++l)
        {
            const size_t beginning = i + 1;
            const size_t ending = i + l + 1;

            const double next_seg_llh = 
                log_sum_log(
                    log(trans[SAME][GAIN])  // transite to GAINing segments
                    + gain_segment_log_likelihood(beginning, ending)
                    + gain_duration.log_likelihood(l)
                    + backward[ending - 1].gain,
 
                    log(trans[SAME][LOSS]) // transite to LOSSing segments 
                    + loss_segment_log_likelihood(beginning, ending)
                    + loss_duration.log_likelihood(l)
                    + backward[ending - 1].loss);
            backward[i].same = log_sum_log(backward[i].same, next_seg_llh);
        }


        assert(isfinite(backward[i].gain));
        assert(isfinite(backward[i].same));
        assert(isfinite(backward[i].loss));

    }
    
    // whole likelihood
    double llh = 0;

    // the segment is no-change
    const double same_seg_llh =
        lp_start.same  
        + same_segment_log_likelihood(start, start + 1)
        + backward[start].same;
    llh = log_sum_log(llh, same_seg_llh);

    const size_t max_len = min(end - start, MAX_LEN) + 1;
    for (size_t l = 1; l < max_len; ++l)
    {
        const size_t beginning = start;
        const size_t ending = start + l;

        // the segment is gaining methylation
        const double gain_seg_llh =
            lp_start.gain  
            + gain_segment_log_likelihood(beginning, ending)
            + gain_duration.log_likelihood(l)
            + backward[ending - 1].gain;
        llh = log_sum_log(llh, gain_seg_llh);

        // the segment is lossing methylation
        const double loss_seg_llh =
            lp_start.loss  
            + loss_segment_log_likelihood(beginning, ending)
            + loss_duration.log_likelihood(l)
            + backward[ending - 1].loss;
        llh = log_sum_log(llh, loss_seg_llh);
    }

#ifdef DEBUG
    cerr << "check exit backward_algorithm: "<< "OK" << endl;
#endif
    return llh;

}

//////////////////////////////////////////////
//////       Baum-Welch Training        //////
//////////////////////////////////////////////
// Expectation
void
ThreeStateHDHMM::estimate_state_posterior(const size_t start, const size_t end)  
{
#ifdef DEBUG
    cerr << "check enter estimate_state_posterior: "<< "OK" << endl;
#endif

    vector<double> gain_evidence(end - start, 0), same_evidence(end - start, 0),
        loss_evidence(end - start, 0);

    // Evidence for GAIN state
    for (size_t s = start; s < end; ++s)
    {
        // foreground
        double accu_evidence = 0;
        for (size_t e =  min(s + MAX_LEN, end); e > s; --e)
        {
            const double evidence = (s == start) ?
                lp_start.gain
                + gain_duration.log_likelihood(e - s)
                + gain_segment_log_likelihood(s, e)
                + backward[e - 1].gain
                :
                forward[s - 1].same + log(trans[SAME][GAIN]) 
                + gain_duration.log_likelihood(e - s)
                + gain_segment_log_likelihood(s, e)
                + backward[e - 1].gain;
            
            accu_evidence = log_sum_log(accu_evidence, evidence);
            gain_evidence[e - 1 - start] =
                log_sum_log(gain_evidence[e - 1 - start], accu_evidence);
        }
        assert(isfinite(gain_evidence[s - start]));
    }

    // Evidence for no-change state
    for (size_t i = start; i < end; ++i)
    {
        same_evidence[i - start] = forward[i].same + backward[i].same;
        assert(isfinite(same_evidence[i - start]));
    }
    
    // Evidence for LOSS state
    for (size_t s = start; s < end; ++s)
    {
        // foreground
        double accu_evidence = 0;
        for (size_t e =  min(s + MAX_LEN, end); e > s; --e)
        {
            const double evidence = (s == start) ?
                lp_start.loss
                + loss_duration.log_likelihood(e - s)
                + loss_segment_log_likelihood(s, e)
                + backward[e - 1].loss
                :
                forward[s - 1].same + log(trans[SAME][LOSS]) 
                + loss_duration.log_likelihood(e - s)
                + loss_segment_log_likelihood(s, e)
                + backward[e - 1].loss;
            
            accu_evidence = log_sum_log(accu_evidence, evidence);
            loss_evidence[e - 1 - start] =
                log_sum_log(loss_evidence[e - 1 - start], accu_evidence);
        }
        assert(isfinite(loss_evidence[s - start]));
    }

    // state posterior
    for (size_t i = start; i < end; ++i)
    {
        const double denom =
            log_sum_log(gain_evidence[i - start],
                        same_evidence[i - start],
                        loss_evidence[i - start]);

        gain_posteriors[i] = exp(gain_evidence[i - start] - denom);
        same_posteriors[i] = exp(same_evidence[i - start] - denom);
        loss_posteriors[i] = exp(loss_evidence[i - start] - denom);

        assert(fabs(gain_posteriors[i] + same_posteriors[i]
                    + loss_posteriors[i] - 1.0) < 1e-6);
    }

#ifdef DEBUG
    cerr << "check exit estimate_state_posterior: "<< "OK" << endl;
#endif
}

void 
ThreeStateHDHMM::estimate_parameters()
{
#ifdef DEBUG
    cerr << "check enter estimate_parameters: "<< "OK" << endl;
#endif

    // estimate emission parameters 
    gain_emission.estimate_params_ml(meth_lp, unmeth_lp, gain_posteriors);
    same_emission.estimate_params_ml(meth_lp, unmeth_lp, same_posteriors);
    loss_emission.estimate_params_ml(meth_lp, unmeth_lp, loss_posteriors);
    update_observation_likelihood();

    // estimate transition matrix 
    vector<double> gain_lengths, same_lengths, loss_lengths;
    for (size_t idx = 0; idx < reset_points.size() - 1; ++idx)
    {
        const size_t start = reset_points[idx];
        const size_t end = reset_points[idx + 1];
        
        STATE_LABELS prev_state =
            get_state(gain_posteriors[start], same_posteriors[start],
                 loss_posteriors[start]);
        size_t len = 1;
        for (size_t i = start + 1; i < end; ++i)
        {
            const STATE_LABELS state = 
                get_state(gain_posteriors[i], same_posteriors[i], loss_posteriors[i]);
            if (state == prev_state)
                ++len;
            else
            {
                switch(prev_state)
                {
                case GAIN: gain_lengths.push_back(len); break;
                case SAME: same_lengths.push_back(len); break;
                case LOSS: loss_lengths.push_back(len); break;
                }

                prev_state = state;
                len = 1;
            }
        }
    }

    if (gain_lengths.size() > 0)
        gain_duration.estimate_params_ml(gain_lengths);

    if (same_lengths.size() > 0)
        same_duration.estimate_params_ml(same_lengths);

    if (loss_lengths.size() > 0)
        loss_duration.estimate_params_ml(loss_lengths);


    // estiamting transition probabilities
    trans[SAME][SAME] = 1 - same_duration.get_params().front();
    assert(isfinite(log(trans[SAME][SAME])));

    trans[SAME][GAIN] = same_duration.get_params().front()
        * gain_lengths.size()
        / (gain_lengths.size() + loss_lengths.size());
    assert(isfinite(log(trans[SAME][GAIN])));

    trans[SAME][LOSS] = same_duration.get_params().front()
        * loss_lengths.size()
        / (gain_lengths.size() + loss_lengths.size());
    assert(isfinite(log(trans[SAME][LOSS])));

#ifdef DEBUG
    cerr << "check exit estimate_parameters: "<< "OK" << endl;
#endif
}

double
ThreeStateHDHMM::single_iteration()
{
#ifdef DEBUG
    cerr << "check enter single_iteration: "<< "OK" << endl;
#endif

    double total_score = 0;
    
    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double forward_score =
            forward_algorithm(reset_points[i], reset_points[i + 1]);
        const double backward_score =
            backward_algorithm(reset_points[i], reset_points[i + 1]);
        
        if (fabs((forward_score - backward_score)
                 / max(forward_score, backward_score))
            > 1e-10)
        {
            cerr << "[WARNING] single_iteration: forward_score="
                 << forward_score << ", "
                 << "backward_score=" << backward_score << ", "
                 << fabs((forward_score - backward_score)
                         / max(forward_score, backward_score)) << endl;
            assert(fabs((forward_score - backward_score)
                        / max(forward_score, backward_score))
                   < 1e-6);
        }
        estimate_state_posterior(reset_points[i], reset_points[i + 1]);
        total_score += forward_score;
    }

    estimate_parameters();

#ifdef DEBUG
    cerr << "check exit single_iteration: "<< "OK" << endl;
#endif

    return total_score;
}

double
ThreeStateHDHMM::BaumWelchTraining()  
{
#ifdef DEBUG
// /////
    cerr << "check enter BaumWelchTraining: "<< "OK" << endl;
// /////
#endif

    if (VERBOSE) cerr << "Model training:" << endl;
  
    double prev_total = -std::numeric_limits<double>::max();
  
    for (size_t i = 0; i < max_iterations; ++i) 
    {
        const Distro old_gain_emission = gain_emission;
        const Distro old_same_emission = same_emission;
        const Distro old_loss_emission = loss_emission;

        const Distro old_gain_duration = gain_duration;
        const Distro old_same_duration = same_duration;
        const Distro old_loss_duration = loss_duration;

        const vector<vector<double> > old_trans(trans);
        
        double total = single_iteration();
    
        if (VERBOSE)
        {
            cerr << "Iteration: " << i + 1 << ";\t" 
                 << "Log-Likelihood: " << total << ";\t"
                 << "Change: " << (total - prev_total)/std::fabs(total) 
                 << endl
                 << "GAIN:\t" << setw(18) << old_gain_emission.tostring()
                 << setw(18) << old_gain_duration.tostring()
                 << endl
                 << "SAME:\t" << setw(18) << old_same_emission.tostring()
                 << setw(18) << old_same_duration.tostring()
                 << endl
                 << "LOSS:\t" << setw(18) << old_loss_emission.tostring()
                 << setw(18) << old_loss_duration.tostring()
                 << endl;

            cerr << setw(4) << "" << setw(10) << "GAIN"
                 << setw(10) << "SAME" << setw(10) << "LOSS" << endl;
            for (size_t r = 0; r < 3; ++r)
            {
                switch (r)
                {
                case 0: cerr << setw(4) << "GAIN"; break;
                case 1: cerr << setw(4) << "SAME"; break;
                case 2: cerr << setw(4) << "LOSS"; break;
                } 

                for (size_t c = 0; c < 3; ++c)
                    cerr << setw(10) << old_trans[r][c];
                cerr << endl;
            }
            cerr << endl;
        }
        
        if ((total - prev_total)/std::fabs(total) < tolerance)
        {
            gain_emission = old_gain_emission;
            same_emission = old_same_emission;
            loss_emission = old_loss_emission;
            update_observation_likelihood();
            
            gain_duration = old_gain_duration;
            same_duration = old_same_duration;
            loss_duration = old_loss_duration;

            if (VERBOSE)
                cerr << "CONVERGED" << endl << endl;
            break;
        }
        prev_total = total;
    }

#ifdef DEBUG
    cerr << "check BaumWelchTraining: "<< "OK" << endl;
#endif

    return prev_total;
}

//////////////////////////////////////////////
//////          Posterior decoding      //////
//////////////////////////////////////////////
double
ThreeStateHDHMM::PosteriorDecoding()
{
#ifdef DEBUG
    cerr << "check enter PosteriorDecoding: "<< "OK" << endl;
#endif

    double total_score = 0;
    
    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double forward_score =
            forward_algorithm(reset_points[i], reset_points[i + 1]);
        const double backward_score =
            backward_algorithm(reset_points[i], reset_points[i + 1]);
        
        if (fabs((forward_score - backward_score)
                 / max(forward_score, backward_score))
            >=  1e-10)
        {
            cerr << "[WARNING] PosteriorDecoding: forward_score="
                 << forward_score << ", "
                 << "backward_score=" << backward_score << ", "
                 << fabs((forward_score - backward_score)
                         / max(forward_score, backward_score)) << endl;
            assert(fabs((forward_score - backward_score)
                        / max(forward_score, backward_score))
                   < 1e-6);
        }
        estimate_state_posterior(reset_points[i], reset_points[i + 1]);
        total_score += forward_score;
    }

#ifdef DEBUG
    cerr << "check exit PosteriorDecoding: "<< "OK" << endl;
#endif
    return total_score;
}

//////////////////////////////////////////////
//////          export result           //////
//////////////////////////////////////////////

void
ThreeStateHDHMM::get_posterior_scores(
    std::vector<Triplet> & scores,
    std::vector<STATE_LABELS> & classes)
{
    scores.resize(gain_posteriors.size());
    for (size_t i = 0; i < gain_posteriors.size(); ++i)
    {
        scores[i].gain = gain_posteriors[i];
        scores[i].same = same_posteriors[i];
        scores[i].loss = loss_posteriors[i];
    }

    classes.resize(observations.size());
    for (size_t i = 0; i < observations.size(); ++i)
        classes[i] = get_state(gain_posteriors[i],
                          same_posteriors[i],
                          loss_posteriors[i]);
}



