/*
  Copyright (C) 2011 University of Southern California
  Authors: Song Qiang, Andrew D. Smith

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

#include "ThreeStateHMM.hpp"
#include "numerical_utils.hpp"
#include "BetaBin.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>
#include <stdexcept>

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
max_state(const Triplet &likelihoods)
{
    if (likelihoods.hypo >= std::max(likelihoods.HYPER, likelihoods.HYPO))
        return hypo;
    else if (likelihoods.HYPER >= std::max(likelihoods.hypo, likelihoods.HYPO))
        return HYPER;
    else
        return HYPO;
}

static double
max_value(const Triplet &likelihoods)
{
    if (likelihoods.hypo >= std::max(likelihoods.HYPER, likelihoods.HYPO))
        return likelihoods.hypo;
    else if (likelihoods.HYPER >= std::max(likelihoods.hypo, likelihoods.HYPO))
        return likelihoods.HYPER;
    else
        return likelihoods.HYPO;
}

// static STATE_LABELS
// max_state(const double hypo_v, const double HYPER_v, const double HYPO_v)
// {
//     if (hypo_v >= std::max(HYPER_v, HYPO_v))
//         return hypo;
//     else if (HYPER_v >= std::max(hypo_v, HYPO_v))
//         return HYPER;
//     else
//         return HYPO;
// }

// static double
// max_value(const double hypo_v, const double HYPER_v, const double HYPO_v)
// {
//     if (hypo_v >= std::max(HYPER_v, HYPO_v))
//         return hypo_v;
//     else if (HYPER_v >= std::max(hypo_v, HYPO_v))
//         return HYPER_v;
//     else
//         return HYPO_v;
// }

ThreeStateHMM::ThreeStateHMM(
    const std::vector<std::pair<double, double> > &_observations,
    const std::vector<size_t> &_reset_points,
    const double tol, const size_t max_itr, const bool v) :
    observations(_observations), reset_points(_reset_points),
    meth_lp(_observations.size()), unmeth_lp(_observations.size()),
    hypo_log_likelihood(_observations.size()),
    HYPER_log_likelihood(_observations.size()),
    HYPO_log_likelihood(_observations.size()),
    forward(_observations.size()), backward(_observations.size()),
    hypo_posteriors(_observations.size()),
    HYPER_posteriors(_observations.size()),
    HYPO_posteriors(_observations.size()),
    hypo_hypo(_observations.size()), hypo_HYPER(_observations.size()),
    HYPER_hypo(_observations.size()), HYPER_HYPER(_observations.size()),
    HYPER_HYPO(_observations.size()), HYPO_HYPER(_observations.size()),
    HYPO_HYPO(_observations.size()), classes(_observations.size()),
    state_posteriors(_observations.size()),
    tolerance(tol), max_iterations(max_itr),
    VERBOSE(v)
{
    for (size_t i = 0; i < observations.size(); ++i)
    {
        const double m = observations[i].first;
        const double u = observations[i].second;

        meth_lp[i] =
            log(std::min(std::max(m/(m + u), 1e-2), 1.0 - 1e-2));
        unmeth_lp[i] =
            log(std::min(std::max(u/(m + u), 1e-2), 1.0 - 1e-2));
    }
}

void
ThreeStateHMM::set_parameters(
    const betabin & _hypo_emission,
    const betabin & _HYPER_emission,
    const betabin & _HYPO_emission,
    const vector<vector<double> > &_trans)
{
    hypo_emission = _hypo_emission;
    HYPER_emission = _HYPER_emission;
    HYPO_emission = _HYPO_emission;
    update_observation_likelihood();

    lp_start.hypo = log(0.5);
    lp_start.HYPER = log(0.25);
    lp_start.HYPO = log(0.25);

    lp_end.hypo = log(1e-10);
    lp_end.HYPER = log(1e-10);
    lp_end.HYPO = log(1e-10);

    trans = _trans;
}

void
ThreeStateHMM::get_parameters(
    betabin & _hypo_emission,
    betabin & _HYPER_emission,
    betabin & _HYPO_emission,
    vector<vector<double> > &_trans) const
{
    _hypo_emission = hypo_emission;
    _HYPER_emission = HYPER_emission;
    _HYPO_emission = HYPO_emission;
    _trans = trans;
}

//////////////////////////////////////////////
////// forward and backward algorithms  //////
//////////////////////////////////////////////
void
ThreeStateHMM::update_observation_likelihood()
{
    hypo_log_likelihood.front() = hypo_emission(observations.front());
    HYPER_log_likelihood.front() = HYPER_emission(observations.front());
    HYPO_log_likelihood.front() = HYPO_emission(observations.front());

    for (size_t i = 1; i < observations.size(); ++i)
    {
        hypo_log_likelihood[i] = hypo_log_likelihood[i - 1]
            + hypo_emission(observations[i]);
        HYPER_log_likelihood[i] = HYPER_log_likelihood[i - 1]
            + HYPER_emission(observations[i]);
        HYPO_log_likelihood[i] = HYPO_log_likelihood[i - 1]
            + HYPO_emission(observations[i]);
    }
}

double
ThreeStateHMM::hypo_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? hypo_log_likelihood[end - 1]
        : hypo_log_likelihood[end - 1] - hypo_log_likelihood[start - 1];
}

double
ThreeStateHMM::HYPER_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? HYPER_log_likelihood[end - 1]
        : HYPER_log_likelihood[end - 1] - HYPER_log_likelihood[start - 1];
}

double
ThreeStateHMM::HYPO_segment_log_likelihood(
    const size_t start, const size_t end)
{
    return
        (start == 0)
        ? HYPO_log_likelihood[end - 1]
        : HYPO_log_likelihood[end - 1] - HYPO_log_likelihood[start - 1];
}

double
ThreeStateHMM::forward_algorithm(const size_t start, const size_t end)
{
/////
    // cerr << "check enter forward_algorithm: "<< "OK" << endl;
/////
    for (size_t i = start; i < end; ++i)
        forward[i].hypo = forward[i].HYPER = forward[i].HYPO = 0.0;

    forward[start].hypo =
        lp_start.hypo + hypo_segment_log_likelihood(start, start + 1);
    forward[start].HYPER =
        lp_start.HYPER + HYPER_segment_log_likelihood(start, start + 1);
    forward[start].HYPO =
        lp_start.HYPO + HYPO_segment_log_likelihood(start, start + 1);

    for (size_t i = start + 1; i < end; ++i)
    {
        // hypomethylated CpG in HypoMR segment
        forward[i].hypo =
            log_sum_log(forward[i - 1].hypo + log(trans[hypo][hypo]),
                        forward[i - 1].HYPER + log(trans[HYPER][hypo]))
            + hypo_segment_log_likelihood(i, i + 1);

        // hypermethylated CpG in HyperMR segment
        forward[i].HYPER =
            log_sum_log(forward[i - 1].hypo + log(trans[hypo][HYPER]),
                        forward[i - 1].HYPER + log(trans[HYPER][HYPER]),
                        forward[i - 1].HYPO + log(trans[HYPO][HYPER]))
            + HYPER_segment_log_likelihood(i, i + 1);

        // hypomethylated CpG in HyperMR segment
        forward[i].HYPO =
            log_sum_log(forward[i - 1].HYPER + log(trans[HYPER][HYPO]),
                        forward[i - 1].HYPO + log(trans[HYPO][HYPO]))
            + HYPO_segment_log_likelihood(i, i + 1);
    }

    return log_sum_log(forward[end - 1].hypo + lp_end.hypo,
                       forward[end - 1].HYPER + lp_end.HYPER,
                       forward[end - 1].HYPO + lp_end.HYPO);

//     /////
//     cerr << "check forward_algorithm: "<< "OK" << endl;
// /////
}

double
ThreeStateHMM::backward_algorithm(const size_t start, const size_t end)
{
// /////
//     cerr << "check backward_algorithm: "<< "OK" << endl;
// /////

    const int start_int(start), end_int(end);

    for (size_t i = start; i < end; ++i)
        backward[i].hypo = backward[i].HYPER = backward[i].HYPO = 0.0;

    backward[end - 1].hypo = lp_end.hypo;
    backward[end - 1].HYPER = lp_end.HYPER;
    backward[end - 1].HYPO = lp_end.HYPO;

    for (int i = end_int - 2; i >= start_int; --i)
    {
        //  i in hypo-methylated state of HypoMR
        backward[i].hypo =
            log_sum_log(log(trans[hypo][hypo])
                        + hypo_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].hypo,
                        log(trans[hypo][HYPER])
                        + HYPER_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].HYPER);

        //  i in hyper-methylated state of HyperMR
        backward[i].HYPER =
            log_sum_log(log(trans[HYPER][hypo])
                        + hypo_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].hypo,
                        log(trans[HYPER][HYPER])
                        + HYPER_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].HYPER,
                        log(trans[HYPER][HYPO])
                        + HYPO_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].HYPO);

        //  i in hypo-methylated state of HyperMR
        backward[i].HYPO =
            log_sum_log(log(trans[HYPO][HYPER])
                        + HYPER_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].HYPER,
                        log(trans[HYPO][HYPO])
                        + HYPO_segment_log_likelihood(i + 1, i + 2)
                        + backward[i + 1].HYPO);
    }

    return log_sum_log(lp_start.hypo
                       + hypo_segment_log_likelihood(start, start + 1)
                       + backward[start].hypo,
                       lp_start.HYPER
                       + HYPER_segment_log_likelihood(start, start + 1)
                       + backward[start].HYPER,
                       lp_start.HYPO
                       + HYPO_segment_log_likelihood(start, start + 1)
                       + backward[start].HYPO);

// /////
//     cerr << "check backward_algorithm: "<< "OK" << endl;
// /////

}

//////////////////////////////////////////////
//////       Baum-Welch Training        //////
//////////////////////////////////////////////
// Expectation
void
ThreeStateHMM::estimate_state_posterior(const size_t start, const size_t end)
{
// /////
//     cerr << "check enter estimate_state_posterior: "<< "OK" << endl;
// /////

    vector<double> hypo_evidence(end - start, 0), HYPER_evidence(end - start, 0),
        HYPO_evidence(end - start, 0);

    double prev_denom(0), denom(0);
    for (size_t i = start; i < end; ++i)
    {
        hypo_evidence[i - start] = forward[i].hypo + backward[i].hypo;
        HYPER_evidence[i - start] = forward[i].HYPER + backward[i].HYPER;
        HYPO_evidence[i - start] = forward[i].HYPO + backward[i].HYPO;

        denom = log_sum_log(hypo_evidence[i - start],
                            HYPER_evidence[i - start],
                            HYPO_evidence[i - start]);

        if (i > start) assert(fabs(exp(prev_denom - denom) -1) < 1e-6);
        prev_denom = denom;
    }

    for (size_t i = start; i < end; ++i)
    {
        hypo_posteriors[i] = exp(hypo_evidence[i - start] - denom);
        HYPER_posteriors[i] = exp(HYPER_evidence[i - start] - denom);
        HYPO_posteriors[i] = exp(HYPO_evidence[i - start] - denom);

        //renormalize the probabilities
        double sum = hypo_posteriors[i] + HYPER_posteriors[i]
                     + HYPO_posteriors[i];
        hypo_posteriors[i] /= sum;
        HYPER_posteriors[i] /= sum;
        HYPO_posteriors[i] /= sum;
        // if (fabs(hypo_posteriors[i] + HYPER_posteriors[i]
        //          + HYPO_posteriors[i] - 1.0) > 1e-6)
        //     cerr << fabs(hypo_posteriors[i] + HYPER_posteriors[i]
        //                  + HYPO_posteriors[i] - 1.0) << endl;

        assert(fabs(hypo_posteriors[i] + HYPER_posteriors[i]
                    + HYPO_posteriors[i] - 1.0) < 1e-3);
    }

// /////
//     cerr << "check estimate_state_posterior: "<< "OK" << endl;
// /////

}

void
ThreeStateHMM::estimate_posterior_trans_prob(
    const size_t start, const size_t end)
{
    const double denom =
        log_sum_log(forward[start].hypo + backward[start].hypo,
                    forward[start].HYPER + backward[start].HYPER,
                    forward[start].HYPO + backward[start].HYPO);

    for (size_t i = start; i < end - 1; ++i)
    {
        hypo_hypo[i] =
            forward[i].hypo + log(trans[hypo][hypo])
            + hypo_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].hypo - denom;
        hypo_HYPER[i] =
            forward[i].hypo + log(trans[hypo][HYPER])
            + HYPER_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].HYPER - denom;

        HYPER_hypo[i] =
            forward[i].HYPER + log(trans[HYPER][hypo])
            + hypo_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].hypo - denom;
        HYPER_HYPER[i] =
            forward[i].HYPER + log(trans[HYPER][HYPER])
            + HYPER_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].HYPER - denom;
        HYPER_HYPO[i] =
            forward[i].HYPER + log(trans[HYPER][HYPO])
            + HYPO_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].HYPO - denom;

        HYPO_HYPER[i] =
            forward[i].HYPO + log(trans[HYPO][HYPER])
            + HYPER_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].HYPER - denom;
        HYPO_HYPO[i] =
            forward[i].HYPO + log(trans[HYPO][HYPO])
            + HYPO_segment_log_likelihood(i + 1, i + 2)
            + backward[i + 1].HYPO - denom;

        double sum = exp(hypo_hypo[i]) + exp(hypo_HYPER[i]) +
                     exp(HYPER_hypo[i]) + exp(HYPER_HYPER[i]) +
                     exp(HYPER_HYPO[i]) + exp(HYPO_HYPER[i]) +
                     exp(HYPO_HYPO[i]);

        hypo_hypo[i] -= log(sum);
        hypo_HYPER[i] -= log(sum) ;
        HYPER_hypo[i] -= log(sum);
        HYPER_HYPER[i] -= log(sum);
        HYPER_HYPO[i] -= log(sum) ;
        HYPO_HYPER[i] -= log(sum);
        HYPO_HYPO[i] -= log(sum);

        assert(fabs(exp(hypo_hypo[i]) + exp(hypo_HYPER[i])
                    + exp(HYPER_hypo[i]) + exp(HYPER_HYPER[i])
                    + exp(HYPER_HYPO[i]) + exp(HYPO_HYPER[i])
                    + exp(HYPO_HYPO[i]) - 1) < 1e-3);
    }
}

void
ThreeStateHMM::estimate_parameters()
{
    // hypo_emission.fit(meth_lp, unmeth_lp, hypo_posteriors);
    // HYPER_emission.fit(meth_lp, unmeth_lp, HYPER_posteriors);
    // HYPO_emission.fit(meth_lp, unmeth_lp, HYPO_posteriors);

    // assume the emission distribution for short and long HypoMR are the same
    vector<double> combined_hypo_posteriors(hypo_posteriors.size());
    std::transform(hypo_posteriors.begin(), hypo_posteriors.end(),
                   HYPO_posteriors.begin(), combined_hypo_posteriors.begin(),
                   std::plus<double>());
    hypo_emission.fit(meth_lp, unmeth_lp, combined_hypo_posteriors);
    HYPO_emission = hypo_emission;
    HYPER_emission.fit(meth_lp, unmeth_lp, HYPER_posteriors);

    const double sum_hypo_hypo =
        exp(log_sum_log(hypo_hypo.begin(), hypo_hypo.end()));
    const double sum_hypo_HYPER =
        exp(log_sum_log(hypo_HYPER.begin(), hypo_HYPER.end()));
    const double sum_hypo = sum_hypo_hypo + sum_hypo_HYPER;
    trans[hypo][hypo] = sum_hypo_hypo / sum_hypo;
    trans[hypo][HYPER] = sum_hypo_HYPER / sum_hypo;

    const double sum_HYPER_hypo =
        exp(log_sum_log(HYPER_hypo.begin(), HYPER_hypo.end()));
    const double sum_HYPER_HYPER =
        exp(log_sum_log(HYPER_HYPER.begin(), HYPER_HYPER.end()));
    const double sum_HYPER_HYPO =
        exp(log_sum_log(HYPER_HYPO.begin(), HYPER_HYPO.end()));
    const double sum_HYPER = sum_HYPER_hypo + sum_HYPER_HYPER + sum_HYPER_HYPO;
    trans[HYPER][hypo] = sum_HYPER_hypo / sum_HYPER;
    trans[HYPER][HYPER] = sum_HYPER_HYPER / sum_HYPER;
    trans[HYPER][HYPO] = sum_HYPER_HYPO / sum_HYPER;

    const double sum_HYPO_HYPER =
        exp(log_sum_log(HYPO_HYPER.begin(), HYPO_HYPER.end()));
    const double sum_HYPO_HYPO =
        exp(log_sum_log(HYPO_HYPO.begin(), HYPO_HYPO.end()));
    const double sum_HYPO = sum_HYPO_HYPER + sum_HYPO_HYPO;
    trans[HYPO][HYPER] = sum_HYPO_HYPER / sum_HYPO;
    trans[HYPO][HYPO] = sum_HYPO_HYPO / sum_HYPO;

    update_observation_likelihood();
}

double
ThreeStateHMM::single_iteration()
{
    double total_score = 0;

    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double forward_score =
            forward_algorithm(reset_points[i], reset_points[i + 1]);
        const double backward_score =
            backward_algorithm(reset_points[i], reset_points[i + 1]);

        assert(fabs((forward_score - backward_score)
                    / max(forward_score, backward_score))
                    < 1e-10);
        estimate_state_posterior(reset_points[i], reset_points[i + 1]);
        estimate_posterior_trans_prob(reset_points[i], reset_points[i + 1]);
        total_score += forward_score;
    }

    estimate_parameters();
    return total_score;
}

double
ThreeStateHMM::BaumWelchTraining()
{
// /////
//     cerr << "check enter BaumWelchTraining: "<< "OK" << endl;
// /////

    if (VERBOSE)
        cerr << "Baum-Welch Training" << endl;

    double prev_total = -std::numeric_limits<double>::max();

    for (size_t i = 0; i < max_iterations; ++i)
    {
        const betabin old_hypo_emission = hypo_emission;
        const betabin old_HYPER_emission = HYPER_emission;
        const betabin old_HYPO_emission = HYPO_emission;
        const vector<vector<double> > old_trans(trans);

        double total = single_iteration();

        if (VERBOSE)
        {
            cerr << "Itration: " << setw(2) << i + 1 << ";\t"
                 << "Log-Likelihood: " << total << ";\t"
                 << "Delta: " << (total - prev_total)/std::fabs(total)
                 << endl
                 << "hypo: " << old_hypo_emission.tostring() << ";\t"
                 << "HYPER: " << old_HYPER_emission.tostring() << ";\t"
                 << "HYPO: " << old_HYPO_emission.tostring() << endl;

            cerr << setw(5) << "" << setw(10) << "hypo"
                 << setw(10) << "HYPER" << setw(10) << "HYPO" << endl;
            for (size_t r = 0; r < 3; ++r)
            {
                switch (r)
                {
                case 0: cerr << setw(5) << "hypo"; break;
                case 1: cerr << setw(5) << "HYPER"; break;
                case 2: cerr << setw(5) << "HYPO"; break;
                }

                for (size_t c = 0; c < 3; ++c)
                    cerr << setw(10) << old_trans[r][c];
                cerr << endl;
            }
            cerr << endl;
        }


        if ((total - prev_total)/std::fabs(total) < tolerance)
        {
            hypo_emission = old_hypo_emission;
            HYPER_emission = old_HYPER_emission;
            HYPO_emission = old_HYPO_emission;
            update_observation_likelihood();
            trans = old_trans;

            if (VERBOSE)
                cerr << "CONVERGED" << endl << endl;
            break;
        }
        prev_total = total;
    }

// /////
//     cerr << "check BaumWelchTraining: "<< "OK" << endl;
// /////

    return prev_total;
}


//////////////////////////////////////////////
//////          export result           //////
//////////////////////////////////////////////
double
ThreeStateHMM::PosteriorDecoding()
{
    double total_score = 0;

    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double forward_score =
            forward_algorithm(reset_points[i], reset_points[i + 1]);
        const double backward_score =
            backward_algorithm(reset_points[i], reset_points[i + 1]);

        assert(fabs((forward_score - backward_score)
                    / max(forward_score, backward_score))
                    < 1e-10);
        estimate_state_posterior(reset_points[i], reset_points[i + 1]);
        estimate_posterior_trans_prob(reset_points[i], reset_points[i + 1]);
        total_score = log_sum_log(total_score, forward_score);
    }


    for (size_t i = 0; i < observations.size(); ++i)
    {
        state_posteriors[i].hypo = hypo_posteriors[i];
        state_posteriors[i].HYPER = HYPER_posteriors[i];
        state_posteriors[i].HYPO = HYPO_posteriors[i];
        classes[i] = max_state(state_posteriors[i]);
    }

    return total_score;
}


double
ThreeStateHMM::ViterbiDecoding(const size_t start, const size_t end)
{
    if (start >= end)
      throw std::runtime_error("Invalid HMM sequence indices");

    const size_t lim = end - start;

    vector<Triplet> llh(lim);
    vector<STATE_LABELS> trace_hypo(lim);
    vector<STATE_LABELS> trace_HYPER(lim);
    vector<STATE_LABELS> trace_HYPO(lim);

    llh.front().hypo =
        lp_start.hypo + hypo_segment_log_likelihood(start, start + 1);
    llh.front().HYPER =
        lp_start.HYPER + HYPER_segment_log_likelihood(start, start + 1);
    llh.front().HYPO =
        lp_start.HYPO + HYPO_segment_log_likelihood(start, start + 1);

    for (size_t i = 1; i < lim; ++i)
    {
        // hypo:
        const double hypo_hypo = llh[i - 1].hypo + log(trans[hypo][hypo]);
        const double HYPER_hypo = llh[i - 1].HYPER + log(trans[HYPER][hypo]);
        if (hypo_hypo > HYPER_hypo)
        {
            llh[i].hypo = hypo_hypo
                + hypo_segment_log_likelihood(start + i, start + i + 1);
            trace_hypo[i] = hypo;
        }
        else
        {
            llh[i].hypo = HYPER_hypo
                + hypo_segment_log_likelihood(start + i, start + i + 1);
            trace_hypo[i] = HYPER;
        }

        // HYPER
        const double hypo_HYPER = llh[i - 1].hypo + log(trans[hypo][HYPER]);
        const double HYPER_HYPER = llh[i - 1].HYPER + log(trans[HYPER][HYPER]);
        const double HYPO_HYPER = llh[i - 1].HYPER + log(trans[HYPO][HYPER]);
        if (hypo_HYPER >= std::max(HYPER_HYPER, HYPO_HYPER))
        {
            llh[i].HYPER = hypo_HYPER
                + HYPER_segment_log_likelihood(start + i, start + i + 1);
            trace_HYPER[i] = hypo;
        }
        else  if (HYPER_HYPER >= std::max(hypo_HYPER, HYPO_HYPER))
        {
            llh[i].HYPER = HYPER_HYPER
                + HYPER_segment_log_likelihood(start + i, start + i + 1);
            trace_HYPER[i] = HYPER;
        }
        else
        {
            llh[i].HYPER = HYPO_HYPER
                + HYPER_segment_log_likelihood(start + i, start + i + 1);
            trace_HYPER[i] = HYPO;
        }

        // HYPO
        const double HYPER_HYPO = llh[i - 1].HYPER + log(trans[HYPER][HYPO]);
        const double HYPO_HYPO = llh[i - 1].HYPO + log(trans[HYPO][HYPO]);
        if (HYPER_HYPO > HYPO_HYPO)
        {
            llh[i].HYPO = HYPER_HYPO
                + HYPO_segment_log_likelihood(start + i, start + i + 1);
            trace_HYPO[i] = HYPER;
        }
        else
        {
            llh[i].HYPO = HYPO_HYPO
                + HYPO_segment_log_likelihood(start + i, start + i + 1);
            trace_HYPO[i] = HYPO;
        }
    }

    vector<STATE_LABELS> inner_ml_classes;

    // do the traceback
    STATE_LABELS curr = max_state(llh.back());
    for (size_t i = 0; i < trace_hypo.size(); ++i)
    {
        inner_ml_classes.push_back(curr);
        switch (curr)
        {
        case hypo:
            curr = trace_hypo[lim - i - 1];
            break;
        case HYPER:
            curr = trace_HYPER[lim - i - 1];
            break;
        case HYPO:; break;
            curr = trace_HYPO[lim - i - 1];
            break;
        }
    }

    reverse(inner_ml_classes.begin(), inner_ml_classes.end());
    std::copy(inner_ml_classes.begin(), inner_ml_classes.end(),
              classes.begin() + start);

    return max_value(llh.back());
}

double
ThreeStateHMM::ViterbiDecoding()
{
    // ml_classes = vector<bool>(values.size());
    double total = 0;
    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const double llh = ViterbiDecoding(reset_points[i], reset_points[i+1]);
        total = log_sum_log(total, llh);
    }
    return total;
}

void
ThreeStateHMM::get_state_posteriors(std::vector<Triplet> & scores) const
{
    scores = state_posteriors;
}

void
ThreeStateHMM::get_classes(std::vector<STATE_LABELS> & ml_classes) const
{
    ml_classes = classes;
}
