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

#include <iomanip>
#include <cmath>
#include <cassert>
#include <numeric>
#include <limits>
#include <utility>
#include <algorithm>
#include <functional>

#include "TwoStateScaleResolveMixture.hpp"
#include "numerical_utils.hpp"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::min;
using std::pair;
using std::make_pair;
using std::isfinite;

static void
initialize_values(const vector<double> &values,
                  const vector<double> &scales,
                  Distro &fg_distro,
                  Distro &bg_distro) {

    vector<double> fg_part, bg_part, fg_scales, bg_scales;

    const double val_weight = accumulate(values.begin(), values.end(), 0.0);
    assert(val_weight > 0);

    const double val_per_part = std::ceil(val_weight/2);

    vector<pair<double, size_t> > helper;
    for (size_t i = 0; i < values.size(); ++i)
        helper.push_back(make_pair(values[i], i));
    sort(helper.begin(), helper.end(), std::greater<pair<double, size_t> >());

    double sum = 0;
    size_t i = 0;
    for (; i < helper.size() && sum < val_per_part; ++i) {
        fg_part.push_back(values[helper[i].second]);
        fg_scales.push_back(scales[helper[i].second]);
        sum += fabs(helper[i].first);
    }
    for (; i < helper.size(); ++i) {
        bg_part.push_back(values[helper[i].second]);
        bg_scales.push_back(scales[helper[i].second]);
    }
    fg_distro.estimate_params_ml(fg_part, fg_scales);
    bg_distro.estimate_params_ml(bg_part, bg_scales);
}

static double
expectation_step(const vector<double> &values, const vector<double> &scales,
                 const double &mixing,
                                 const Distro &fg_distro, const Distro &bg_distro,
                                 vector<double> &fg_probs, vector<double> &bg_probs)
{
    double score = 0;

    const double fg_log_mixing = log(mixing);
    assert(isfinite(fg_log_mixing));
    const double bg_log_mixing = log(1 - mixing);
    assert(isfinite(bg_log_mixing));

    for (size_t i = 0; i < values.size(); ++i)
    {

        const double fg_part = fg_log_mixing
            + fg_distro.log_likelihood(values[i], scales[i]);
        assert(isfinite(fg_part));

        const double bg_part = bg_log_mixing
            + bg_distro.log_likelihood(values[i], scales[i]);
        assert(isfinite(fg_part));

        const double denom = ((fg_part > bg_part) ?
                              fg_part + log(1.0 + exp(bg_part - fg_part)) :
                              bg_part + log(1.0 + exp(fg_part - bg_part)));
        assert(isfinite(denom));

        fg_probs[i] = exp(fg_part - denom);
        bg_probs[i] = exp(bg_part - denom);

        score += denom;
    }
    return score;
}

static void
maximization_step(const vector<double> &values,
                  const vector<double> &scales,
                                  const vector<double> &fg_probs,
                                  const vector<double> &bg_probs,
                                  double &mixing,
                                  Distro &fg_distro, Distro &bg_distro)
{
    fg_distro.estimate_params_ml(values, scales, fg_probs);
    bg_distro.estimate_params_ml(values, scales, bg_probs);

    vector<double> log_fg_probs(fg_probs.size());
    vector<double> log_bg_probs(bg_probs.size());
    for (size_t i = 0; i < log_fg_probs.size(); ++i)
    {
        log_fg_probs[i] = log(fg_probs[i]);
        log_bg_probs[i] = log(bg_probs[i]);
    }

    mixing = log_sum_log_vec(log_fg_probs, log_fg_probs.size());
    const double bg_mixing = log_sum_log_vec(log_bg_probs, log_bg_probs.size());
    const double mix_sum = ((mixing > bg_mixing) ?
                            mixing + log(1 + exp(bg_mixing - mixing)) :
                            bg_mixing + log(1 + exp(mixing - bg_mixing)));
    mixing = exp(mixing - mix_sum);
    // assert(false && "This code should not be used");
}

void
TwoStateResolveMixture(
    const std::vector<double> &values,
    const std::vector<double> &scales,
    const size_t max_iterations,
    const double tolerance, int VERBOSE,
    Distro &fg_distro, Distro &bg_distro, double &mixing)
{
    // partition the observations to get the initial guess
    initialize_values(values, scales, fg_distro, bg_distro);

    vector<double> fg_probs(values.size(), 0), bg_probs(values.size(), 0);
    mixing = 0.5;

    if (VERBOSE)
        cout << endl << std::setw(10) << "DELTA"
             << std::setw(14) << "(PARAMS,MIX)" << endl;

    // Do the expectation maximization
    double prev_score = std::numeric_limits<double>::max();
    for (size_t itr = 0; itr < max_iterations; ++itr)
    {
        const double score =
            expectation_step(values, scales, mixing,
                             fg_distro, bg_distro, fg_probs, bg_probs);
        maximization_step(values, scales, fg_probs, bg_probs, mixing,
                          fg_distro, bg_distro);
        if (VERBOSE)
        {
            cout << std::setw(10) << std::setprecision(4)
                 << (prev_score - score)/prev_score << "\t"
                 << std::setw(14) << fg_distro.tostring() << " "
                 << std::setw(10) << mixing << " "
                 << std::setw(14) << bg_distro.tostring() << " "
                 << std::setw(10) << 1 - mixing << " "
                 << endl;
        }
        if ((prev_score - score)/prev_score < tolerance)
            break;
        prev_score = score;
    }
}
