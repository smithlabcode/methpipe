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

#include <cassert>

#include "EvaluateBoundaries.hpp"
#include "numerical_utils.hpp"

using std::vector;
using std::pair;
using std::make_pair;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::string;

void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
					const vector<vector<bool> > &classes,
					const vector<vector<double> > &scores,
					vector<vector<GenomicRegion> > &boundaries) const
{
    // Separate the class values for each bin into domains of contiguous
    // bins having the same class
    const size_t n_regions = classes.size();
    vector<vector<Domain> > domains(n_regions);
    for (size_t i = 0; i < n_regions; ++i)
    {
        size_t prev_end = 0;
        for (size_t j = 0; j < classes[i].size(); ++j)
            if (j > 0 && classes[i][j] != classes[i][j - 1])
            {
                domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
                prev_end = j;
            }
        domains[i].push_back(Domain(scores[i], prev_end, classes[i].size(),
                                    classes[i].back()));
    }

    boundaries.resize(n_regions, vector<GenomicRegion>());
    for (size_t i = 0; i < domains.size(); ++i)
    {
        const string chrom(bin_bounds[i].front().get_chrom());
        boundaries[i].push_back(GenomicRegion(chrom,
                                              bin_bounds[i].front().get_start(),
                                              bin_bounds[i].front().get_start() + 1,
                                              "END", 0, '+'));
        size_t offset = 0;
        for (size_t j = 0; j < domains[i].size() - 1; ++j)
        {
            offset += domains[i][j].vals.size();
            const double peak_score = domains[i][j + 1].vals.front();
            const string peak_name("B:" +
                                   smithlab::toa(static_cast<size_t>(peak_score*1000)));
            boundaries[i].push_back(GenomicRegion(chrom,
                                                  bin_bounds[i][offset].get_start(),
                                                  bin_bounds[i][offset].get_start() + 1,
                                                  peak_name, peak_score, '+'));
        }
        boundaries[i].push_back(boundaries[i].back());
    }
}

void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
					const vector<vector<size_t> > &classes,
					const vector<vector<double> > &scores,
					vector<vector<GenomicRegion> > &boundaries) const
{

    // Separate the class values for each bin into domains of contiguous
    // bins having the same class
    const size_t n_regions = classes.size();
    vector<vector<Domain> > domains(n_regions);
    for (size_t i = 0; i < n_regions; ++i)
    {
        size_t prev_end = 0;
        for (size_t j = 0; j < classes[i].size(); ++j)
            if (j > 0 && classes[i][j] != classes[i][j - 1])
            {
                domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
                prev_end = j;
            }
        domains[i].push_back(Domain(scores[i], prev_end, classes[i].size(),
                                    classes[i].back()));
    }

    boundaries.resize(n_regions, vector<GenomicRegion>());
    for (size_t i = 0; i < domains.size(); ++i)
    {
        const string chrom(bin_bounds[i].front().get_chrom());
        boundaries[i].push_back(GenomicRegion(chrom,
                                              bin_bounds[i].front().get_start(),
                                              bin_bounds[i].front().get_start() + 1,
                                              "END", 0, '+'));
        size_t offset = 0;
        for (size_t j = 0; j < domains[i].size() - 1; ++j)
        {
            offset += domains[i][j].vals.size();
            const double peak_score = domains[i][j + 1].vals.front();
            const string peak_name("B:" +
                                   smithlab::toa(static_cast<size_t>(peak_score*1000)));
            boundaries[i].push_back(GenomicRegion(chrom,
                                                  bin_bounds[i][offset].get_start(),
                                                  bin_bounds[i][offset].get_start() + 1,
                                                  peak_name, peak_score, '+'));
        }
        boundaries[i].push_back(boundaries[i].back());
    }
}



BoundEval::BoundEval(const size_t b, const double bw) :
    boundary_size(b), bandwidth(bw)
{}

string
Domain::tostring() const
{
    std::ostringstream ss;
    for (size_t i = 0; i < vals.size(); ++i)
        ss << vals[i] << endl;
    return ss.str();
}

void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
					const vector<vector<bool> > &classes,
					const vector<vector<double> > &scores,
					vector<vector<GenomicRegion> > &boundaries,
					vector<vector<GenomicRegion> > &boundary_peaks,
					vector<vector<size_t> > &boundary_sizes) const
{
    const size_t n_regions = classes.size();
    vector<vector<Domain> > domains(n_regions);
    for (size_t i = 0; i < n_regions; ++i)
    {
        size_t prev_end = 0;
        for (size_t j = 0; j < classes[i].size(); ++j)
            if (j > 0 && classes[i][j] != classes[i][j - 1])
            {
                domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
                prev_end = j;
            }
        domains[i].push_back(Domain(scores[i], prev_end,
                                    classes[i].size(), classes[i].back()));
    }

    boundaries.resize(n_regions, vector<GenomicRegion>());
    boundary_peaks.resize(n_regions, vector<GenomicRegion>());
    boundary_sizes.resize(n_regions, vector<size_t>());
    for (size_t i = 0; i < domains.size(); ++i)
    {
        const string chrom(bin_bounds[i].front().get_chrom());
        boundaries[i].push_back(GenomicRegion(chrom,
                                              bin_bounds[i].front().get_start(),
                                              bin_bounds[i].front().get_start() + 1,
                                              "END", 0, '+'));
        boundary_peaks[i].push_back(GenomicRegion(chrom,
                                                  bin_bounds[i].front().get_start(),
                                                  bin_bounds[i].front().get_start() + 1,
                                                  "END", 0, '+'));
        boundary_sizes[i].push_back(0);

        size_t offset = 0;
        for (size_t j = 0; j < domains[i].size() - 1; ++j)
        {
            double area_under_curve = 0;

            size_t left_index = 0;
            for (size_t k = 0; k < domains[i][j].vals.size(); ++k)
            {
                const size_t index = domains[i][j].vals.size() - 1 - k;
                if (domains[i][j].vals[index] < 0.01)
                {
                    left_index = k;
                    break;
                }
                area_under_curve += domains[i][j].vals[index];
            }

            // TODO: make sure the bin actually at the boundary is always counted,
            // even if it has a low score.
            size_t right_index = 0;
            for (size_t k = 0; k < domains[i][j + 1].vals.size(); ++k)
            {
                if (domains[i][j + 1].vals[k] < 0.01)
                {
                    right_index = k;
                    break;
                }
                area_under_curve += domains[i][j + 1].vals[k];
            }
            // TODO: see above todo
            // area_under_curve = max(0.01, area_under_curve);

            offset += domains[i][j].vals.size();
            const size_t bound_start = bin_bounds[i][offset - left_index].get_start();
            const size_t bound_end = bin_bounds[i][offset + right_index].get_end();
            const size_t bound_bins = left_index + right_index;
            const double peak_score = domains[i][j + 1].vals.front();
            const double bound_score = peak_score -
                (area_under_curve - peak_score)/(max(1, int(bound_bins - 1)));

            const string bound_name("B:" +
                                    smithlab::toa(static_cast<size_t>(bound_score*1000)));
            const string peak_name("B:" +
                                   smithlab::toa(static_cast<size_t>(peak_score*1000)));
            boundaries[i].push_back(GenomicRegion(chrom, bound_start, bound_end,
                                                  bound_name, bound_score, '+'));
            boundary_peaks[i].push_back(GenomicRegion(chrom,
                                                      bin_bounds[i][offset].get_start(),
                                                      bin_bounds[i][offset].get_end(),
                                                      peak_name, peak_score, '+'));
            boundary_sizes[i].push_back(bound_bins);
        }
        boundaries[i].push_back(GenomicRegion(chrom,
                                              bin_bounds[i].back().get_start(),
                                              bin_bounds[i].back().get_end(),
                                              "END", 0, '+'));
        boundary_peaks[i].push_back(boundaries[i].back());
        boundary_sizes[i].push_back(0);
    }
}




void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
					const vector<vector<size_t> > &classes,
					const vector<vector<double> > &scores,
					vector<vector<GenomicRegion> > &boundaries,
					vector<vector<GenomicRegion> > &boundary_peaks,
					vector<vector<size_t> > &boundary_sizes) const
{

    const size_t n_regions = classes.size();
    vector<vector<Domain> > domains(n_regions);
    for (size_t i = 0; i < n_regions; ++i)
    {
        size_t prev_end = 0;
        for (size_t j = 0; j < classes[i].size(); ++j)
            if (j > 0 && classes[i][j] != classes[i][j - 1])
            {
                domains[i].push_back(Domain(scores[i], prev_end, j, classes[i][j - 1]));
                prev_end = j;
            }
        domains[i].push_back(Domain(scores[i], prev_end,
                                    classes[i].size(), classes[i].back()));
    }

    boundaries.resize(n_regions, vector<GenomicRegion>());
    boundary_peaks.resize(n_regions, vector<GenomicRegion>());
    boundary_sizes.resize(n_regions, vector<size_t>());
    for (size_t i = 0; i < domains.size(); ++i)
    {
        const string chrom(bin_bounds[i].front().get_chrom());
        boundaries[i].push_back(GenomicRegion(chrom,
                                              bin_bounds[i].front().get_start(),
                                              bin_bounds[i].front().get_end(),
                                              "END", 0, '+'));
        boundary_peaks[i].push_back(GenomicRegion(chrom,
                                                  bin_bounds[i].front().get_start(),
                                                  bin_bounds[i].front().get_end(),
                                                  "END", 0, '+'));
        boundary_sizes[i].push_back(0);

        size_t offset = 0;
        for (size_t j = 0; j < domains[i].size() - 1; ++j)
        {
            double area_under_curve = 0;

            size_t left_index = 0;
            for (size_t k = 0; k < domains[i][j].vals.size(); ++k)
            {
                const size_t index = domains[i][j].vals.size() - 1 - k;
                if (domains[i][j].vals[index] < 0.01)
                {
                    left_index = k;
                    break;
                }
                area_under_curve += domains[i][j].vals[index];
            }

            // TODO: make sure the bin actually at the boundary is always counted,
            // even if it has a low score.
            size_t right_index = 0;
            for (size_t k = 0; k < domains[i][j + 1].vals.size(); ++k)
            {
                if (domains[i][j + 1].vals[k] < 0.01)
                {
                    right_index = k;
                    break;
                }
                area_under_curve += domains[i][j + 1].vals[k];
            }
            // TODO: see above todo
            // area_under_curve = max(0.01, area_under_curve);

            offset += domains[i][j].vals.size();
            const size_t bound_start = bin_bounds[i][offset - left_index].get_start();
            const size_t bound_end = bin_bounds[i][offset + right_index].get_end();
            const size_t bound_bins = left_index + right_index;
            const double peak_score = domains[i][j + 1].vals.front();
            const double bound_score = peak_score -
                (area_under_curve - peak_score)/(max(1, int(bound_bins - 1)));

            const string bound_name("B:" +
                                    smithlab::toa(static_cast<size_t>(bound_score*1000)));
            const string peak_name("B:" +
                                   smithlab::toa(static_cast<size_t>(peak_score*1000)));
            boundaries[i].push_back(GenomicRegion(chrom, bound_start, bound_end,
                                                  bound_name, bound_score, '+'));
            boundary_peaks[i].push_back(GenomicRegion(chrom,
                                                      bin_bounds[i][offset].get_start(),
                                                      bin_bounds[i][offset].get_end(),
                                                      peak_name, peak_score, '+'));
            boundary_sizes[i].push_back(bound_bins);
        }
        boundaries[i].push_back(GenomicRegion(chrom,
                                              bin_bounds[i].back().get_start(),
                                              bin_bounds[i].back().get_end(),
                                              "END", 0, '+'));
        boundary_peaks[i].push_back(boundaries[i].back());
        boundary_sizes[i].push_back(0);
    }
}

double
bound_trans_prob(const vector<double> &ff_probs,
                 const vector<double> &fb_probs,
                 const vector<double> &bf_probs,
                 const vector<double> &bb_probs,
                 const size_t start,
                 const size_t end)
{
//     static size_t c = 0;

    const size_t sz = end - start;

    vector<double> log_ff_probs(ff_probs.begin() + start, ff_probs.begin() + end);
    vector<double> log_fb_probs(fb_probs.begin() + start, fb_probs.begin() + end);
    vector<double> log_bf_probs(bf_probs.begin() + start, bf_probs.begin() + end);
    vector<double> log_bb_probs(bb_probs.begin() + start, bb_probs.begin() + end);


    for (size_t i = 0; i < sz; ++i)
    {
        log_ff_probs[i] = log(log_ff_probs[i]);
        log_fb_probs[i] = log(log_fb_probs[i]);
        log_bf_probs[i] = log(log_bf_probs[i]);
        log_bb_probs[i] = log(log_bb_probs[i]);
    }

    double sum = 0;

    // transition occurs in the start
    double prod = log_fb_probs[0];
    for (size_t i = 1; i < sz; ++i)
        prod += log_bb_probs[i] - log_sum_log(log_bb_probs[i], log_bf_probs[i]);

    sum = prod;

    // transition occurs afterwards
    prod -= log_fb_probs[0];
    prod += log_ff_probs[0];
    for (size_t i = 1; i < sz; ++i)
    {
        if (i > 1)
        {
            prod -= log_fb_probs[i-1] - log_sum_log(log_ff_probs[i-1], log_fb_probs[i-1]);
            prod += log_ff_probs[i-1] - log_sum_log(log_ff_probs[i-1], log_fb_probs[i-1]);
        }
        prod -= log_bb_probs[i] - log_sum_log(log_bb_probs[i], log_bf_probs[i]);
        prod += log_fb_probs[i] - log_sum_log(log_ff_probs[i], log_fb_probs[i]);
        sum = log_sum_log(sum, prod);
    }

//     cerr << "check" << start << "\t" << end << "\t" << prod << "\t" << sum << endl;
//     ++c;
//     if (c > 10) exit(0);

    return exp(sum);
}

void
make_boundary(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
              const vector<double>  &trans_scores,
              const vector<double> &fg_to_fg_trans_score,
              const vector<double> &fg_to_bg_trans_score,
              const vector<double> &bg_to_fg_trans_score,
              const vector<double> &bg_to_bg_trans_score,
              const size_t i,
              const size_t offset,
              const size_t left_index,
              const size_t right_index,
              GenomicRegion &bound)
{
    const string bound_chrom =
        bin_bounds[i][left_index - offset].get_chrom();
    const size_t bound_start =
        bin_bounds[i][left_index - offset].get_start();
    const size_t bound_end =
        bin_bounds[i][right_index - offset - 1].get_end();

    const size_t peak_index =
        std::max_element(trans_scores.begin() + left_index,
                         trans_scores.begin() + right_index)
        - trans_scores.begin();
    const double peak_score = trans_scores[peak_index];
    const size_t peak_loc = bin_bounds[i][peak_index - offset].get_start();

    const double bound_score = std::max(
        bound_trans_prob(bg_to_bg_trans_score,
                         bg_to_fg_trans_score,
                         fg_to_bg_trans_score,
                         fg_to_fg_trans_score,
                         left_index,
                         right_index),
        bound_trans_prob(fg_to_fg_trans_score,
                         fg_to_bg_trans_score,
                         bg_to_fg_trans_score,
                         bg_to_bg_trans_score,
                         left_index,
                         right_index));
    const string bound_name =
        "B:" + smithlab::toa(right_index - left_index) + ":" +
        smithlab::toa(peak_loc) + ":" +
        smithlab::toa(peak_score);

    bound.set_chrom(bound_chrom);
    bound.set_start(bound_start);
    bound.set_end(bound_end);
    bound.set_name(bound_name);
    bound.set_score(bound_score);
    bound.set_strand('+');
}

void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
                    const vector<size_t> &reset_points,
                    const vector<bool> &classes,
                    const vector<double>  &trans_scores,
                    const vector<double> &fg_to_fg_trans_score,
                    const vector<double> &fg_to_bg_trans_score,
                    const vector<double> &bg_to_fg_trans_score,
                    const vector<double> &bg_to_bg_trans_score,
                    const double cutoff,
                    const bool Both_Domain_Ends,
                    vector<GenomicRegion> &boundaries) const
{
    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {

        const size_t offset = reset_points[i];
        const size_t start = reset_points[i];
        const size_t end = reset_points[i + 1];

        // find domains and domain ends
        vector<pair<size_t, size_t> > domain_ends;
        bool prev_class = false;
        if (classes[start])
        {
            prev_class = true;
            domain_ends.push_back(make_pair(start, 0));
        }
        for (size_t j = start + 1; j < end; ++j)
            if (classes[j] != classes[j-1])
            {
                if (prev_class)
                {
                    domain_ends.back().second = j;
                    prev_class = false;
                }
                else
                {
                    domain_ends.push_back(make_pair(j, 0));
                    prev_class = true;
                }
            }
        if (prev_class)
            domain_ends.back().second = end;

        // build and evaluate boundaries
        for (size_t j = 0; j < domain_ends.size(); ++j)
        {
            size_t first_left_index = domain_ends[j].first;
            while (first_left_index > start &&
                   trans_scores[first_left_index - 1] > cutoff)
                --first_left_index;

            size_t first_right_index = domain_ends[j].first;
            while (first_right_index < end &&
                   trans_scores[first_right_index] > cutoff)
                ++first_right_index;

            size_t second_left_index = domain_ends[j].second;
            while (second_left_index > start &&
                   trans_scores[second_left_index - 1] > cutoff)
                --second_left_index;

            size_t second_right_index = domain_ends[j].second;
            while (second_right_index < end &&
                   trans_scores[second_right_index] > cutoff)
                ++second_right_index;

            if (Both_Domain_Ends &&
                (first_left_index == first_right_index ||
                 second_left_index == second_right_index))
                continue;

            // output first bound
            if (first_left_index < first_right_index)
            {
                GenomicRegion bound;
                make_boundary(bin_bounds, trans_scores,
                              fg_to_fg_trans_score, fg_to_bg_trans_score,
                              bg_to_fg_trans_score, bg_to_bg_trans_score,
                              i, offset, first_left_index, first_right_index,
                              bound);
                boundaries.push_back(bound);
            }

            // output second bound
            if (second_left_index < second_right_index)
            {
                GenomicRegion bound;
                make_boundary(bin_bounds, trans_scores,
                              fg_to_fg_trans_score, fg_to_bg_trans_score,
                              bg_to_fg_trans_score, bg_to_bg_trans_score,
                              i, offset, second_left_index, second_right_index,
                              bound);
                boundaries.push_back(bound);
            }
        }
    }
}


void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
                    const vector<size_t> &reset_points,
                    const vector<bool> &classes,
                    const vector<double> &trans_scores,
                    const vector<double> &fg_to_fg_trans_score,
                    const vector<double> &fg_to_bg_trans_score,
                    const vector<double> &bg_to_fg_trans_score,
                    const vector<double> &bg_to_bg_trans_score,
                    const double cutoff,
                    vector<GenomicRegion> &boundaries) const
{
    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const size_t offset = reset_points[i];
        const size_t start = reset_points[i];
        const size_t end = reset_points[i + 1];

        size_t left_index = start;
        size_t right_index = left_index;
        while (left_index < end)
        {
            while (left_index < end && trans_scores[left_index] < cutoff)
                ++left_index;

            right_index = left_index;
            while (right_index < end && trans_scores[right_index] >= cutoff)
                ++right_index;
            if (left_index < end)
            {
                const string bound_chrom =
                    bin_bounds[i][left_index - offset].get_chrom();
                const size_t bound_start =
                    bin_bounds[i][left_index - offset].get_start();
                const size_t bound_end =
                    bin_bounds[i][right_index - offset - 1].get_end();

                const size_t peak_index =
                    std::max_element(trans_scores.begin() + left_index,
                                      trans_scores.begin() + right_index)
                    - trans_scores.begin();
                const double peak_score = trans_scores[peak_index];
                const size_t peak_loc = bin_bounds[i][peak_index - offset].get_start();

                const double bound_score = std::max(
                    bound_trans_prob(bg_to_bg_trans_score,
                                     bg_to_fg_trans_score,
                                     fg_to_bg_trans_score,
                                     fg_to_fg_trans_score,
                                     left_index,
                                     right_index),
                    bound_trans_prob(fg_to_fg_trans_score,
                                     fg_to_bg_trans_score,
                                     bg_to_fg_trans_score,
                                     bg_to_bg_trans_score,
                                     left_index,
                                     right_index));
                const string bound_name =
                    "B:" + smithlab::toa(right_index - left_index) + ":" +
                    smithlab::toa(peak_loc) + ":" +
                    smithlab::toa(peak_score);
                boundaries.push_back(GenomicRegion(bound_chrom, bound_start, bound_end,
                                                   bound_name, bound_score, '+'));
            }
            left_index = right_index;
        }
    }
}


double
bound_trans_prob(const vector<vector<vector<double> > > &post_trans,
                 const size_t from_state, const size_t to_state,
                 const size_t start, const size_t end)
{
    const size_t other_state = 3 - from_state - to_state;

    double sum = 0;

    // transition occurs in the start
    double prod = post_trans[from_state][to_state][start];
    for (size_t i = start + 1; i < end; ++i)
        prod *= post_trans[to_state][to_state][i] /
            (post_trans[to_state][from_state][i] +
             post_trans[to_state][to_state][i] +
             post_trans[to_state][other_state][i]);

    sum = prod;

    // transition occurs afterwards
    prod /= post_trans[from_state][to_state][start];
    prod *= post_trans[from_state][from_state][start];
    for (size_t i = start + 1; i < end; ++i)
    {
        if (i > start + 1)
        {
            prod /= post_trans[from_state][to_state][i-1] /
            (post_trans[from_state][from_state][i-1] +
             post_trans[from_state][to_state][i-1] +
             post_trans[from_state][other_state][i-1]);

            prod *= post_trans[from_state][from_state][i-1] /
            (post_trans[from_state][from_state][i-1] +
             post_trans[from_state][to_state][i-1] +
             post_trans[from_state][other_state][i-1]);
        }

        prod /= post_trans[to_state][to_state][i] /
            (post_trans[to_state][from_state][i] +
             post_trans[to_state][to_state][i] +
             post_trans[to_state][other_state][i]);

        prod *= post_trans[from_state][to_state][i] /
            (post_trans[from_state][from_state][i] +
             post_trans[from_state][to_state][i] +
             post_trans[from_state][other_state][i]);

        sum += prod;
    }

    return sum;
}

// void
// BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
//                     const vector<size_t> &reset_points,
//                     const vector<size_t> &classes,
//                     const vector<double> &trans_scores,
//                     const vector<vector<vector<double> > > &post_trans,
//                     const double cutoff,
//                     vector<GenomicRegion> &boundaries) const
// {
//     for (size_t i = 0; i < reset_points.size() - 1; ++i)
//     {
//         const size_t offset = reset_points[i];
//         const size_t start = reset_points[i];
//         const size_t end = reset_points[i + 1];

//         size_t left_index = start;
//         size_t right_index = left_index;
//         while (left_index < end)
//         {
//             while (left_index < end && trans_scores[left_index] < cutoff)
//                 ++left_index;

//             right_index = left_index;
//             while (right_index < end && trans_scores[right_index] >= cutoff)
//                 ++right_index;
//             if (left_index < end)
//             {
//                 const string bound_chrom =
//                     bin_bounds[i][left_index - offset].get_chrom();
//                 const size_t bound_start =
//                     bin_bounds[i][left_index - offset].get_start();
//                 const size_t bound_end =
//                     bin_bounds[i][right_index - offset - 1].get_end();

//                 const size_t peak_index =
//                     std::max_element(trans_scores.begin() + left_index,
//                                       trans_scores.begin() + right_index)
//                     - trans_scores.begin();
//                 const double peak_score = trans_scores[peak_index];
//                 const size_t peak_loc = bin_bounds[i][peak_index - offset].get_start();

//                 const size_t state_a = classes[left_index];
//                 size_t state_b, state_c;
//                 if (state_a == 0)
//                 {
//                     state_b = 1;
//                     state_c = 2;
//                 }
//                 else if (state_a == 1)
//                 {
//                     state_b = 0;
//                     state_c = 2;
//                 }
//                 else // if (state_a == 2)
//                 {
//                     state_b = 0;
//                     state_c = 1;
//                 }

//                 const double bound_score = std::max(
//                     std::max(bound_trans_prob(post_trans,
//                                               state_a, state_b,
//                                               left_index, right_index),
//                              bound_trans_prob(post_trans,
//                                               state_b, state_a,
//                                               left_index, right_index)),
//                     std::max(bound_trans_prob(post_trans,
//                                               state_a, state_c,
//                                               left_index, right_index),
//                              bound_trans_prob(post_trans,
//                                               state_c, state_a,
//                                               left_index, right_index)));
//                 const string bound_name =
//                     "B:" + smithlab::toa(right_index - left_index) + ":" +
//                     smithlab::toa(peak_loc) + ":" +
//                     smithlab::toa(peak_score);
//                 boundaries.push_back(GenomicRegion(bound_chrom, bound_start, bound_end,
//                                                    bound_name, bound_score, '+'));
//             }
//             left_index = right_index;
//         }
//     }
// }

void
BoundEval::evaluate(const vector<vector<SimpleGenomicRegion> > &bin_bounds,
                    const vector<size_t> &reset_points,
                    const vector<size_t> &classes,
                    const vector<double> &trans_scores,
                    const vector<vector<vector<double> > > &post_trans,
                    const double cutoff,
                    vector<GenomicRegion> &boundaries) const
{
    for (size_t i = 0; i < reset_points.size() - 1; ++i)
    {
        const size_t offset = reset_points[i];
        const size_t start = reset_points[i];
        const size_t end = reset_points[i + 1];

        // find domains and domain ends
        vector<size_t> domain_ends;
        domain_ends.push_back(start);
        for (size_t j = start + 1; j < end; ++j)
            if (classes[j] != classes[j-1])
                domain_ends.push_back(j);

        domain_ends.push_back(end);

        for (size_t j = 0; j < domain_ends.size(); ++j)
        {
            size_t left_index = domain_ends[j];
            while (left_index > start &&
                   trans_scores[left_index - 1] > cutoff)
                --left_index;

            size_t right_index = domain_ends[j];
            while (right_index < end &&
                   trans_scores[right_index] > cutoff)
                ++right_index;

            if (left_index < right_index)
            {
                const string bound_chrom =
                    bin_bounds[i][left_index - offset].get_chrom();
                const size_t bound_start =
                    bin_bounds[i][left_index - offset].get_start();
                const size_t bound_end =
                    bin_bounds[i][right_index - offset - 1].get_end();

                const size_t peak_index =
                    std::max_element(trans_scores.begin() + left_index,
                                     trans_scores.begin() + right_index)
                    - trans_scores.begin();
                const double peak_score = trans_scores[peak_index];
                const size_t peak_loc = bin_bounds[i][peak_index - offset].get_start();

                const size_t state_a = classes[left_index];
                size_t state_b, state_c;
                if (state_a == 0)
                {
                    state_b = 1;
                    state_c = 2;
                }
                else if (state_a == 1)
                {
                    state_b = 0;
                    state_c = 2;
                }
                else // if (state_a == 2)
                {
                    state_b = 0;
                    state_c = 1;
                }

                const double bound_score = std::max(
                    std::max(bound_trans_prob(post_trans,
                                              state_a, state_b,
                                              left_index, right_index),
                             bound_trans_prob(post_trans,
                                              state_b, state_a,
                                              left_index, right_index)),
                    std::max(bound_trans_prob(post_trans,
                                              state_a, state_c,
                                              left_index, right_index),
                             bound_trans_prob(post_trans,
                                              state_c, state_a,
                                              left_index, right_index)));
                const string bound_name =
                    "B:" + smithlab::toa(right_index - left_index) + ":" +
                    smithlab::toa(peak_loc) + ":" +
                    smithlab::toa(peak_score);
                boundaries.push_back(GenomicRegion(bound_chrom, bound_start, bound_end,
                                                   bound_name, bound_score, '+'));
            }
        }
    }
}



