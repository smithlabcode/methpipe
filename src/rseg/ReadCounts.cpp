/*
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 *
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

#include <cmath>
#include <cassert>
#include <numeric>
#include <sys/types.h>
#include <unistd.h>

#include "ReadCounts.hpp"

using std::vector;
using std::string;
using std::min;
using std::max;

using std::cerr;
using std::endl;

void
AdjustBinSize(const vector<double> &old_read_bins,
              const vector<double> &old_nondead_scales,
              const vector<size_t> &old_reset_points,
              const size_t old_bin_size,
              vector<double>  &read_bins,
              vector<double>  &nondead_scales,
              vector<size_t> &reset_points,
              const size_t bin_size)
{
    assert(bin_size % old_bin_size == 0);
    read_bins.clear();
    nondead_scales.clear();
    reset_points.clear();

    const size_t n_steps = bin_size / old_bin_size;
    reset_points.push_back(0);
    for (size_t i = 0; i < old_reset_points.size() - 1; ++i)
    {
        const size_t start = old_reset_points[i];
        const size_t end = old_reset_points[i + 1];
        size_t j = start;
        while (j +  n_steps <= end)
        {
            read_bins.push_back(
                std::accumulate(old_read_bins.begin() + j,
                                old_read_bins.begin() + j + n_steps, 0.0));

            nondead_scales.push_back(
                std::accumulate(old_nondead_scales.begin() + j,
                                old_nondead_scales.begin() + j + n_steps,
                                0.0) / n_steps);
            j += n_steps;
        }
        reset_points.push_back(read_bins.size());
    }
    assert(*std::max_element(nondead_scales.begin(), nondead_scales.end())<=1.0);
    assert(*std::min_element(nondead_scales.begin(), nondead_scales.end())>=0.0);
}

static double runif() {
  return static_cast<double>(rand())/RAND_MAX;
}

void
GetCorrectedReadCounts(const vector<double> &read_bins,
                       const vector<double> &nondead_scales,
                       vector<double> &corrected_read_bins,
                       const double max_dead_proportion)
{
    assert(read_bins.size() == nondead_scales.size());
    corrected_read_bins.clear();

    for (size_t i = 0; i < read_bins.size(); ++i)
        if (nondead_scales[i] > 1 - max_dead_proportion)
        {
            const double corrected = read_bins[i] / nondead_scales[i];
            const double floor_count = std::floor(corrected);
            const double frac_part = corrected - floor_count;
            corrected_read_bins.push_back(
                (1.1*runif() < frac_part) ?
                std::ceil(corrected) : floor_count);
        }
}


void
AdjustBinSize(vector<SimpleGenomicRegion> &old_bin_boundaries,
              vector<double> &old_read_bins,
              vector<double> &old_nondead_scales,
              vector<size_t> &old_reset_points,
              const size_t old_bin_size,
              const size_t bin_size)
{
    assert(bin_size % old_bin_size == 0);

    vector<SimpleGenomicRegion> bin_boundaries;
    vector<double> read_bins;
    vector<double> nondead_scales;
    vector<size_t> reset_points;

    const size_t n_steps = bin_size / old_bin_size;

    reset_points.push_back(0);
    for (size_t i = 0; i < old_reset_points.size() - 1; ++i)
    {
        const size_t start = old_reset_points[i];
        const size_t end = old_reset_points[i + 1];
        size_t j = start;
        while (j + n_steps <= end )
        {
            assert(j + n_steps - 1 < old_bin_boundaries.size());
            bin_boundaries.push_back(old_bin_boundaries[j]);
            bin_boundaries.back().set_end(
                old_bin_boundaries[j + n_steps - 1].get_end());

            read_bins.push_back(
                std::accumulate(old_read_bins.begin() + j,
                                old_read_bins.begin() + j + n_steps, 0.0));

            nondead_scales.push_back(
                std::accumulate(old_nondead_scales.begin() + j,
                                old_nondead_scales.begin() + j + n_steps,
                                0.0) / n_steps);
            j += n_steps;
        }
        reset_points.push_back(bin_boundaries.size());
    }

    assert(bin_boundaries.size() == reset_points.back());
    assert(read_bins.size() == reset_points.back());
    assert(nondead_scales.size() == reset_points.back());

    assert(*std::max_element(nondead_scales.begin(), nondead_scales.end())<=1.0);
    assert(*std::min_element(nondead_scales.begin(), nondead_scales.end())>=0.0);

    std::swap(old_bin_boundaries, bin_boundaries);
    std::swap(old_read_bins, read_bins);
    std::swap(old_nondead_scales, nondead_scales);
    std::swap(old_reset_points, reset_points);

}

void
RemoveDeserts(vector<SimpleGenomicRegion> &old_bin_boundaries,
              vector<double> &old_read_bins,
              vector<double> &old_nondead_scales,
              vector<size_t> &old_reset_points,
              const size_t bin_size,
              const size_t desert_size,
              const double max_dead_proportion)
{
    vector<SimpleGenomicRegion> bin_boundaries;
    vector<double> read_bins;
    vector<double> nondead_scales;
    vector<size_t> reset_points;

    reset_points.push_back(0);
    for (size_t i = 0; i < old_reset_points.size() - 1; ++i)
    {
        const size_t start = old_reset_points[i];
        const size_t end = old_reset_points[i + 1];
        size_t n_gaps = 0;
        for (size_t j = start; j < end; ++j)
            if (old_nondead_scales[j] > 1 - max_dead_proportion)
            {
                if (n_gaps * bin_size > desert_size
                    && reset_points.back() < bin_boundaries.size())
                    reset_points.push_back(bin_boundaries.size());

                bin_boundaries.push_back(old_bin_boundaries[j]);
                read_bins.push_back(old_read_bins[j]);
                nondead_scales.push_back(old_nondead_scales[j]);

                n_gaps = 0;
            }
            else
            {
                ++n_gaps;
            }
        if (reset_points.back() < bin_boundaries.size())
            reset_points.push_back(bin_boundaries.size());
    }


    std::swap(old_bin_boundaries, bin_boundaries);
    std::swap(old_read_bins, read_bins);
    std::swap(old_nondead_scales, nondead_scales);
    std::swap(old_reset_points, reset_points);
}

void
AdjustBinSize(vector<SimpleGenomicRegion> &old_bin_boundaries,
              vector<double> &old_read_bins_a,
              vector<double> &old_read_bins_b,
              vector<double> &old_nondead_scales,
              vector<size_t> &old_reset_points,
              const size_t old_bin_size,
              const size_t bin_size)
{
    assert(bin_size % old_bin_size == 0);
    vector<SimpleGenomicRegion> bin_boundaries;
    vector<double> read_bins_a;
    vector<double> read_bins_b;
    vector<double> nondead_scales;
    vector<size_t> reset_points;

    const size_t n_steps = bin_size / old_bin_size;
    reset_points.push_back(0);
    for (size_t i = 0; i < old_reset_points.size() - 1; ++i)
    {
        const size_t start = old_reset_points[i];
        const size_t end = old_reset_points[i + 1];
        size_t j = start;
        while (j + n_steps <= end )
        {
            bin_boundaries.push_back(old_bin_boundaries[j]);
            bin_boundaries.back().set_end(
                old_bin_boundaries[j + n_steps - 1].get_end());

            read_bins_a.push_back(
                std::accumulate(old_read_bins_a.begin() + j,
                                old_read_bins_a.begin() + j + n_steps,
                                0.0));

            read_bins_b.push_back(
                std::accumulate(old_read_bins_b.begin() + j,
                                old_read_bins_b.begin() + j + n_steps,
                                0.0));

            nondead_scales.push_back(
                std::accumulate(old_nondead_scales.begin() + j,
                                old_nondead_scales.begin() + j + n_steps,
                                0.0) / n_steps);
            j += n_steps;
        }
        reset_points.push_back(bin_boundaries.size());
    }

    assert(bin_boundaries.size() == reset_points.back());
    assert(read_bins_a.size() == reset_points.back());
    assert(read_bins_b.size() == reset_points.back());
    assert(nondead_scales.size() == reset_points.back());

    assert(*std::max_element(nondead_scales.begin(), nondead_scales.end())<=1.0);
    assert(*std::min_element(nondead_scales.begin(), nondead_scales.end())>=0.0);

    std::swap(old_bin_boundaries, bin_boundaries);
    std::swap(old_read_bins_a, read_bins_a);
    std::swap(old_read_bins_b, read_bins_b);
    std::swap(old_nondead_scales, nondead_scales);
    std::swap(old_reset_points, reset_points);
}

void
RemoveDeserts(vector<SimpleGenomicRegion> &old_bin_boundaries,
              vector<double> &old_read_bins_a,
              vector<double> &old_read_bins_b,
              vector<double> &old_nondead_scales,
              vector<size_t> &old_reset_points,
              const size_t bin_size,
              const size_t desert_size,
              const double max_dead_proportion)
{
    vector<SimpleGenomicRegion> bin_boundaries;
    vector<double> read_bins_a;
    vector<double> read_bins_b;
    vector<double> nondead_scales;
    vector<size_t> reset_points;

    reset_points.push_back(0);
    for (size_t i = 0; i < old_reset_points.size() - 1; ++i)
    {
        const size_t start = old_reset_points[i];
        const size_t end = old_reset_points[i + 1];
        size_t n_gaps = 0;
        for (size_t j = start; j < end; ++j)
            if (old_nondead_scales[j] > 1 - max_dead_proportion)
            {
                if (n_gaps * bin_size > desert_size
                    && reset_points.back() < bin_boundaries.size())
                    reset_points.push_back(bin_boundaries.size());

                bin_boundaries.push_back(old_bin_boundaries[j]);
                read_bins_a.push_back(old_read_bins_a[j]);
                read_bins_b.push_back(old_read_bins_b[j]);
                nondead_scales.push_back(old_nondead_scales[j]);

                n_gaps = 0;
            }
            else
            {
                ++n_gaps;
            }
        if (reset_points.back() < bin_boundaries.size())
            reset_points.push_back(bin_boundaries.size());
    }

    assert(bin_boundaries.size() == reset_points.back());
    assert(read_bins_a.size() == reset_points.back());
    assert(read_bins_b.size() == reset_points.back());
    assert(nondead_scales.size() == reset_points.back());

    std::swap(old_bin_boundaries, bin_boundaries);
    std::swap(old_read_bins_a, read_bins_a);
    std::swap(old_read_bins_b, read_bins_b);
    std::swap(old_nondead_scales, nondead_scales);
    std::swap(old_reset_points, reset_points);
}
