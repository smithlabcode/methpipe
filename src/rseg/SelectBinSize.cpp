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

#include <cmath>
#include <cstdlib>
#include <numeric>

#include "SelectBinSize.hpp"
#include "ReadCounts.hpp"
#include "rseg_utils.hpp"

using std::cerr;
using std::endl;
using std::vector;

static const size_t MIN_BIN_SIZE = 100;
static const size_t MAX_BIN_SIZE = 3000;

static size_t
discretize_bin_size(const double bin_size, const double bin_size_step)
{
    const double n_steps = bin_size / bin_size_step;
    return static_cast<size_t>(
        bin_size_step * (floor(n_steps) + (n_steps - floor(n_steps) > 0.5)));

}

size_t
select_bin_size_waterman(const vector<double> &read_bins,
                         const vector<double> &nondead_scales,
                         const size_t bin_size_step,
                         const bool smooth)
{

  static const double ref_genome_size = 2287968180;
  static const double ref_read_num = 10000000;
  static const double ref_bin_size = smooth ? 600 : 400;
  static const double exponent = smooth ? 1.0/6.0 : 1.0/4.0;

  const double genome_size = bin_size_step
      * std::accumulate(nondead_scales.begin(), nondead_scales.end(), 0.0);

  const double n_reads = std::accumulate(read_bins.begin(), read_bins.end(),
                                         0.0);

  const double b =
      ref_bin_size / pow(n_reads/ref_read_num*ref_genome_size/genome_size,
                         exponent);

  return std::max(MIN_BIN_SIZE,
                  std::min(discretize_bin_size(b, bin_size_step),
                           MAX_BIN_SIZE));
}

size_t
select_bin_size_hideaki(const vector<double> &read_bins,
                        const vector<double> &nondead_scales,
                        const size_t bin_size_step,
                        const bool smooth)
{

  const double ref_genome_size = 2287968180;
  const double ref_read_num = 10000000;
  const double ref_bin_size = smooth ? 600 : 400;
  const double exponent = smooth ? 1.0/3.0 : 1.0/2.0;

  const double genome_size = bin_size_step
      * std::accumulate(nondead_scales.begin(), nondead_scales.end(),
                        0.0);

  const double n_reads = std::accumulate(read_bins.begin(), read_bins.end(),
                                         0.0);

  const double b =
      ref_bin_size / pow(n_reads/ref_read_num*ref_genome_size/genome_size,
                         exponent);
  return std::max(MIN_BIN_SIZE,
                  std::min(discretize_bin_size(b, bin_size_step),
                           MAX_BIN_SIZE));
}

static void
get_mean_var(const vector<double>  &vals, double &mean, double &var)
{
  mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();

  var = 0.0;
  for (size_t i = 0; i < vals.size(); ++i)
    var += pow(vals[i] - mean, 2.0);
  var /= vals.size();
}

size_t
select_bin_size_hideaki_emp(const vector<double> &read_bins,
                            const vector<double> &nondead_scale,
                            const vector<size_t> &reset_points,
                            const size_t bin_size_step,
                            const double max_dead_proportion)
{

  size_t bin_size_low = MIN_BIN_SIZE;
  size_t bin_size_high = MAX_BIN_SIZE;
  size_t range  = bin_size_high > bin_size_low
      ?  bin_size_high - bin_size_low : 0;

  while (range > bin_size_step )
  {
      vector<double> tmp_read_bins;
      vector<double> tmp_nondead_scales;
      vector<size_t> tmp_reset_points;
      vector<double> corrected_read_bins;

      size_t b_one_third(bin_size_low + (bin_size_high - bin_size_low) / 3);
      b_one_third = discretize_bin_size(b_one_third, bin_size_step);
      AdjustBinSize(read_bins, nondead_scale, reset_points, bin_size_step,
                    tmp_read_bins, tmp_nondead_scales, tmp_reset_points,
                    b_one_third);
      GetCorrectedReadCounts(tmp_read_bins, tmp_nondead_scales,
                             corrected_read_bins, max_dead_proportion);

      double mean_one_third, var_one_third, cost_one_third;
      get_mean_var(corrected_read_bins, mean_one_third, var_one_third);
      cost_one_third = (2*mean_one_third - var_one_third)
          / pow(static_cast<double>(b_one_third), 2.0);

      size_t b_two_third(bin_size_high - (bin_size_high - bin_size_low) / 3);
      b_two_third = discretize_bin_size(b_two_third, bin_size_step);
      AdjustBinSize(read_bins, nondead_scale, reset_points, bin_size_step,
                    tmp_read_bins, tmp_nondead_scales, tmp_reset_points,
                    b_two_third);
      GetCorrectedReadCounts(tmp_read_bins, tmp_nondead_scales,
                             corrected_read_bins, max_dead_proportion);

      double mean_two_third, var_two_third, cost_two_third;
      get_mean_var(corrected_read_bins, mean_two_third, var_two_third);
      cost_two_third = (2*mean_two_third - var_two_third)
          / pow(static_cast<double>(b_two_third), 2.0);

      if (cost_one_third > cost_two_third)
          bin_size_low = b_one_third;
      else if (cost_one_third < cost_two_third)
          bin_size_high = b_two_third;
      else if (fabs(cost_one_third - cost_two_third) < 1e-20)
      {
          bin_size_low = b_one_third;
          bin_size_high = b_two_third;
      }
      size_t tmp_range = bin_size_high > bin_size_low
          ?  bin_size_high - bin_size_low : 0;
      if (tmp_range == range)
      {
          bin_size_low += bin_size_step;
          tmp_range = bin_size_high > bin_size_low
          ?  bin_size_high - bin_size_low : 0;
      }
      range = tmp_range;
    }

  return discretize_bin_size((bin_size_low+bin_size_high)/2, bin_size_step);
}













