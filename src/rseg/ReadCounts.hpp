/* Copyright (C) 2011 University of Southern California
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

#ifndef READ_COUNTS_HPP
#define READ_COUNTS_HPP

#include <vector>
#include "GenomicRegion.hpp"


void
AdjustBinSize(const std::vector<double> &old_read_bins,
              const std::vector<double> &old_nondead_scales,
              const std::vector<size_t> &old_reset_points,
              const size_t old_bin_size,
              std::vector<double>  &read_bins,
              std::vector<double>  &nondead_scales,
              std::vector<size_t> &reset_points,
              const size_t bin_size);

void
GetCorrectedReadCounts(const std::vector<double> &read_bins,
                       const std::vector<double> &nondead_scales,
                       std::vector<double> &corrected_read_bins,
                       const double max_dead_proportion);

void
AdjustBinSize(std::vector<SimpleGenomicRegion> &old_bin_boundaries,
              std::vector<double> &old_read_bins,
              std::vector<double> &old_nondead_scales,
              std::vector<size_t> &old_reset_points,
              const size_t old_bin_size,
              const size_t bin_size);

void
RemoveDeserts(std::vector<SimpleGenomicRegion> &old_bin_boundaries,
              std::vector<double> &old_read_bins,
              std::vector<double> &old_nondead_scales,
              std::vector<size_t> &old_reset_points,
              const size_t bin_size,
              const size_t desert_size,
              const double max_dead_proportion);

void
AdjustBinSize(std::vector<SimpleGenomicRegion> &old_bin_boundaries,
              std::vector<double> &old_read_bins_a,
              std::vector<double> &old_read_bins_b,
              std::vector<double> &old_nondead_scales,
              std::vector<size_t> &old_reset_points,
              const size_t old_bin_size,
              const size_t bin_size);

void
RemoveDeserts(std::vector<SimpleGenomicRegion> &old_bin_boundaries,
              std::vector<double> &old_read_bins_a,
              std::vector<double> &old_read_bins_b,
              std::vector<double> &old_nondead_scales,
              std::vector<size_t> &old_reset_points,
              const size_t bin_size,
              const size_t desert_size,
              const double max_dead_proportion);
#endif
