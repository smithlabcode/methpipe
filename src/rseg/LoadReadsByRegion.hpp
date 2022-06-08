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

#ifndef LOAD_READS_BY_REGION_HPP
#define LOAD_READS_BY_REGION_HPP

#include <string>
#include <vector>

#include "GenomicRegion.hpp"

void
LoadReadsByRegion(const bool VERBOSE,
                  const std::string &regions_file,
                  const std::string &reads_file,
                  const std::string &deads_file,
                  const size_t bin_size,
                  std::vector<SimpleGenomicRegion> &boundaries,
                  std::vector<double> &read_bins,
                  std::vector<double> &nondead_scales,
                  std::vector<size_t> &reset_points,
                  const size_t FRAGMENT_LEN = 0,
                  const bool BAM_FORMAT = false,
                  const bool REMOVE_JACKPOT = true);

void
LoadReadsByRegion(const bool VERBOSE,
                  const std::string &regions_file,
                  const std::string &reads_file_a,
                  const std::string &reads_file_b,
                  const std::string &deads_file,
                  const size_t bin_size,
                  std::vector<SimpleGenomicRegion> &boundaries,
                  std::vector<double> &read_bins_a,
                  std::vector<double> &read_bins_b,
                  std::vector<double> &nondead_scales,
                  std::vector<size_t> &reset_points,
                  const size_t FRAGMENT_LEN = 0,
                  const bool BAM_FORMAT = false,
                  const bool REMOVE_JACKPOT = true);
#endif
