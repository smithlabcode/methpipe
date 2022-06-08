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

#ifndef RSEG_UTILS_HPP
#define RSEG_UTILS_HPP

#include <vector>
#include <string>

#include "GenomicRegion.hpp"
#include "Distro.hpp"
#include "SplitDistro.hpp"

void
pick_training_sample(const std::vector<double> &read_bins,
                     const std::vector<double> &read_bins_a,
                     const std::vector<double> &read_bins_b,
                     const std::vector<double> &scales,
                     const std::vector<size_t> &reset_points,
                     const size_t training_sample_size,
                     std::vector<double> &read_bins_sample,
                     std::vector<double> &read_bins_a_sample,
                     std::vector<double> &read_bins_b_sample,
                     std::vector<double> &scales_sample,
                     std::vector<size_t> &reset_points_sample);
void
clear_training_sample(std::vector<double> &read_bins_sample,
		      std::vector<double> &read_bins_a_sample,
		      std::vector<double> &read_bins_b_sample,
		      std::vector<double> &scales_sample,
              std::vector<size_t> &reset_points_sample);

void
set_transitions(const size_t bin_size, const double fg_size,
		const std::vector<double> &mixing,
		const bool VERBOSE,
		std::vector<double> &start_trans,
		std::vector<std::vector<double> > &trans,
		std::vector<double> &end_trans);
void
set_transitions(const size_t bin_size, const double fg_size,
		const double mixing, const bool VERBOSE,
		std::vector<double> &start_trans,
		std::vector<std::vector<double> > &trans,
		std::vector<double> &end_trans);

void
report_final_values(const std::vector<Distro> &distros,
		    const std::vector<double> &start_trans,
		    const std::vector<std::vector<double> > &trans,
		    const std::vector<double> &end_trans);
void
report_final_values(const std::vector<SplitDistro> &distros,
		    const std::vector<double> &start_trans,
		    const std::vector<std::vector<double> > &trans,
		    const std::vector<double> &end_trans);

void
chk_and_mk_dirs(const std::string & path);

void
write_read_counts_by_bin(const std::vector< std::vector<SimpleGenomicRegion> > &bin_boundaries,
                         const std::vector<double> &read_bins,
                         const std::vector<double> &scales,
                         const std::vector<bool> &classes,
                         const std::string &file_name,
                         const bool VERBOSE = false);
void
write_read_counts_by_bin(const std::vector< std::vector<SimpleGenomicRegion> > &bin_boundaries,
                         const std::vector<double> &read_bins,
                         const std::vector<double> &read_bins_a,
                         const std::vector<double> &read_bins_b,
                         const std::vector<bool> &classes,
                         const std::string &file_name,
                         const bool VERBOSE = false);
void
write_read_counts_by_bin(const std::vector< std::vector<SimpleGenomicRegion> > &bin_boundaries,
                         const std::vector<double> &read_bins,
                         const std::vector<double> &read_bins_a,
                         const std::vector<double> &read_bins_b,
                         const std::vector<size_t> &classes,
                         const std::string &file_name,
                         const bool VERBOSE = false);

std::string
strip_path_and_bed_suffix(const std::string &full_path);

void
write_wigfile(const std::vector<std::vector<double> > &scores,
	      const std::vector<std::vector<SimpleGenomicRegion> > &bin_bounds,
	      const std::string &wigfile_name);
void
write_wigfile(const std::vector<double > &fg_scores,
              const std::vector<double > &bg_scores,
              const std::vector<std::vector<SimpleGenomicRegion> > &bin_bounds,
              const std::string &wigfile_name);

void
write_bed_file(const std::vector<std::vector<GenomicRegion> > &regions,
	       const std::string &bed_file);

template <class T> void
expand_bins(const std::vector<T> &tmp_bins, const std::vector<size_t> &resets,
	    std::vector<std::vector<T> > &bins)
{
  bins.clear();
  size_t j = 0;
  for (size_t i = 0; i < tmp_bins.size(); ++i)
    {
      if (i == resets[j])
        {
	  bins.push_back(std::vector<T>());
	  ++j;
        }
      bins.back().push_back(tmp_bins[i]);
    }
}

void
build_domains(const std::vector<std::vector<SimpleGenomicRegion> > &bins,
              const std::vector<std::vector<bool> > &classes,
              const std::vector<std::vector<double> > &scores, // posterior score of classes[i]
              const double score_cutoff,
              std::vector<std::vector<GenomicRegion> > &domains,
              const size_t undef_domain_cutoff);

void
build_domains(const std::vector<std::vector<SimpleGenomicRegion> > &bins,
              const std::vector<std::vector<size_t> > &classes,
              const std::vector<std::vector<double> > &scores, // posterior score of classes[i]
              const double score_cutoff,
              std::vector<std::vector<GenomicRegion> > &domains,
              const size_t undef_domain_cutoff);

void
pick_domains(const std::vector<std::vector<SimpleGenomicRegion> > &bins,
             const std::vector<std::vector<double> > &read_counts,
             const std::vector<std::vector<double> > &scales,
             const std::vector<Distro> &distros,
             std::vector<std::vector<GenomicRegion> > &domains,
             const double p_value);

void
pick_domains(const std::vector<std::vector<SimpleGenomicRegion> > &bins,
             const std::vector<std::vector<double> > &read_counts,
             const std::vector<std::vector<double> > &scales,
             const std::vector<SplitDistro> &distros,
             std::vector<std::vector<GenomicRegion> > &domains,
             const double cdf_cutoff);

void
pick_domains_3s(const std::vector<std::vector<SimpleGenomicRegion> > &bins,
                const std::vector<std::vector<double> > &read_counts,
                const std::vector<std::vector<double> > &scales,
                const std::vector<SplitDistro> &distros,
                std::vector<std::vector<GenomicRegion> > &domains,
                const double cdf_cutoff);
#endif
