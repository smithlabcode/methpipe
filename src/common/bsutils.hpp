/*
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BSUTILS_HPP
#define BSUTILS_HPP

#include <string>
#include <vector>
#include <cmath>

#include <GenomicRegion.hpp>
#include <smithlab_utils.hpp>

inline bool
is_cytosine(char c) {return (c == 'c' || c == 'C');}

inline bool
is_guanine(char c)  {return (c == 'g' || c == 'G');}

inline bool
is_thymine(char c)  {return (c == 't' || c == 'T');}

inline bool
is_adenine(char c)  {return (c == 'a' || c == 'A');}


//// CONFIDENCE INTERVALS //**************////////////////////////
#include <gsl/gsl_cdf.h>
void
wilson_ci_for_binomial(const double alpha, const double n,
                       const double p_hat, double &lower, double &upper) {
  const double z = gsl_cdf_ugaussian_Pinv(1 - alpha/2);
  const double denom = 1 + z*z/n;
  const double first_term = p_hat + z*z/(2*n);
  const double discriminant = p_hat*(1 - p_hat)/n + z*z/(4*n*n);
  lower = std::max(0.0, (first_term - z*std::sqrt(discriminant))/denom);
  upper = std::min(1.0, (first_term + z*std::sqrt(discriminant))/denom);
}
//////////////////////////**************////////////////////////


bool
is_cpg(const std::string &s, size_t i) {
  return (i < (s.length() - 1)) &&
    is_cytosine(s[i]) && is_guanine(s[i + 1]);
}


void
adjust_region_ends(const std::vector<std::vector<GenomicRegion> > &clusters,
                   std::vector<GenomicRegion> &regions) {
  assert(clusters.size() == regions.size());
  for (size_t i = 0; i < regions.size(); ++i) {
    size_t max_pos = regions[i].get_end();
    size_t min_pos = regions[i].get_start();
    for (size_t j = 0; j < clusters[i].size(); ++j) {
      max_pos = std::max(clusters[i][j].get_end(), max_pos);
      min_pos = std::min(clusters[i][j].get_start(), min_pos);
    }
    regions[i].set_end(max_pos);
    regions[i].set_start(min_pos);
  }
}


void
relative_sort(const std::vector<GenomicRegion> &mapped_locations,
              const std::vector<std::string> &names,
              std::vector<size_t> &lookup) {

  std::unordered_map<std::string, size_t> names_map;
  for (size_t i = 0; i < names.size(); ++i)
    names_map[names[i]] = i;

  for (size_t i = 0; i < mapped_locations.size(); ++i) {
    const std::unordered_map<std::string, size_t>::const_iterator
      j(names_map.find(mapped_locations[i].get_name()));
    if (j == names_map.end())
      throw SMITHLABException("read sequence not found for: " + names[i]);
    lookup.push_back(j->second);
  }
}


template <class T, class U, class V> static void
separate_regions(const std::vector<T> &big_regions,
                 const std::vector<U> &regions,
                 const std::vector<V> &seqs,
                 std::vector<std::vector<U> > &sep_regions,
                 std::vector<std::vector<V> > &sep_seqs) {
  size_t rr_id = 0;
  const size_t n_regions = regions.size();
  assert(n_regions <= seqs.size());

  const size_t n_big_regions = big_regions.size();
  sep_regions.resize(n_big_regions);
  sep_seqs.resize(n_big_regions);
  for (size_t i = 0; i < n_big_regions; ++i) {
    const std::string current_chrom(big_regions[i].get_chrom());
    const size_t current_start = big_regions[i].get_start();
    const size_t current_end = big_regions[i].get_end();
    while (rr_id < n_regions &&
           (regions[rr_id].get_chrom() < current_chrom ||
            (regions[rr_id].get_chrom() == current_chrom &&
             regions[rr_id].get_end() <= current_start)))
      ++rr_id;
    while (rr_id < n_regions &&
           (regions[rr_id].get_chrom() == current_chrom &&
            regions[rr_id].get_start() < current_end)) {
      sep_regions[i].push_back(regions[rr_id]);
      sep_seqs[i].push_back(seqs[rr_id]);
      ++rr_id;
    }
  }
}


#endif
