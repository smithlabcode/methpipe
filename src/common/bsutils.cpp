/*
 *    Copyright (C) 2018 Andrew D. Smith
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

#include "bsutils.hpp"

#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

#include <GenomicRegion.hpp>
#include <smithlab_utils.hpp>

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
      throw std::runtime_error("read sequence not found for: " + names[i]);
    lookup.push_back(j->second);
  }
}
