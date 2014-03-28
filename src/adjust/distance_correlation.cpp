/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
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
 */
 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

#include "distance_correlation.hpp"

using std::cout;  using std::endl;
using std::cerr;  using std::vector;
using std::string;  using std::pair;
using std::vector;

static bool 
valid_score(const vector<LocusIterator>::const_iterator &it) {
  return 0 <= (*it)->score() && (*it)->score() <= 1;
}

static double 
to_zscore(double pval) {
  if (pval == 1)
    pval = 0.9999;
    
  if (pval == 0)
    pval = 1 - 0.9999;
    
  return gsl_cdf_ugaussian_Pinv(1 - pval);
}

void
DistanceCorrelation::bin(const vector<LocusIterator> &loci) {
  
  x_pvals_for_bin_.clear();
  y_pvals_for_bin_.clear();
  
  vector<LocusIterator>::const_iterator it = loci.begin();
  
  while (it != loci.end()) {
    //Skip invalid scores
    if (valid_score(it)) {
      
      vector<LocusIterator>::const_iterator forward_it = it + 1;
      
      bool too_far = false;
      
      while (forward_it != loci.end() && !too_far) {
        const size_t dist = (*forward_it)->begin() - (*it)->end();
        const size_t bin = bin_for_dist_.which_bin(dist);
        
        //check if the appropriate bin exists
        if (bin != bin_for_dist_.invalid_bin()) {
          //Skip invalid scores
          if (valid_score(forward_it)) {
            x_pvals_for_bin_[bin].push_back(to_zscore((*it)->score()));
            y_pvals_for_bin_[bin].push_back(to_zscore((*forward_it)->score()));
          }
        }
        
        if (dist > bin_for_dist_.max_dist())
          too_far = true;

        ++forward_it;
      }
      
    }
    ++it;
  }
}

double
DistanceCorrelation::correlation(const vector<double> &x,
                                    const vector<double> &y) {

  //correlation is 0 for all empty bins
  if (x.size() <= 1)
    return 0;
    
  gsl_vector_const_view gsl_x = gsl_vector_const_view_array(&x[0], x.size());
  gsl_vector_const_view gsl_y = gsl_vector_const_view_array(&y[0], y.size());

  const size_t stride = 1;
  double corr = gsl_stats_correlation( gsl_x.vector.data, stride,
                                         gsl_y.vector.data, stride,
                                         x.size());
  return corr;
}

vector<double>
DistanceCorrelation::correlation_table(const vector<LocusIterator> &loci) {

    const size_t num_bins = bin_for_dist_.num_bins();
    x_pvals_for_bin_.resize(num_bins);
    y_pvals_for_bin_.resize(num_bins);

    bin(loci);

    vector<double> correlation_table;

    for (size_t bin = 0; bin < num_bins; ++bin) {

      const double corr = correlation(x_pvals_for_bin_[bin],
                                            y_pvals_for_bin_[bin]);

      correlation_table.push_back(corr);
    }

    return correlation_table;
}

