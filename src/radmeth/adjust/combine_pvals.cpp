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
 
#include <vector>

//#include <gsl/gsl_rng.h>

#include "smithlab_utils.hpp"

#include "combine_pvals.hpp"
#include "distance_correlation.hpp"
#include "proximal_loci.hpp"
#include "stoufferliptak.hpp"

using std::vector; using std::cerr;
using std::endl;


void 
distance_corr_matrix(BinForDistance bin_for_dist,
                     const std::vector<double> &acor_for_bin,
                     const std::vector<PvalLocus> &neighbors,
                     std::vector< std::vector<double> > &corr_matrix) {
  corr_matrix.clear();
  const size_t num_neighbors = neighbors.size();

  corr_matrix.resize(num_neighbors);

  for (std::vector<std::vector<double> >::iterator row = corr_matrix.begin();
        row != corr_matrix.end(); ++row)
    row->resize(num_neighbors);

  for (size_t row = 0; row < num_neighbors; ++row) {
    corr_matrix[row][row] = 1.0;
    const size_t row_locus = neighbors[row].pos + 1;

    for (size_t col = row + 1; col < num_neighbors; ++col) {
      const size_t col_locus = neighbors[col].pos;
      const size_t dist = col_locus - row_locus;
      
      const size_t bin = bin_for_dist.which_bin(dist);

      if (bin == bin_for_dist.invalid_bin()) {
        corr_matrix[row][col] = 0;
      } else
        corr_matrix[row][col] = corr_matrix[col][row] = acor_for_bin[bin];
    }
  }
}

void 
combine_pvals(vector<PvalLocus> &loci, const BinForDistance &bin_for_distance) {

  DistanceCorrelation distance_correlation(bin_for_distance);
  

  vector<double> correlation_for_bin = 
                                  distance_correlation.correlation_table(loci);
  
  //output correlation in each bin
  //for (vector<double>::const_iterator it = correlation_for_bin.begin(); 
  //        it != correlation_for_bin.end(); ++it)
  //std::cerr << *it << std::endl;

  
  ProximalLoci proximal_loci(loci, bin_for_distance.max_dist());
  
  vector<double> combined_pvalues;
  
  vector<PvalLocus> neighbors;
  
  size_t i = 0; 

  
  while (proximal_loci.get(neighbors)) {
   vector< vector<double> > correlation_matrix;
   vector<double> p_vals;
  
   for (vector<PvalLocus>::const_iterator it = neighbors.begin(); 
        it != neighbors.end(); ++it) { 
      double pval = it->raw_pval;
      p_vals.push_back(pval);
    }
    
    distance_corr_matrix(bin_for_distance, correlation_for_bin,
                         neighbors, correlation_matrix);
    
    double combined_pval = stouffer_liptak(p_vals, correlation_matrix);
    
    loci[i].combined_pval = combined_pval;
    //PvalLocus &cur_region = loci[i];

    //combined_pvalues.push_back(combined_pval);
    i++;
  }

  //gsl_rng_free(r);
  
  //assert(loci_iterators.size() == combined_pvalues.size());
  
  /*
  for (size_t ind = 0; ind < loci_iterators.size(); ++ind) {
    std::stringstream ss;
    PvalLocus &cur_region = loci_iterators[ind];
    ss << cur_region.combined_pval = cur_region.score();
    cur_region.set_name(ss.str());
    cur_region.set_score(combined_pvalues[ind]); 
  }*/
}
