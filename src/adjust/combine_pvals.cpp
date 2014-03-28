
#include <vector>

#include <gsl/gsl_rng.h>

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
                     const std::vector<LocusIterator> &neighbors,
                     std::vector< std::vector<double> > &corr_matrix) {
  corr_matrix.clear();
  const size_t num_neighbors = neighbors.size();

  corr_matrix.resize(num_neighbors);

  for (std::vector<std::vector<double> >::iterator row = corr_matrix.begin();
        row != corr_matrix.end(); ++row)
    row->resize(num_neighbors);

  for (size_t row = 0; row < num_neighbors; ++row) {
    corr_matrix[row][row] = 1.0;
    const size_t row_locus = neighbors[row]->end();

    for (size_t col = row + 1; col < num_neighbors; ++col) {
      const size_t col_locus = neighbors[col]->begin();
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
combine_pvals(vector<LocusIterator> &loci_iterators, 
                    const BinForDistance &bin_for_distance) {
                      
  DistanceCorrelation distance_correlation(bin_for_distance);
  vector<double> correlation_for_bin = 
                        distance_correlation.correlation_table(loci_iterators);
  
  //output correlation in each bin
  for (vector<double>::const_iterator it = correlation_for_bin.begin(); 
          it != correlation_for_bin.end(); ++it)
  std::cerr << *it << std::endl;
  
  
  ProximalLoci proximal_loci(loci_iterators, bin_for_distance.max_dist());
  
  vector<double> combined_pvalues;
  
  vector<LocusIterator> neighbors;
  
  size_t i = 0; 

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  while (proximal_loci.get(neighbors)) {
            
   vector< vector<double> > correlation_matrix;
   vector<double> p_vals;
  
   for (vector<LocusIterator>::const_iterator it = neighbors.begin(); 
        it != neighbors.end(); ++it) { 
      double pval = (*it)->score();
      if (pval == 1) {
        //cerr << "replaced " << pval << " with ";
        //pval = 0.999;
        //pval = gsl_rng_uniform(r);
        //cerr << pval << std::endl;
      }
      p_vals.push_back(pval);
    }
    
    distance_corr_matrix(bin_for_distance, correlation_for_bin,
                                neighbors, correlation_matrix);
    
    double combined_pval = stouffer_liptak(p_vals, correlation_matrix);
    Locus &cur_region = *(loci_iterators[i]);
    
    //cerr << cur_region.begin() << " " << cur_region.score() << " " << combined_pval << endl;
    
    if (cur_region.score() > 0.90 && combined_pval < 0.01) {
      cerr << "-----------" << endl;
      cerr << cur_region << " <-- cur_score = " << cur_region.score() << " combined pval = " << combined_pval << endl;
      cerr << "p-vals = ";
      for (vector<double>::const_iterator it = p_vals.begin(); it != p_vals.end(); ++it)
        cerr << *it << " ";
      cerr << endl;
      for (vector<LocusIterator>::const_iterator it = neighbors.begin(); 
           it != neighbors.end(); ++it) 
        cerr << *(*it) << endl;
      cerr << "----correlation matrix----" << endl;
      for (vector<vector<double> >::const_iterator row_it = correlation_matrix.begin(); 
           row_it != correlation_matrix.end(); ++row_it) {
        for (vector<double>::const_iterator it = row_it->begin(); it != row_it->end(); ++it)
          cerr << *it << " ";
        cerr << endl;
      }
      cerr << "-----------" << endl;
    } 
    combined_pvalues.push_back(combined_pval);
    i++;
  }

  gsl_rng_free(r);
  
  assert(loci_iterators.size() == combined_pvalues.size());
  
  for (size_t ind = 0; ind < loci_iterators.size(); ++ind) {
    std::stringstream ss;
    Locus &cur_region = *(loci_iterators[ind]);
    ss << cur_region.name() << ":" << cur_region.score();
    cur_region.set_name(ss.str());
    cur_region.set_score(combined_pvalues[ind]); 
  }
}
