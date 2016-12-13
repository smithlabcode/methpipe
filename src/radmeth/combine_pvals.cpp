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
#include <iterator>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <numeric>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

#include "smithlab_utils.hpp"
#include "combine_pvals.hpp"

using std::vector;
using std::cerr;
using std::endl;
using std::copy;
using std::vector;
using std::istream;
using std::string;

static double
to_zscore(double pval) {
  static const double local_epsilon = 1e-6;

  if (pval > 1.0 - local_epsilon)
    pval = 1.0 - local_epsilon;
  else if (pval < local_epsilon)
    pval = local_epsilon;

  return gsl_cdf_ugaussian_Pinv(1.0 - pval);
}

static double
stouffer_liptak(const vector<vector<double> > &corr_mat, vector<double> &pvals) {

  double correction = 0;
  for (size_t row_ind = 0; row_ind < corr_mat.size(); ++row_ind)
    for (size_t col_ind = row_ind + 1; col_ind < corr_mat.size(); ++col_ind)
      correction += corr_mat[row_ind][col_ind];

  vector<double> zscores;
  transform(pvals.begin(), pvals.end(), std::back_inserter(zscores), to_zscore);

  const double sum = std::accumulate(zscores.begin(), zscores.end(), 0.0);
  const double test_stat =
    sum/sqrt(static_cast<double>(pvals.size()) + 2.0*correction);

  return 1.0 - gsl_cdf_gaussian_P(test_stat, 1.0);
}

void
update_pval_loci(std::istream &input_encoding,
                 const vector<PvalLocus> &pval_loci,
                 std::ostream &output_encoding) {

  string record, chrom, name, sign;
  size_t position, coverage_factor, meth_factor, coverage_rest, meth_rest;
  double pval;

  vector<PvalLocus>::const_iterator cur_locus_iter = pval_loci.begin();

  while (getline(input_encoding, record)) {
    // ADS: this seems not to be done well; the code should exit in a
    // "normal" state if bad parse
    try {
      std::istringstream iss(record);
      iss.exceptions(std::ios::failbit);
      iss >> chrom >> position >> sign >> name >> pval
          >> coverage_factor >> meth_factor >> coverage_rest >> meth_rest;
    }
    catch (std::exception const & err) {
      cerr << err.what() << endl << "could not parse line:\n"
           << record << endl;
      std::terminate();
    }

    output_encoding << chrom << "\t" << position << "\t" << sign << "\t"
                    << name << "\t" << pval << "\t";

    if (0.0 <= pval && pval <= 1.0) {
      output_encoding << cur_locus_iter->combined_pval << "\t"
                      << cur_locus_iter->corrected_pval << "\t";
      cur_locus_iter++;
    }
    else output_encoding << -1 << "\t" << -1 << pval << "\t"; // MAGIC??

    output_encoding << coverage_factor << "\t" << meth_factor << "\t"
                    << coverage_rest << "\t" << meth_rest << endl;
  }
}

BinForDistance::BinForDistance(std::string spec_string) {
  std::replace(spec_string.begin(), spec_string.end(), ':', ' ');

  std::istringstream iss(spec_string);
  iss >> min_dist_ >> max_dist_ >> bin_size_;

  num_bins_ = (max_dist_ - min_dist_) / bin_size_;
  invalid_bin_ = num_bins_ + 1;
}

size_t
BinForDistance::which_bin(size_t value) const {
  if (value < min_dist_)
    return invalid_bin_;

  const size_t bin = (value - min_dist_)/bin_size_;

  //Bin numbering is 0 based.
  if (bin >= num_bins_)
    return invalid_bin_;

  return bin;
}

bool
ProximalLoci::get(vector<PvalLocus> &neighbors) {

  if (next_pos_ == loci_.end())
    return false;

  vector<PvalLocus>::const_iterator cur_pos = next_pos_;
  neighbors.clear();
  neighbors.push_back(*cur_pos);

  if ( cur_pos != loci_.begin() ) {
    vector<PvalLocus>::const_iterator up_pos = cur_pos;
    bool too_far = false;

    do {
      --up_pos;
      size_t up_dist = cur_pos->pos - (up_pos->pos + 1);

      if(up_dist <= max_distance_) {
          neighbors.push_back(*up_pos);
      } else
        too_far = true;

    } while (!too_far && up_pos != loci_.begin());
  }

  std::reverse(neighbors.begin(), neighbors.end());

  if (cur_pos != loci_.end() - 1) {
    bool too_far = false;
    vector<PvalLocus>::const_iterator down_pos = cur_pos;

    do {
      ++down_pos;
      size_t down_dist = down_pos->pos - (cur_pos->pos + 1);

      if (down_dist <= max_distance_) {
          neighbors.push_back(*down_pos);
      }
      else too_far = true;

    } while (!too_far && down_pos != loci_.end() - 1);
  }

  ++next_pos_;
  return true;
}

void
DistanceCorrelation::bin(const vector<PvalLocus> &loci) {
  x_pvals_for_bin_.clear();
  y_pvals_for_bin_.clear();
  vector<PvalLocus>::const_iterator it = loci.begin();

  while (it != loci.end()) {
    vector<PvalLocus>::const_iterator forward_it = it + 1;
    bool too_far = false;

    while (forward_it != loci.end() && !too_far) {
      const size_t dist = forward_it->pos - (it->pos + 1);
      const size_t bin = bin_for_dist_.which_bin(dist);

      //check if the appropriate bin exists
      if (bin != bin_for_dist_.invalid_bin()) {
        x_pvals_for_bin_[bin].push_back(to_zscore(it->raw_pval));
        y_pvals_for_bin_[bin].push_back(to_zscore(forward_it->raw_pval));
      }

      if (dist > bin_for_dist_.max_dist())
        too_far = true;

      ++forward_it;
    }

    ++it;
  }
}

double
DistanceCorrelation::correlation(const vector<double> &x,
                                    const vector<double> &y) {
  //Correlation is 0 when all bins are empty.
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
DistanceCorrelation::correlation_table(const vector<PvalLocus> &loci) {
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
    double combined_pval = stouffer_liptak(correlation_matrix, p_vals);
    loci[i].combined_pval = combined_pval;

    i++;
  }
}
