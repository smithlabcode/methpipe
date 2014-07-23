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
 
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <gsl/gsl_cdf.h>

#include "stoufferliptak.hpp"

using std::for_each;  using std::vector;
using std::cout;  using std::endl;
using std::back_inserter;

static double 
pval_to_zscores(double x) {
  if (x == 1)
    x = 0.9999;
  if (x == 0)
    x = 1 - 0.9999;
  return gsl_cdf_ugaussian_Pinv(1 - x);
}

double 
stouffer_liptak(std::vector<double> &pvals, 
                const std::vector< std::vector<double> > &cor_matrix) {
  double correction = 0;
  size_t num_pvals = pvals.size();
  for (size_t row_ind = 0; row_ind < cor_matrix.size(); ++row_ind)
    for (size_t col_ind = row_ind + 1; col_ind < cor_matrix.size(); ++col_ind)
      correction += cor_matrix[row_ind][col_ind];
  
  vector<double> zscores;
  
  transform(pvals.begin(), pvals.end(), back_inserter(zscores), 
            pval_to_zscores);
  
  double sum = 0;
  for (size_t ind = 0; ind < num_pvals; ++ind)
    sum += zscores[ind];
  
  double test_statistic = sum/sqrt(double(num_pvals) + 2*correction);
  
  return 1 - gsl_cdf_gaussian_P(test_statistic, 1);
}