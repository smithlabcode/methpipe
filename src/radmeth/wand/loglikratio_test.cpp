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

// This file contains the implementation of log-likelihood ratio test.

#include <vector>

#include <gsl/gsl_cdf.h>

#include "design.hpp"
#include "regression.hpp"

#include "loglikratio_test.hpp"

using std::vector; using std::size_t;

double loglikratio_test(double null_loglik, double full_loglik) {
  
  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2*(null_loglik - full_loglik);
  
  // It is assumed that null model has one fewer factor than the full model.
  // Hence the number of degrees of freedom is 1.
  const size_t degrees_of_freedom = 1;
  
  // Log-likelihood ratio statistic has a chi-sqare distribution.
  double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double pval = 1.0 - chisq_p;
  
  return pval;
}
