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

// STD headers.
#include <sstream>
#include <vector>

// GSL headers.
#include <gsl/gsl_cdf.h>

// Google Mock headers.
#include "gmock/gmock.h"

// Local headers.
#include "design.hpp"
#include "regression.hpp"
#include "gsl_fitter.hpp"
#include "loglikratio_test.hpp"

using std::istringstream; using std::vector;

using testing::Eq; using testing::DoubleNear;

const double max_abs_error = 0.001;

TEST(a_loglikratio_test, gives_low_pvalue) {

  istringstream design_encoding(
    "f0\tf1\ns0\t1\t0\ns1\t1\t0\ns2\t1\t0\ns3\t1\t0\n"
            "s4\t1\t1\ns5\t1\t1\ns6\t1\t1\ns7\t1\t1");

  vector<size_t> response_total(8, 15);
  vector<size_t> response_meth = { 8, 5, 2, 2, 13, 14, 8, 5};
  
  const size_t test_factor = 1;

  Design design(design_encoding);
  Regression regression(design);

  Design null_design = design;
  null_design.remove_factor(test_factor);
  Regression null_regression(null_design);
    
  regression.set_response(response_total, response_meth);
  gsl_fitter(regression);

  null_regression.set_response(response_total, response_meth);
  gsl_fitter(null_regression);

  double pval = loglikratio_test(null_regression.maximum_likelihood(), 
                                regression.maximum_likelihood());
  
  ASSERT_THAT(pval, DoubleNear(0.0302161 , max_abs_error));
}

TEST(a_loglik_ratio_test, gives_high_pvalue) {
  istringstream design_encoding(
    "f0\tf1\ns0\t1\t0\ns1\t1\t0\ns2\t1\t0\ns3\t1\t0\n"
            "s4\t1\t1\ns5\t1\t1\ns6\t1\t1\ns7\t1\t1");
  
  vector<size_t> response_total(8, 15);
  vector<size_t> response_meth = { 7, 11,  1,  6,  5,  4,  2,  4};
  
  const size_t test_factor = 1;
  
  Design design(design_encoding);
  Regression regression(design);

  Design null_design = design;
  null_design.remove_factor(test_factor);
  Regression null_regression(null_design);
    
  regression.set_response(response_total, response_meth);
  gsl_fitter(regression);

  null_regression.set_response(response_total, response_meth);
  gsl_fitter(null_regression);

  double pval = loglikratio_test(null_regression.maximum_likelihood(), 
                                regression.maximum_likelihood());
                                  
  ASSERT_THAT(pval, DoubleNear(0.276009 , max_abs_error));
}
