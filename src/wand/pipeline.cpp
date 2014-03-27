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
 
// Contains the implementation of the wand function.

// STD headers. 
#include <iostream>
#include <algorithm>

// Smithlab headers.
#include "smithlab_utils.hpp"

// Local headers.
#include "design.hpp"
#include "pipeline.hpp"
#include "table_row.hpp"
#include "loglikratio_test.hpp"
#include "gsl_fitter.hpp"
#include "regression.hpp"

using std::istream;
using std::istringstream;
using std::ostream;
using std::string;
using std::vector;

vector<string> 
split(string input) {
  istringstream iss(input);
  string token;
  vector<string> tokens;
  
  while (iss >> token)
    tokens.push_back(token);
  
  return tokens;
}

void
wand(istream &design_encoding, istream &table_encoding, string test_factor_name, 
      ostream &out) {
  
  Design full_design(design_encoding);
  
  vector<string> factor_names = full_design.factor_names();
  
  vector<string>::const_iterator test_factor_it = 
      std::find(factor_names.begin(), factor_names.end(), test_factor_name);
  
  // Checking that the provided test factor names exitst.    
  if (test_factor_it == factor_names.end())
    throw SMITHLABException(test_factor_name + " is not a part of the design" 
                            " specification.");
  
  // Factors are identified with their indexes to simplify naming.
  size_t test_factor = test_factor_it - factor_names.begin();
  
  Regression full_regression(full_design);
  Design null_design = full_design;
  null_design.remove_factor(test_factor);
  Regression null_regression(null_design);
  
  // Read the first line of the count table which must contain names of the
  // samples.
  string sample_names_encoding;
  getline(table_encoding, sample_names_encoding);
  
  if (  full_design.sample_names() != split(sample_names_encoding) )
    throw SMITHLABException(sample_names_encoding + " does not match factor "
                            "names (or their order) in the design matrix.");
  
  string row_encoding;
  
  // Perform the log-likelihood ratio on proportions from every row of the 
  // proportion table.
  while( getline(table_encoding, row_encoding) ) {
    
    TableRow row;
    read_row(row_encoding, row);
    
    assert_compatibility(full_design, row);
                                    
    out << row.chrom << "\t" 
        << row.begin << "\t"
        << row.end   << "\t";
    
    if (has_low_coverage(full_design, test_factor, row)) {
      // Don't test if there's no coverage in either case or control samples. 
      out << "c:0:0\t" << -1;
    } else if (has_extreme_counts(full_design, row)) {
      // Don't test if the site is completely methylated or completely 
      // unmethylated across all samples.
      out << "c:0:0\t" << -1;
    } else {
      full_regression.set_response(row.total_counts, row.meth_counts);
      gsl_fitter(full_regression);

      null_regression.set_response(row.total_counts, row.meth_counts);
      gsl_fitter(null_regression);

      double pval = loglikratio_test(null_regression.maximum_likelihood(), 
                                     full_regression.maximum_likelihood());
      
      // If error occured in the fitting algorithm (i.e. p-val is nan or -nan).
      if (pval != pval) {
        out << "c:0:0" << "\t" << "-1";
      } else {
        out << "c:" << full_regression.log_fold_change(test_factor)
            << ":"  << full_regression.min_methdiff(test_factor)
            << "\t" << pval;
      }
    }
    out << std::endl;
  }
}
