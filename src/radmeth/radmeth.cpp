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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// GSL headers
#include <gsl/gsl_cdf.h>

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

// Local headers.
#include "gsl_fitter.hpp"
#include "regression.hpp"

using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;

// Splits a string using white-space characters as delimeters.
static vector<string>
split(string input) {
  istringstream iss(input);
  string token;
  vector<string> tokens;

  while (iss >> token)
    tokens.push_back(token);

  return tokens;
}

// Parses a natural number from its string representation. Throws exception if
// the string does not encode one.
static size_t
parse_natural_number(string encoding) {
  istringstream iss(encoding);
  size_t number;
  iss >> number;
  if (!iss)
    throw SMITHLABException("The token \"" +encoding + "\" "
                            "does not encode a natural number");
  return number;
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
double
loglikratio_test(double null_loglik, double full_loglik) {

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

// Stores a row of a proportion table.
struct TableRow {
  std::string chrom;
  size_t begin;
  size_t end;
  std::vector<size_t> meth_counts;
  std::vector<size_t> total_counts;

  friend std::ostream& operator<<(std::ostream& os, const TableRow& row) {
    os << row.chrom << ":" << row.begin << ":" << row.end;

    for (size_t sample = 0; sample < row.total_counts.size(); ++sample)
      os << " " << row.total_counts[sample]
         << " " << row.meth_counts[sample];

    return os;
  }
};

// Populates a TableRow object from a string represetnation.
static void
parse_row(const string &row_encoding, TableRow &row) {
  // Get the row name (which must be specified like this: "chr:start:end") and
  // parse it.
  istringstream row_stream(row_encoding);
  string row_name_encoding;
  row_stream >> row_name_encoding;

  // Every row must start an identifier consisiting of genomic loci of the
  // corresponding site. Here we check this identifier has the correct number
  // of colons.
  const size_t num_colon =
            std::count(row_name_encoding.begin(), row_name_encoding.end(), ':');

  if (num_colon != 2)
    throw SMITHLABException("Each row in the count table must start with "
                            "a line chromosome:start:end. Got \"" +
                            row_name_encoding + "\" instead." );

  // First parse the row identifier.
  istringstream name_stream(row_name_encoding);
  getline(name_stream, row.chrom, ':');

  if (row.chrom.empty())
    throw SMITHLABException("Error parsing " + row_name_encoding +
                            ": chromosome name is missing.");

  string coordinate_encoding;

  getline(name_stream, coordinate_encoding, ':');
  row.begin = parse_natural_number(coordinate_encoding);

  name_stream >> coordinate_encoding;
  row.end = parse_natural_number(coordinate_encoding);

  // After parsing the row identifier, parse count proportions.
  size_t total_count, meth_count;

  while (row_stream >> total_count >> meth_count) {
    row.total_counts.push_back(total_count);
    row.meth_counts.push_back(meth_count);
  }

  if (!row_stream.eof())
    throw SMITHLABException("Some row entries are not natural numbers: " +
                            row_stream.str());

  if (row.total_counts.size() != row.meth_counts.size())
    throw SMITHLABException("This row does not encode proportions"
                            "correctly:\n" + row_encoding);
}

bool
has_low_coverage(const Design &design, size_t test_factor,
              const TableRow &row) {

  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;

  for (size_t sample = 0; sample < design.num_samples(); ++sample) {
    if (design(sample, test_factor) == 1) {
      if (row.total_counts[sample] != 0)
        is_covered_in_test_factor_samples = true;
    } else {
      if (row.total_counts[sample] != 0)
        is_covered_in_other_samples = true;
    }
  }

  return !is_covered_in_test_factor_samples || !is_covered_in_other_samples;
}

bool
has_extreme_counts(const Design &design, const TableRow &row) {

  bool is_maximally_methylated = true;
  bool is_unmethylated = true;

  for (size_t sample = 0; sample < design.num_samples(); ++sample) {
    if (row.total_counts[sample] != row.meth_counts[sample])
      is_maximally_methylated = false;

    if (row.meth_counts[sample] != 0)
      is_unmethylated = false;
  }

  return is_maximally_methylated || is_unmethylated;
}

int
main(int argc, const char **argv) {

  try {
    string outfile;
    string test_factor_name;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Produces multi-factor "
          "differential methylation scores", "<design-matrix> <data-matrix>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);

    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    opt_parse.add_opt("factor", 'f', "a factor to test",
                      true, test_factor_name);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string design_filename(leftover_args.front());
    const string table_filename(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream design_file(design_filename.c_str());
    if (!design_file)
      throw SMITHLABException("could not open file: " + design_filename);

    std::ifstream table_file(table_filename.c_str());
    if (!table_file)
      throw SMITHLABException("could not open file: " + table_filename);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  Design full_design(design_file); // Initialize the full design matrix from
                                   // file.

  vector<string> factor_names = full_design.factor_names();

  // Check that the provided test factor name exists and find it's index. Here
  // we identify with their indexes to simplify naming.
  vector<string>::const_iterator test_factor_it =
    std::find(factor_names.begin(), factor_names.end(), test_factor_name);

  if (test_factor_it == factor_names.end())
    throw SMITHLABException(test_factor_name + " is not a part of the design"
                            " specification.");

  size_t test_factor = test_factor_it - factor_names.begin();

  Regression full_regression(full_design);
  Design null_design = full_design;
  null_design.remove_factor(test_factor);
  Regression null_regression(null_design);

  // Make sure that the first line of the proportion table file contains names
  // of the samples. Throw an exception if the names or their order in
  // the proportion table does not match those in the full design matrix.
  string sample_names_encoding;
  getline(table_file, sample_names_encoding);

  if (full_design.sample_names() != split(sample_names_encoding))
    throw SMITHLABException(sample_names_encoding + " does not match factor "
                            "names or their order in the design matrix. "
                            "Please verify that the design matrix and the "
                            "proportion table are correctly formatted.");

  string row_encoding;

  // Performing the log-likelihood ratio test on proportions from each row of
  // the proportion table.
  while (getline(table_file, row_encoding)) {

    TableRow row;
    parse_row(row_encoding, row);

    if (full_design.num_samples() != row.total_counts.size())
      throw SMITHLABException("This row has incorrect number of proportions:\n"
                              + row_encoding);

    out << row.chrom << "\t"
        << row.begin << "\t"
        << row.end   << "\t";

    // Do not perform the test if there's no coverage in either all case or all
    // control samples. Also do not test if the site is completely methylated
    // or completely unmethylated across all samples.
    if (has_low_coverage(full_design, test_factor, row)) {
      out << "c:0:0\t" << -1;
    }
    else if (has_extreme_counts(full_design, row)) {
      out << "c:0:0\t" << -1;
    }
    else {
      full_regression.set_response(row.total_counts, row.meth_counts);
      gsl_fitter(full_regression);

      null_regression.set_response(row.total_counts, row.meth_counts);
      gsl_fitter(null_regression);

      const double pval = loglikratio_test(null_regression.maximum_likelihood(),
                                     full_regression.maximum_likelihood());

      // If error occured in the fitting algorithm (i.e. p-val is nan or -nan).
      if (pval != pval) {
        out << "c:0:0" << "\t" << "-1";
      }
      else {
        out << "c:" << full_regression.log_fold_change(test_factor)
            << ":"  << full_regression.min_methdiff(test_factor)
            << "\t" << pval;
      }
    }
    out << endl;
  }

  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
