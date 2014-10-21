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
};

/*
std::ostream&
operator<<(std::ostream& os, const TableRow& row) {
  os << row.chrom << ":" << row.begin << ":" << row.end;

  for (size_t sample = 0; sample < row.total_counts.size(); ++sample)
    os << " " << row.total_counts[sample]
       << " " << row.meth_counts[sample];

  return os;
}*/


// Populates a TableRow object from a string represetnation.
std::istream&
operator>>(std::istream &table_encoding, TableRow &row) {
  row.chrom.clear();
  row.begin = 0;
  row.end = 0;
  row.meth_counts.clear();
  row.total_counts.clear();

  string row_encoding;
  getline(table_encoding, row_encoding);

  // Skip lines contining only the newline character (e.g. the last line of the
  // proportion table).
  if(row_encoding.empty())
    return table_encoding;

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
  return table_encoding;
}

bool
has_low_coverage(const Design &design, size_t test_factor,
              const TableRow &row) {

  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;

  for (size_t sample = 0; sample < design.sample_names.size(); ++sample) {
    if (design.matrix[sample][test_factor] == 1) {
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

  for (size_t sample = 0; sample < design.sample_names.size(); ++sample) {
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

  Regression full_regression;         // Initialize the full design matrix from
  design_file >> full_regression.design; // file.

  //cerr << full_regression.design << endl;

  // Check that the provided test factor name exists and find it's index. Here
  // we identify with their indexes to simplify naming.
  vector<string>::const_iterator test_factor_it =
    std::find(full_regression.design.factor_names.begin(),
              full_regression.design.factor_names.end(), test_factor_name);

  if (test_factor_it == full_regression.design.factor_names.end())
    throw SMITHLABException("Error: " + test_factor_name +
                            " is not a part of the design specification.");

  size_t test_factor = test_factor_it -
                              full_regression.design.factor_names.begin();

  Regression null_regression;
  null_regression.design = full_regression.design;
  remove_factor(null_regression.design, test_factor);

  //cerr << null_regression.design << endl;

  // Make sure that the first line of the proportion table file contains names
  // of the samples. Throw an exception if the names or their order in
  // the proportion table does not match those in the full design matrix.
  string sample_names_encoding;
  getline(table_file, sample_names_encoding);

  if (full_regression.design.sample_names != split(sample_names_encoding))
    throw SMITHLABException(sample_names_encoding + " does not match factor "
                            "names or their order in the design matrix. "
                            "Please verify that the design matrix and the "
                            "proportion table are correctly formatted.");

  // Performing the log-likelihood ratio test on proportions from each row of
  // the proportion table.
  TableRow row;
  while (table_file >> row) {

    if (full_regression.design.sample_names.size() != row.total_counts.size())
      throw SMITHLABException("There is a row with"
                              "incorrect number of proportions.");

    out << row.chrom << "\t"
        << row.begin << "\t"
        << row.end   << "\t";

    // Do not perform the test if there's no coverage in either all case or all
    // control samples. Also do not test if the site is completely methylated
    // or completely unmethylated across all samples.
    if (has_low_coverage(full_regression.design, test_factor, row)) {
      out << "c:0:0\t" << -1;
    }
    else if (has_extreme_counts(full_regression.design, row)) {
      out << "c:0:0\t" << -1;
    }
    else {
      full_regression.props.total = row.total_counts;
      full_regression.props.meth  = row.meth_counts;
      fit(full_regression);

      null_regression.props.total = row.total_counts;
      null_regression.props.meth = row.meth_counts;
      fit(null_regression);

      const double pval = loglikratio_test(null_regression.max_loglik,
                                     full_regression.max_loglik);

      // If error occured in the fitting algorithm (i.e. p-val is nan or -nan).
      if (pval != pval) {
        out << "c:0:0" << "\t" << "-1";
      }
      else {
        const double
              log_fold_change = full_regression.fitted_parameters[test_factor];
        out << "c:" << log_fold_change
            << ":"  << min_methdiff(full_regression, test_factor)
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
