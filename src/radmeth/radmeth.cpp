/* Copyright (C) 2013 University of Southern California and
 *                    Egor Dolzhenko
 *                    Andrew D Smith
 *
 * Authors: Andrew D. Smith and Egor Dolzhenko
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// GSL headers
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::cout;
using std::endl;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::runtime_error;
using std::sort;
using std::min;
using std::begin;
using std::end;

/***************** REGRESSION *****************/
struct Design {
  vector<string> factor_names;
  vector<string> sample_names;
  vector<vector<double> > matrix;
};

struct SiteProportions {
  string chrom;
  size_t position;
  string strand;
  string context;
  vector<size_t> total;
  vector<size_t> meth;
};

struct Regression {
  Design design;
  SiteProportions props;
  double max_loglik;
};

istream&
operator>>(istream &is, Design &design) {
  string header_encoding;
  getline(is, header_encoding);

  istringstream header_is(header_encoding);
  string header_name;
  while (header_is >> header_name)
    design.factor_names.push_back(header_name);

  string row;
  while (getline(is, row)) {

    if (row.empty())
      continue;

    istringstream row_is(row);
    string token;
    row_is >> token;
    design.sample_names.push_back(token);

    vector<double> matrix_row;
    while (row_is >> token) {
      if (token.length() == 1 && (token == "0" || token == "1"))
        matrix_row.push_back(token == "1");
      else
        throw runtime_error("only binary factor levels are allowed:\n"
                            + row);
    }

    if (matrix_row.size() != design.factor_names.size())
      throw runtime_error("each row must have as many columns as "
                          "factors:\n" + row);

    design.matrix.push_back(vector<double>());
    swap(design.matrix.back(), matrix_row);
  }
  return is;
}

ostream&
operator<<(ostream &os, const Design &design) {
  for(size_t factor = 0; factor < design.factor_names.size(); ++factor) {
    os << design.factor_names[factor];
    if (factor + 1 != design.factor_names.size())
      os << "\t";
  }
  os << endl;

  for(size_t sample = 0; sample < design.sample_names.size(); ++sample) {
    os << design.sample_names[sample] << "\t";
    for(size_t factor = 0; factor < design.factor_names.size(); ++factor) {
      os << design.matrix[sample][factor];
      if (factor + 1 != design.factor_names.size())
        os << "\t";
    }
    os << "\n";
  }
  return os;
}

void
remove_factor(Design &design, size_t factor) {
  design.factor_names.erase(design.factor_names.begin() + factor);
  for (size_t sample = 0; sample < design.sample_names.size(); ++sample)
    design.matrix[sample].erase(design.matrix[sample].begin() + factor);
}

// Parses a natural number from its string representation. Throws exception if
// the string does not encode one.
static size_t
parse_natural_number(string encoding) {
  istringstream iss(encoding);
  size_t number;
  iss >> number;
  if (!iss)
    throw runtime_error("The token \"" +encoding + "\" "
                        "does not encode a natural number");
  return number;
}

istream&
operator>>(istream &table_encoding, SiteProportions &props) {
  props.chrom.clear();
  props.position = 0;
  props.strand.clear();
  props.context.clear();
  props.meth.clear();
  props.total.clear();

  string row_encoding;
  getline(table_encoding, row_encoding);

  // Skip lines contining only the newline character (e.g. the last line of the
  // proportion table).
  if(row_encoding.empty())
    return table_encoding;

  // Get the row name (which must be specified like this: "chr:position") and
  // parse it.
  istringstream row_stream(row_encoding);
  string row_name_encoding;
  row_stream >> row_name_encoding;

  // Every row must start an identifier consisiting of genomic loci of the
  // corresponding site. Here we check this identifier has the correct number
  // of colons.
  const size_t num_colon =
    count(row_name_encoding.begin(), row_name_encoding.end(), ':');

  if (num_colon != 3)
    throw runtime_error("Each row in the count table must start with "
                        "a line chromosome:position:strand:context."
                        "Got \"" + row_name_encoding + "\" instead." );

  // First parse the row identifier.
  istringstream name_stream(row_name_encoding);
  getline(name_stream, props.chrom, ':');

  if (props.chrom.empty())
    throw runtime_error("Error parsing " + row_name_encoding +
                        ": chromosome name is missing.");

  string position_encoding;

  getline(name_stream, position_encoding, ':');
  props.position = parse_natural_number(position_encoding);
  getline(name_stream, props.strand, ':');
  getline(name_stream, props.context, ':');

  // After parsing the row identifier, parse count proportions.
  size_t total_count, meth_count;

  while (row_stream >> total_count >> meth_count) {
    props.total.push_back(total_count);
    props.meth.push_back(meth_count);
  }

  if (!row_stream.eof())
    throw runtime_error("Some row entries are not natural numbers: " +
                        row_stream.str());

  if (props.total.size() != props.meth.size())
    throw runtime_error("This row does not encode proportions"
                        "correctly:\n" + row_encoding);
  return table_encoding;
}

static double
pi(Regression *reg, size_t sample, const gsl_vector *parameters) {
  double dot_prod = 0;
  // ADS: this function doesn't have a very helpful name

  for (size_t factor = 0; factor < reg->design.factor_names.size(); ++factor)
    dot_prod +=
      reg->design.matrix[sample][factor]*gsl_vector_get(parameters, factor);

  double p = exp(dot_prod)/(1 + exp(dot_prod));

  return p;
}

static double
neg_loglik(const gsl_vector *parameters, void *object) {
  Regression *reg = (Regression *)(object);
  const size_t num_parameters = reg->design.factor_names.size() + 1;

  double log_lik = 0;

  //dispersion parameter phi is the last element of parameter vector
  const double dispersion_param = gsl_vector_get(parameters,
                                                 num_parameters - 1);
  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t s = 0; s < reg->design.sample_names.size(); ++s) {
    const double n_s = reg->props.total[s];
    const double y_s = reg->props.meth[s];
    const double p_s = pi(reg, s, parameters);

    for(int k = 0; k < y_s; ++k) {
      log_lik += log((1 - phi)*p_s + phi*k);
    }

    for(int k = 0; k < n_s - y_s; ++k) {
      log_lik += log((1 - phi)*(1 - p_s) + phi*k);
    }

    for(int k = 0; k < n_s; ++k) {
      log_lik -= log(1 + phi*(k - 1));
    }
  }

  return (-1)*log_lik;
}

static void
neg_gradient(const gsl_vector *parameters, void *object,
             gsl_vector *output) {

  Regression *reg = (Regression *)(object);
  const size_t num_parameters = reg->design.factor_names.size() + 1;

  const double dispersion_param = gsl_vector_get(parameters,
                                                 num_parameters - 1);

  const double phi = exp(dispersion_param)/(1 + exp(dispersion_param));

  for(size_t f = 0; f < num_parameters; ++f) {

    double deriv = 0;

    for(size_t s = 0; s < reg->design.sample_names.size(); ++s) {
      int n_s = reg->props.total[s];
      int y_s = reg->props.meth[s];
      double p_s = pi(reg, s, parameters);

      double term = 0;

      //a parameter linked to p
      if(f < reg->design.factor_names.size()) {
        double factor = (1 - phi)*p_s*(1 - p_s)*reg->design.matrix[s][f];
        if (factor == 0) continue;

        for(int k = 0; k < y_s; ++k)
          term += 1/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term -= 1/((1 - phi)*(1 - p_s) + phi*k);

        deriv += term*factor;
      } else { // the parameter linked to phi
        for(int k = 0; k < y_s; ++k)
          term += (k - p_s)/((1 - phi)*p_s + phi*k);

        for(int k = 0; k < n_s - y_s; ++k)
          term += (k - (1 - p_s))/((1 - phi)*(1 - p_s) + phi*k);

        for(int k = 0; k < n_s; ++k) {
          term -= (k - 1)/(1 + phi*(k - 1));
        }

        deriv += term * phi * (1 - phi);
      }
    }

    gsl_vector_set(output, f, deriv);
  }

  gsl_vector_scale(output, -1.0);
}

static void
neg_loglik_and_grad(const gsl_vector *parameters,
                    void *object,
                    double *loglik_val,
                    gsl_vector *d_loglik_val) {

  *loglik_val = neg_loglik(parameters, object);
  neg_gradient(parameters, object, d_loglik_val);
}

bool
fit(Regression &r, vector<double> initial_parameters) {
  const size_t num_parameters = r.design.factor_names.size() + 1;

  if (initial_parameters.empty()) {
    for(size_t ind = 0; ind < num_parameters - 1; ++ind)
      initial_parameters.push_back(0.0);
    initial_parameters.push_back(-2.5);
  }

  if (initial_parameters.size() != num_parameters)
    throw runtime_error("Wrong number of initial parameters.");

  int status = 0;

  size_t iter = 0;

  gsl_multimin_function_fdf loglik_bundle;

  loglik_bundle.f = &neg_loglik;
  loglik_bundle.df = &neg_gradient;
  loglik_bundle.fdf = &neg_loglik_and_grad;
  loglik_bundle.n = num_parameters;
  loglik_bundle.params = (void *)&r;

  gsl_vector *parameters = gsl_vector_alloc(num_parameters);

  for (size_t parameter = 0; parameter < initial_parameters.size();
       ++parameter) {
    gsl_vector_set(parameters, parameter, initial_parameters[parameter]);
  }

  const gsl_multimin_fdfminimizer_type *T;

  //can also try gsl_multimin_fdfminimizer_conjugate_pr;
  T = gsl_multimin_fdfminimizer_conjugate_fr;

  gsl_multimin_fdfminimizer *s;
  s = gsl_multimin_fdfminimizer_alloc (T, num_parameters);

  gsl_multimin_fdfminimizer_set (s, &loglik_bundle, parameters, 0.001, 1e-4);

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    if (status)
      break;

    status = gsl_multimin_test_gradient (s->gradient, 1e-4);
  }
  while (status == GSL_CONTINUE && iter < 700);
  //It it reasonable to reduce the number of iterations to 500?
  // ADS: 700 vs. 500? what's the difference?

  r.max_loglik = (-1)*neg_loglik(s->x, &r);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(parameters);

  return status == GSL_SUCCESS;
}
bool
fit(Regression &r) {
  return fit(r, vector<double>());
}

/***************** RADMETH ALGORITHM *****************/
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

bool
has_low_coverage(const Regression &reg, const size_t test_factor) {

  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.design.matrix[sample][test_factor] == 1) {
      if (reg.props.total[sample] != 0)
        is_covered_in_test_factor_samples = true;
    }
    else {
      if (reg.props.total[sample] != 0)
        is_covered_in_other_samples = true;
    }
  }

  return !is_covered_in_test_factor_samples || !is_covered_in_other_samples;
}

bool
has_extreme_counts(const Regression &reg) {

  bool is_maximally_methylated = true;
  bool is_unmethylated = true;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.props.total[sample] != reg.props.meth[sample])
      is_maximally_methylated = false;

    if (reg.props.meth[sample] != 0)
      is_unmethylated = false;
  }

  return is_maximally_methylated || is_unmethylated;
}

/***********************************************************************
 * Run beta-binoimial regression using the specified table with
 * proportions and design matrix
 */
int
main(int argc, const char **argv) {

  try {

    static const string description =
      "calculate differential methylation scores";

    string outfile;
    string test_factor_name;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<design-matrix> <data-matrix>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("factor", 'f', "a factor to test", true, test_factor_name);

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

    ifstream design_file(design_filename);
    if (!design_file)
      throw runtime_error("could not open file: " + design_filename);

    ifstream table_file(table_filename);
    if (!table_file)
      throw runtime_error("could not open file: " + table_filename);

    ofstream of;
    if (!outfile.empty()) of.open(outfile);
    ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // initialize full design matrix from file
    Regression full_regression;
    design_file >> full_regression.design;

    // Check that provided test factor name exists and find its index.
    // Identify test factors with their indexes to simplify naming
    auto test_factor_it = find(begin(full_regression.design.factor_names),
                               end(full_regression.design.factor_names),
                               test_factor_name);

    if (test_factor_it == end(full_regression.design.factor_names))
      throw runtime_error("Error: " + test_factor_name +
                          " is not a part of the design specification.");

    const size_t test_factor = test_factor_it -
      begin(full_regression.design.factor_names);

    Regression null_regression;
    null_regression.design = full_regression.design;
    remove_factor(null_regression.design, test_factor);

    // Make sure that the first line of the proportion table file contains
    // names of the samples. Throw an exception if the names or their order
    // in the proportion table does not match those in the full design matrix.
    string sample_names_encoding;
    getline(table_file, sample_names_encoding);

    if (full_regression.design.sample_names != split(sample_names_encoding))
      throw runtime_error(sample_names_encoding + " does not match factor "
                          "names or their order in the design matrix. "
                          "Please verify that the design matrix and the "
                          "proportion table are correctly formatted.");

    // Performing the log-likelihood ratio test on proportions from each row
    // of the proportion table.
    while (table_file >> full_regression.props) {

      if (full_regression.design.sample_names.size() !=
          full_regression.props.total.size())
        throw runtime_error("found row with wrong number of columns");

      size_t coverage_factor = 0, coverage_rest = 0,
        meth_factor = 0, meth_rest = 0;

      for(size_t s = 0; s < full_regression.design.sample_names.size(); ++s) {
        if (full_regression.design.matrix[s][test_factor] != 0) {
          coverage_factor += full_regression.props.total[s];
          meth_factor += full_regression.props.meth[s];
        }
        else {
          coverage_rest += full_regression.props.total[s];
          meth_rest += full_regression.props.meth[s];
        }
      }

      out << full_regression.props.chrom << "\t"
          << full_regression.props.position << "\t"
          << full_regression.props.strand << "\t"
          << full_regression.props.context << "\t";

      // Do not perform the test if there's no coverage in either all case or
      // all control samples. Also do not test if the site is completely
      // methylated or completely unmethylated across all samples.
      if (has_low_coverage(full_regression, test_factor)) {
        out << -1;
      }
      else if (has_extreme_counts(full_regression)) {
        out << -1;
      }
      else {
        fit(full_regression);
        null_regression.props = full_regression.props;
        fit(null_regression);
        const double pval = loglikratio_test(null_regression.max_loglik,
                                             full_regression.max_loglik);

        // If error occured in fitting (p-val = nan or -nan).
        out << ((pval != pval) ? -1 : pval);
      }
      out << "\t" << coverage_factor << "\t" << meth_factor
          << "\t" << coverage_rest << "\t" << meth_rest << endl;
    }
  }
  catch (const std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
