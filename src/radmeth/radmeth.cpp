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
#include "GenomicRegion.hpp"

// Local headers.
#include "regression.hpp"
#include "combine_pvals.hpp"

using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::runtime_error;
using std::sort;
using std::min;

static bool
lt_locus_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.combined_pval < r2.combined_pval;
}

static bool
ls_locus_position(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.pos < r2.pos;
}

void
fdr(vector<PvalLocus> &loci) {

  sort(begin(loci), end(loci), lt_locus_pval);

  for (size_t idx = 0; idx < loci.size(); ++idx) {
    const double current_score = loci[idx].combined_pval;
    //Assign a new one.
    const double corrected_pval = loci.size()*current_score/(idx + 1);
    loci[idx].corrected_pval = corrected_pval;
  }

  for (vector<PvalLocus>::reverse_iterator it = loci.rbegin() + 1;
       it != loci.rend(); ++it) {

    const PvalLocus &prev_locus = *(it - 1);
    PvalLocus &cur_locus = *(it);

    cur_locus.corrected_pval =
      min(prev_locus.corrected_pval, cur_locus.corrected_pval);
  }

  for (auto it(begin(loci)); it != end(loci); ++it)
    it->corrected_pval = min(it->corrected_pval, 1.0);

  // restore original order
  sort(begin(loci), end(loci), ls_locus_position);
}

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
static int
run_regression(int argc, const char **argv) {

  const string command_name = argv[0];

  string outfile;
  string test_factor_name;
  bool VERBOSE = false;

  OptionParser opt_parse(command_name,
                         "calculate differential methylation scores",
                         "<design-matrix> <data-matrix>");
  opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                    false, outfile);
  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  opt_parse.add_opt("factor", 'f', "a factor to test", true, test_factor_name);

  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 2 || opt_parse.help_requested()) {
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

  std::ifstream design_file(design_filename);
  if (!design_file)
    throw runtime_error("could not open file: " + design_filename);

  std::ifstream table_file(table_filename);
  if (!table_file)
    throw runtime_error("could not open file: " + table_filename);

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  // initialize full design matrix from file
  Regression full_regression;
  design_file >> full_regression.design;

  // Check that provided test factor name exists and find its index.
  // Identify test factors with their indexes to simplify naming
  vector<string>::const_iterator test_factor_it =
    std::find(begin(full_regression.design.factor_names),
              end(full_regression.design.factor_names), test_factor_name);

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

  return EXIT_SUCCESS;
}

static int
run_adjust(int argc, const char **argv) {

  // first argument is name of command
  const string command_name = argv[0];

  string outfile;
  string bin_spec = "1:200:1";
  bool VERBOSE = false;

  /**************** GET COMMAND LINE ARGUMENTS *************************/
  OptionParser opt_parse(command_name,
                         "compute adjusted p-values using autocorrelation",
                         "<regression-output>");
  opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                    false, outfile);
  opt_parse.add_opt("bins", 'b', "corrlation bin specs", false , bin_spec);
  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 2 || opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl;
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
  if (leftover_args.size() != 1) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string bed_filename = leftover_args.front();
  /*********************************************************************/

  BinForDistance bin_for_dist(bin_spec);

  std::ifstream bed_file(bed_filename);
  if (!bed_file)
    throw "could not open file: " + bed_filename;

  if (VERBOSE)
    cerr << "[reading input]" << endl;

  // Read in all p-value loci. The loci that are not correspond to valid
  // p-values (i.e. values in [0, 1]) are skipped.
  vector<PvalLocus> pvals;
  std::string input_line, prev_chrom;

  size_t chrom_offset = 0;

  while (getline(bed_file, input_line)) {

    std::istringstream iss(input_line);
    std::string chrom, sign, name;
    size_t position;
    double pval;
    if (!(iss >> chrom >> position >> sign >> name >> pval))
      throw runtime_error("failed to parse line: " + input_line);

    // Skip loci that do not correspond to valid p-values.
    if (0 <= pval && pval <= 1) {
      // locus is on new chrom.
      if (!prev_chrom.empty() && prev_chrom != chrom)
        chrom_offset += pvals.back().pos;

      PvalLocus plocus;
      plocus.raw_pval = pval;
      plocus.pos = chrom_offset + bin_for_dist.max_dist() + 1 + position;

      pvals.push_back(plocus);
      prev_chrom = chrom;
    }
  }

  if (VERBOSE)
    cerr << "[combining p-values]" << endl;
  combine_pvals(pvals, bin_for_dist);

  if (VERBOSE)
    cerr << "[running multiple test adjustment]" << endl;
  fdr(pvals);

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  std::ifstream original_bed_file(bed_filename);

  update_pval_loci(original_bed_file, pvals, out);

  //TODO: Check that the regions do not overlap & sorted
  return EXIT_SUCCESS;
}

// Attemps to find the next significant CpG site. Returns true if one was found
// and flase otherwise.
static bool
read_next_significant_cpg(istream &cpg_stream, GenomicRegion &cpg,
                          double cutoff, bool &skipped_any, bool &n_sig_sites,
                          size_t &test_cov, size_t &test_meth,
                          size_t &rest_cov, size_t &rest_meth) {
  GenomicRegion region;
  skipped_any = false;
  n_sig_sites = false;
  string cpg_encoding;

  while (getline(cpg_stream, cpg_encoding)) {
    string record, chrom, name, sign;
    size_t position;
    double raw_pval, adjusted_pval, corrected_pval;

    std::istringstream iss(cpg_encoding);
    iss.exceptions(std::ios::failbit);
    iss >> chrom >> position >> sign >> name >> raw_pval
        >> adjusted_pval >> corrected_pval
        >> test_cov >> test_meth >> rest_cov >> rest_meth;

    if (0 <= corrected_pval && corrected_pval < cutoff) {
      cpg.set_chrom(chrom);
      cpg.set_start(position);
      cpg.set_end(position + 1);
      n_sig_sites = (0 <= raw_pval && raw_pval < cutoff);
      return true;
    }
    skipped_any = true;
  }

  return false;
}

static void
merge(istream &cpg_stream, ostream &dmr_stream, double cutoff) {

  GenomicRegion dmr;
  dmr.set_name("dmr");

  size_t dmr_test_cov = 0;
  size_t dmr_test_meth = 0;
  size_t dmr_rest_cov = 0;
  size_t dmr_rest_meth = 0;

  size_t test_cov = 0;
  size_t test_meth = 0;
  size_t rest_cov = 0;
  size_t rest_meth = 0;

  // Find the first significant CpG, or terminate the function if none exist.
  bool skipped_last_cpg, n_sig_sites;
  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg,
                                 n_sig_sites, test_cov, test_meth,
                                 rest_cov, rest_meth))
    return;

  dmr.set_score(n_sig_sites);
  dmr_test_cov += test_cov;
  dmr_test_meth += test_meth;
  dmr_rest_cov += rest_cov;
  dmr_rest_meth += rest_meth;

  GenomicRegion cpg;
  cpg.set_name("dmr");

  while (read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg,
                                   n_sig_sites, test_cov, test_meth,
                                   rest_cov, rest_meth)) {

    if (skipped_last_cpg || cpg.get_chrom() != dmr.get_chrom()) {
      if (dmr.get_score() != 0)
        dmr_stream << dmr.get_chrom() << '\t'
                   << dmr.get_start() << '\t'
                   << dmr.get_end()   << '\t'
                   << dmr.get_name()  << '\t'
                   << dmr.get_score() << '\t'
                   << double(dmr_test_meth)/dmr_test_cov -
          double(dmr_rest_meth)/dmr_rest_cov << endl;
      dmr = cpg;
      dmr.set_score(n_sig_sites);
      dmr_test_cov = test_cov;
      dmr_test_meth = test_meth;
      dmr_rest_cov = rest_cov;
      dmr_rest_meth = rest_meth;
    }
    else {
      dmr.set_end(cpg.get_end());
      dmr.set_score(dmr.get_score() + n_sig_sites);
      dmr_test_cov += test_cov;
      dmr_test_meth += test_meth;
      dmr_rest_cov += rest_cov;
      dmr_rest_meth += rest_meth;
    }
  }
  if (dmr.get_score() != 0) {
    dmr_stream << dmr.get_chrom() << '\t'
               << dmr.get_start() << '\t'
               << dmr.get_end()   << '\t'
               << dmr.get_name()  << '\t'
               << dmr.get_score() << '\t'
               << double(dmr_test_meth)/dmr_test_cov -
      double(dmr_rest_meth)/dmr_rest_cov << endl;
  }
}

static int
run_merge(int argc, const char **argv) {

  // first argument is name of command
  const string command_name = argv[0];

  /* FILES */
  string outfile;
  string bin_spec = "1:200:25";
  double cutoff = 0.01;

  /**************** GET COMMAND LINE ARGUMENTS *************************/
  OptionParser opt_parse(command_name,
                         "merge significantly differentially"
                         " methylated CpGs into DMRs",
                         "<bed-file-in-radmeth-format>");
  opt_parse.add_opt("output", 'o',
                    "output file (default: stdout)", false, outfile);
  opt_parse.add_opt("cutoff", 'p', "P-value cutoff (default: 0.01)",
                    false , cutoff);
  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 1 || opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl;
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
  if (leftover_args.size() != 1) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string bed_filename = leftover_args.front();
  /************************************************************************/

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  std::ifstream in(bed_filename);
  if (!in)
    throw runtime_error("could not open file: " + bed_filename);

  merge(in, out, cutoff);

  return EXIT_SUCCESS;
}


int
main(int argc, const char **argv) {

  try {
    const string prog_name = strip_path(argv[0]);

    const string main_help_message =
      "Uage: " + prog_name + " [COMMAND] [PARAMETERS]\n\n"
      "Available commands: \n"
      "  regression  Calculates multi-factor differential methylation scores.\n"
      "  adjust      Adjusts the p-value of each site based on the p-value of "
      "its neighbors.\n"
      "  merge       Combines significantly differentially methylated CpGs into"
      " DMRs.\n";

    if (argc == 1) {
      cerr << "Analysis of differential methylation in"
           << " multi-factor bisulfite sequencing experiments" << endl
           << main_help_message;
      return EXIT_SUCCESS;
    }

    // first argument is name of command
    const string command_name = argv[1];

    if (command_name == "regression")
      return run_regression(argc - 1, argv + 1);
    else if (command_name == "adjust")
      return run_adjust(argc - 1, argv + 1);
    else if (command_name == "merge")
      return run_merge(argc - 1, argv + 1);
    else {
      cerr << "ERROR: invalid command name: \""
           << command_name << "\"" << endl
           << main_help_message;
      return EXIT_SUCCESS;
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
