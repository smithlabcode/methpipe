/*    levels: a program to compute coverage statistics, mutation rates,
 *    and three different formulas for methylation levels described in
 *    the paper:
 *
 *        'Leveling' the playing field for analyses of single-base
 *         resolution DNA methylomes
 *         Schultz, Schmitz & Ecker (TIG 2012)
 *
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Benjamin E Decato
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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <tr1/cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeFiles.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;



static void
parse_cpg_line(const string &buffer,
	       string &context, size_t &n_meth, size_t &n_unmeth) {

  std::istringstream is(buffer);
  string name, dummy;
  double meth_freq = 0.0;
  is >> dummy >> dummy >> dummy >> name >> meth_freq;
  const size_t sep_pos = name.find_first_of(":");
  const size_t total = atoi(name.substr(sep_pos + 1).c_str());
  context = name.substr(0, sep_pos);
  n_meth = std::tr1::round(meth_freq*total);
  n_unmeth = std::tr1::round((1.0 - meth_freq)*total);
  assert(n_meth + n_unmeth == total);
}



static bool
get_meth_unmeth(const bool IS_METHPIPE_FILE, const bool VERBOSE,
        std::ifstream &in, string &context, size_t &n_meth, size_t &n_unmeth,
        string &prev_chr, size_t &coverage) {

  if (IS_METHPIPE_FILE) {
    string dummy;
    string chr;
    size_t dummy_pos = 0;
    double meth = 0.0;
    if (!methpipe::read_site(in, chr, dummy_pos, dummy,
			     context, meth, coverage))
      return false;
    else {
      if (chr != prev_chr && VERBOSE) {
          cerr << "PROCESSING:\t" << chr << "\n";
          prev_chr = chr;
      }
      n_meth = std::tr1::round(meth*coverage);
      n_unmeth = std::tr1::round((1.0 - meth)*coverage);
      assert(n_meth + n_unmeth == coverage);
    }
  }
  else {
    string buffer;
    if (!getline(in, buffer))
      return false;
    parse_cpg_line(buffer, context, n_meth, n_unmeth);
  }
  return true;
}



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool IS_METHPIPE_FILE = true;
    string outfile;
    double alpha = 0.95;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "compute methylation levels",
			   "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
		      false, outfile);
    opt_parse.add_opt("alpha", 'a', "alpha for confidence interval",
		      false, alpha);
    opt_parse.add_opt("bed", 'b', "file in bed format",
		      false, IS_METHPIPE_FILE);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string meth_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in(meth_file.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + meth_file);

    size_t mapped_sites = 0, total_sites = 0, total_mut = 0;
    size_t total_cov = 0, max_cov = 0, cpg_cov = 0;
    size_t total_cpg_sites = 0, total_cpg_mapped = 0, cpg_mutations = 0;
    size_t total_c_cpg = 0, total_c_chh = 0, total_c_cxg = 0, total_c_ccg = 0;
    size_t total_t_cpg = 0, total_t_chh = 0, total_t_cxg = 0, total_t_ccg = 0;
    size_t called_meth_cpg = 0, called_meth_chh = 0;
    size_t called_meth_cxg = 0, called_meth_ccg = 0;
    size_t called_unmeth_cpg = 0, called_unmeth_chh = 0;
    size_t called_unmeth_cxg = 0, called_unmeth_ccg = 0;
    size_t count_cpg = 0, count_chh = 0, count_cxg = 0, count_ccg = 0;
    double mean_agg_cpg= 0, mean_agg_chh= 0, mean_agg_cxg= 0, mean_agg_ccg= 0;

    string buffer;
    string prev_chr = "";
    size_t n_meth = 0, n_unmeth = 0, coverage = 0;
    string context;
    while (get_meth_unmeth(IS_METHPIPE_FILE, VERBOSE,
			   in, context, n_meth, n_unmeth, prev_chr, coverage)) {
      if(context.substr(0,3)=="CpG") {
        ++total_cpg_sites;
        if(context=="CpGx")
          ++cpg_mutations;
      }
      if (n_meth + n_unmeth > 0) {
        // get info for mean methylation
        const size_t N = n_meth + n_unmeth;
        const double level = static_cast<double>(n_meth)/N;

        // get info for binomial test
        double lower = 0.0, upper = 0.0;
        wilson_ci_for_binomial(alpha, N, level, lower, upper);

        if(context.substr(0,3) == "CpG") {
          total_c_cpg += n_meth;
          total_t_cpg += n_unmeth;
          mean_agg_cpg += level;
          ++count_cpg;
          called_meth_cpg += (lower > 0.5);
          called_unmeth_cpg += (upper < 0.5);
          cpg_cov += coverage;
          ++total_cpg_mapped;
        }
        else if (context.substr(0,3) == "CXG") {
          total_c_cxg += n_meth;
          total_t_cxg+= n_unmeth;
          mean_agg_cxg += level;
          ++count_cxg;
          called_meth_cxg += (lower > 0.5);
          called_unmeth_cxg += (upper < 0.5);
        }
        else if (context.substr(0,3) == "CHH") {
          total_c_chh += n_meth;
          total_t_chh += n_unmeth;
          mean_agg_chh += level;
          ++count_chh;
          called_meth_chh += (lower > 0.5);
          called_unmeth_chh += (upper < 0.5);
        }
        else if (context.substr(0,3) == "CCG") {
          total_c_ccg += n_meth;
          total_t_ccg += n_unmeth;
          mean_agg_ccg += level;
          ++count_ccg;
          called_meth_ccg += (lower > 0.5);
          called_unmeth_ccg += (upper < 0.5);
        }
        else {
          throw SMITHLABException("bad context in input file: " + context);
        }

        total_cov += coverage;
        if (coverage > max_cov)
            max_cov = coverage;

        // in C++11, this should be context.back(): this is
        // a hack until we move all of methpipe to C++11.
        if (*context.rbegin() == 'x')
          ++total_mut;

        ++mapped_sites;
      }
      ++total_sites;
    }

    const double weighted_mean_meth_cpg =
      static_cast<double>(total_c_cpg)/(total_c_cpg + total_t_cpg);

    const double fractional_meth_cpg =
      static_cast<double>(called_meth_cpg)/
                         (called_meth_cpg + called_unmeth_cpg);

    const double mean_meth_cpg = mean_agg_cpg/count_cpg;

    const double weighted_mean_meth_chh =
      static_cast<double>(total_c_chh)/(total_c_chh + total_t_chh);

    const double fractional_meth_chh =
      static_cast<double>(called_meth_chh)/
                         (called_meth_chh + called_unmeth_chh);

    const double mean_meth_chh = mean_agg_chh/count_chh;

    const double weighted_mean_meth_cxg =
      static_cast<double>(total_c_cxg)/(total_c_cxg + total_t_cxg);

    const double fractional_meth_cxg =
      static_cast<double>(called_meth_cxg)/
                         (called_meth_cxg + called_unmeth_cxg);

    const double mean_meth_cxg = mean_agg_cxg/count_cxg;

    const double weighted_mean_meth_ccg =
      static_cast<double>(total_c_ccg)/(total_c_ccg + total_t_ccg);

    const double fractional_meth_ccg =
      static_cast<double>(called_meth_ccg)/
                         (called_meth_ccg + called_unmeth_ccg);

    const double mean_meth_ccg = mean_agg_ccg/count_ccg; 

    const double cpg_mutation_rate = 
      static_cast<double>(cpg_mutations)/total_cpg_mapped;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << "SITES:" << '\t' << total_sites << endl
        << "SITES COVERED:" << '\t' << mapped_sites << endl
        << "FRACTION MUTATED:" << '\t'
            << static_cast<double>(total_mut)/total_sites << endl
        << "FRACTION COVERED:"
            << '\t' << static_cast<double>(mapped_sites)/total_sites << endl
        << "MAX COVERAGE:" << '\t' << max_cov << endl
        << "SYMMETRICAL CpG COVERAGE:" << '\t'
           << static_cast<double>(2.0*cpg_cov)/total_cpg_sites << endl
        << "SYMMETRICAL CpG COVERAGE (WHEN > 0):" << '\t'
           << static_cast<double>(2.0*cpg_cov)/total_cpg_mapped << endl
        << "SYMMETRICAL CpG COVERAGE (minus mutations):" << '\t'
           << static_cast<double>((2.0-cpg_mutation_rate)*cpg_cov)
                                    /total_cpg_sites << endl
        << "SYMMETRICAL CpG COVERAGE (WHEN > 0) (minus mutations):" << '\t'
           << static_cast<double>((2.0-cpg_mutation_rate)*cpg_cov)
                                   /total_cpg_mapped << endl
        << "MEAN COVERAGE:"
            << '\t' << static_cast<double>(total_cov)/total_sites << endl
        << "MEAN COVERAGE (WHEN > 0):" << '\t' 
            << static_cast<double>(total_cov)/mapped_sites << endl;

    if (count_cpg != 0) {
      out << "METHYLATION LEVELS (CpG CONTEXT):" << endl
          << '\t' <<  "mean_meth" << '\t' << mean_meth_cpg << endl
          << '\t' << "w_mean_meth" << '\t' << weighted_mean_meth_cpg << endl
          << '\t' << "frac_meth" << '\t' << fractional_meth_cpg << endl;
    }
    else {
      out << "METHYLATION LEVELS (CpG CONTEXT):" << endl
          << '\t' <<  "mean_meth\tN/A" << endl
          << '\t' << "w_mean_meth\tN/A" << endl
          << '\t' << "frac_meth\tN/A" << endl;
    }
    if (count_chh != 0) {
      out << "METHYLATION LEVELS (CHH CONTEXT):" << endl
          << '\t' <<  "mean_meth" << '\t' << mean_meth_chh << endl
          << '\t' << "w_mean_meth" << '\t' << weighted_mean_meth_chh << endl
          << '\t' << "frac_meth" << '\t' << fractional_meth_chh  << endl;
    }
    else {
      out << "METHYLATION LEVELS (CHH CONTEXT):" << endl
          << '\t' <<  "mean_meth\tN/A" << endl
          << '\t' << "w_mean_meth\tN/A" << endl
          << '\t' << "frac_meth\tN/A" << endl;
    }
    if (count_cxg != 0) {
      out << "METHYLATION LEVELS (CXG CONTEXT):" << endl
          << '\t' <<  "mean_meth" << '\t' << mean_meth_cxg << endl
          << '\t' << "w_mean_meth" << '\t' << weighted_mean_meth_cxg << endl
          << '\t' << "frac_meth" << '\t' << fractional_meth_cxg << endl;
    } else {
      out << "METHYLATION LEVELS (CXG CONTEXT):" << endl
          << '\t' <<  "mean_meth\tN/A" << endl
          << '\t' << "w_mean_meth\tN/A" << endl
          << '\t' << "frac_meth\tN/A" << endl;
    }
    if (count_ccg != 0) {
      out << "METHYLATION LEVELS (CCG CONTEXT):" << endl
          << '\t' <<  "mean_meth" << '\t' << mean_meth_ccg << endl
          << '\t' << "w_mean_meth" << '\t' << weighted_mean_meth_ccg << endl
          << '\t' << "frac_meth" << '\t' << fractional_meth_ccg << endl;
    } else {
      out << "METHYLATION LEVELS (CCG CONTEXT):" << endl
          << '\t' <<  "mean_meth\tN/A" << endl
          << '\t' << "w_mean_meth\tN/A" << endl
          << '\t' << "frac_meth\tN/A" << endl;
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
