/*    levels: a program to compute coverage statistics, mutation rates,
 *    and three different formulas for methylation levels described in
 *    the paper:
 *
 *        'Leveling' the playing field for analyses of single-base
 *         resolution DNA methylomes
 *         Schultz, Schmitz & Ecker (TIG 2012)
 *
 *    Note: the fractional methylation level calculated in this program
 *    is inspired but different from the paper. What we are doing here is
 *    using binomial test to determine significantly hyper/hypomethylated
 *    sites, and only use these subset of sites to calculate methylation
 *    level.
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
	       string &context, size_t &n_meth, size_t &n_unmeth,
               string &chr, size_t &pos, string &strand) {

  std::istringstream is(buffer);
  string name, dummy;
  double meth_freq = 0.0;
  is >> chr >> pos >> dummy >> name >> meth_freq >> strand;
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
        string &prev_chr, int &chr_count, size_t &coverage,
        size_t &pos, string &strand) {

  string chr;
  if (IS_METHPIPE_FILE) {
    double meth = 0.0;
    if (!methpipe::read_site(in, chr, pos, strand,
			     context, meth, coverage))
      return false;
    else {
      n_meth = std::tr1::round(meth*coverage);
      n_unmeth = std::tr1::round((1.0 - meth)*coverage);
      assert(n_meth + n_unmeth == coverage);
    }
  }
  else {
    string buffer;
    if (!getline(in, buffer))
      return false;
    parse_cpg_line(buffer, context, n_meth, n_unmeth, chr, pos, strand);
    coverage = n_meth + n_unmeth; 
  }
  if (chr != prev_chr) {
    ++chr_count;
    if (VERBOSE) {
      cerr << "PROCESSING:\t" << chr << "\n";
    }
  }
  prev_chr = chr;
  return true;
}

static bool
is_complementary_sites( size_t &position, string &strand,
    size_t &position_prev, string &strand_prev) {
  // pos difference is maximum 2, b/c we are looking at 2-mers and 3-mers
  return position - position_prev <= 2 && strand != strand_prev
    && strand == "-";
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
    size_t count_chh = 0, count_cxg = 0, count_ccg = 0;
    double mean_agg_cpg= 0, mean_agg_chh= 0, mean_agg_cxg= 0, mean_agg_ccg= 0;

    string prev_chr = "";
    int chr_count = 0;
    size_t n_meth = 0, n_unmeth = 0, coverage = 0, pos = 0;
    size_t n_meth_prev = 0, n_unmeth_prev = 0, coverage_prev = 0, pos_prev = 0;
    string context, strand;
    string context_prev, strand_prev;
    double level_prev = 0.0, level = 0.0;
    double lower = 0.0, upper = 0.0;
    double lower_prev = 0.0, upper_prev = 0.0;
    bool prev_mutated = false, last_is_complementary = true;
    while (get_meth_unmeth(IS_METHPIPE_FILE, VERBOSE,
          in, context, n_meth, n_unmeth, prev_chr, chr_count, coverage, pos, strand)) {
      // CpG is treated separately
      if (context.substr(0,3) == "CpG") {
        if (is_complementary_sites(pos, strand, pos_prev, strand_prev)) {
          // this site is the complementary site; it needs to be merged
          if (prev_mutated || *context.rbegin() == 'x') {
            // at least one of the two sites is mutated; ignore the whole site
            // note that mutated sites *are* counted in total_cpg_site,
            // but they are *not* counted in total_cpg_mapped
            ++cpg_mutations;
            ++total_cpg_sites;
            ++total_sites;
            total_mut += (*context.rbegin() == 'x');
            mapped_sites += (coverage>0);
            prev_mutated = false;
            last_is_complementary = true;
            continue;
          }
          else {
            // the complementary is not mutated, nor is the previous one
            // so the whole site is normal!
            if (coverage > max_cov)
                max_cov = coverage;
            ++total_cpg_sites;
            mapped_sites += (coverage>0);
            coverage = coverage + coverage_prev;
            if (coverage > 0) {
              ++total_cpg_mapped;
              n_meth = n_meth + n_meth_prev;
              n_unmeth = n_unmeth + n_unmeth_prev;
              cpg_cov += coverage;
              total_c_cpg += n_meth;
              total_t_cpg += n_unmeth;
              level = static_cast<double>(n_meth)/coverage;
              mean_agg_cpg += level;
              wilson_ci_for_binomial(alpha,
                  coverage, level, lower, upper);
              called_meth_cpg += (lower > 0.5);
              called_unmeth_cpg += (upper < 0.5);
              total_cov += coverage;
            }
          }
          last_is_complementary = true;
        }
        else {
          if (!last_is_complementary) {
            // it is likely that the previous CpG is isolated, and
            // it needs to be taken care of
            ++total_cpg_sites;
            if (*context_prev.rbegin() == 'x') {
              ++cpg_mutations;
              ++total_mut;
            }
            else {
              if (coverage_prev > 0) {
                ++total_cpg_mapped;
                ++mapped_sites;
                cpg_cov += coverage_prev;
                total_c_cpg += n_meth_prev;
                total_t_cpg += n_unmeth_prev;
                level_prev =
                  static_cast<double>(n_meth_prev)/coverage_prev;
                mean_agg_cpg += level_prev;
                wilson_ci_for_binomial(alpha,
                    coverage_prev, level_prev, lower_prev, upper_prev);
                called_meth_cpg += (lower_prev > 0.5);
                called_unmeth_cpg += (upper_prev < 0.5);
                total_cov += coverage_prev;
                if (coverage_prev > max_cov)
                    max_cov = coverage_prev;
              }
            }
          }
          else {
            // coupled symmetric CpG
            if (context == "CpGx") {
              // the first site is mutated; mark it and make the next
              // site also ignored
              prev_mutated = true;
              ++total_mut;
            }
            else {
              if (coverage > max_cov)
                  max_cov = coverage;
              prev_mutated = false;
            }
            mapped_sites += (coverage>0);
          }
          last_is_complementary = false;
        }
      }
      else {
        // Non-CpGs sites
        if (!last_is_complementary && context_prev == "CpG") {
          // just in case the input file is mixed with CpG pairs
          // and CpG sigulars, this one deals with sigulars
          if (*context_prev.rbegin() == 'x') {
            ++cpg_mutations;
            ++total_mut;
          }
          else {
            ++total_cpg_sites;
            if (coverage_prev > 0) {
              ++total_cpg_mapped;
              ++mapped_sites;
              cpg_cov += coverage_prev;
              total_c_cpg += n_meth_prev;
              total_t_cpg += n_unmeth_prev;
              level_prev =
                static_cast<double>(n_meth_prev)/coverage_prev;
              mean_agg_cpg += level_prev;
              wilson_ci_for_binomial(alpha,
                  coverage_prev, level_prev, lower_prev, upper_prev);
              called_meth_cpg += (lower_prev > 0.5);
              called_unmeth_cpg += (upper_prev < 0.5);
              total_cov += coverage_prev;
              if (coverage_prev > max_cov)
                  max_cov = coverage_prev;
              }
            }
          }
        if (n_meth + n_unmeth > 0) {
          // get info for mean methylation
          level = static_cast<double>(n_meth)/coverage;
          // get info for binomial test
          wilson_ci_for_binomial(alpha, coverage, level, lower, upper);
          if (context.substr(0,3) == "CXG") {
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

          ++mapped_sites;
        }
        // in C++11, this should be context.back(): this is
        // a hack until we move all of methpipe to C++11.
        if (*context.rbegin() == 'x')
          ++total_mut;

        last_is_complementary = true;
      }
      ++total_sites;
      // store previous site information for complementary CpGs
      n_meth_prev = n_meth;
      n_unmeth_prev = n_unmeth;
      coverage_prev = coverage;
      strand_prev = strand;
      pos_prev = pos;
      context_prev = context;
      lower_prev = lower;
      upper_prev = upper;
      level_prev = level;
    }

    // process the last line, in case it is a isolated CpG
    if (context_prev.substr(0,3) == "CpG" && !last_is_complementary) {
      // it is likely that previous CpG is isolated, and then
      // it needs to be taken care of
      if (context_prev == "CpGx") {
        ++cpg_mutations;
      }
      else {
        ++total_cpg_sites;
        if (coverage_prev > 0) {
          ++total_cpg_mapped;
          cpg_cov += coverage_prev;
          total_c_cpg += n_meth_prev;
          total_t_cpg += n_unmeth_prev;
          mean_agg_cpg += level_prev;
          called_meth_cpg += (lower_prev > 0.5);
          total_cov += coverage_prev;
          if (coverage_prev > max_cov)
              max_cov = coverage_prev;
        }
      }
    }

    const double weighted_mean_meth_cpg =
      static_cast<double>(total_c_cpg)/(total_c_cpg + total_t_cpg);
    const double fractional_meth_cpg =
      static_cast<double>(called_meth_cpg)/(called_meth_cpg
          + called_unmeth_cpg);
    const double mean_meth_cpg = mean_agg_cpg/total_cpg_mapped;

    const double weighted_mean_meth_chh =
      static_cast<double>(total_c_chh)/(total_c_chh + total_t_chh);
    const double fractional_meth_chh =
      static_cast<double>(called_meth_chh)/(called_meth_chh
          + called_unmeth_chh);
    const double mean_meth_chh = mean_agg_chh/count_chh;

    const double weighted_mean_meth_cxg =
      static_cast<double>(total_c_cxg)/(total_c_cxg + total_t_cxg);
    const double fractional_meth_cxg =
      static_cast<double>(called_meth_cxg)/(called_meth_cxg
          + called_unmeth_cxg);
    const double mean_meth_cxg = mean_agg_cxg/count_cxg;

    const double weighted_mean_meth_ccg =
      static_cast<double>(total_c_ccg)/(total_c_ccg + total_t_ccg);
    const double fractional_meth_ccg =
      static_cast<double>(called_meth_ccg)/(called_meth_ccg
          + called_unmeth_ccg);
    const double mean_meth_ccg = mean_agg_ccg/count_ccg; 

    const double weighted_mean_meth_all_c =
      static_cast<double>(total_c_cpg + total_c_chh
          + total_c_cxg + total_c_ccg)/
          (total_c_cpg + total_c_chh + total_c_cxg + total_c_ccg
          + total_t_cpg + total_t_chh + total_t_cxg + total_t_ccg);
    const double fractional_meth_all_c =
      static_cast<double>(called_meth_cpg + called_meth_chh
          + called_meth_cxg + called_meth_ccg)
          /(called_meth_cpg + called_meth_chh + called_meth_cxg + called_meth_ccg
          + called_unmeth_cpg + called_unmeth_chh
          + called_unmeth_cxg + called_unmeth_ccg);
    const double mean_meth_all_c = (mean_agg_cpg + mean_agg_chh
      + mean_agg_cxg + mean_agg_ccg)
      /(total_cpg_mapped + count_chh + count_cxg + count_ccg); 

    const double mean_coverage_cpg_all = total_cpg_sites > 0 ?
      static_cast<double>(cpg_cov)/total_cpg_sites : 0;
    const double mean_coverage_cpg_mapped = total_cpg_mapped > 0 ?
      static_cast<double>(cpg_cov)/total_cpg_mapped : 0;

    //const double cpg_mutation_rate = 
    //  static_cast<double>(cpg_mutations)/total_cpg_mapped;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << "NUMBER OF CHROMOSOMES:" << '\t' << chr_count << endl
        << "SITES:" << '\t' << total_sites << endl
        << "SITES COVERED:" << '\t' << mapped_sites << endl
        << "FRACTION COVERED:"
            << '\t' << static_cast<double>(mapped_sites)/total_sites << endl
        << "MEAN COVERAGE:"
            << '\t' << static_cast<double>(total_cov)/total_sites << endl
        << "MEAN COVERAGE (WHEN > 0):" << '\t' 
            << static_cast<double>(total_cov)/mapped_sites << endl
        << "MAX COVERAGE:" << '\t' << max_cov << endl
        << "SITES MUTATED:" << '\t' << total_mut << endl
        << "FRACTION MUTATED:" << '\t'
            << static_cast<double>(total_mut)/mapped_sites << endl
        << "SYMMETRICAL CpG SITES:" << '\t'
           << total_cpg_sites << endl
        << "SYMMETRICAL CpG SITES COVERED:" << '\t'
           << total_cpg_mapped << endl
        << "SYMMETRICAL CpG FRACTION COVERED:"
            << '\t' << static_cast<double>(total_cpg_mapped)
                  /total_cpg_sites << endl
        << "SYMMETRICAL CpG MEAN COVERAGE:" << '\t'
           << mean_coverage_cpg_all << endl
        << "SYMMETRICAL CpG MEAN COVERAGE (WHEN > 0):" << '\t'
           << mean_coverage_cpg_mapped << endl
        << "SYMMETRICAL CpG MUTATED:" << '\t'
           << cpg_mutations << endl
        << "SYMMETRICAL CpG FRACTION MUTATED:" << '\t'
            << static_cast<double>(cpg_mutations)/total_cpg_mapped << endl;
        // CpG coverage statistics are already minus mutations

    if (total_cpg_mapped != 0) {
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
    if (total_cpg_mapped + count_chh + count_cxg + count_ccg > 0) {
      out << "METHYLATION LEVELS (ALL CONTEXT):" << endl
          << '\t' <<  "mean_meth\t" << mean_meth_all_c << endl
          << '\t' << "w_mean_meth\t" << weighted_mean_meth_all_c << endl
          << '\t' << "frac_meth\t" << fractional_meth_all_c <<endl;
    } else {
      out << "METHYLATION LEVELS (ALL CONTEXT):" << endl
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
