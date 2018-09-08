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
 *    Copyright (C) 2014-2015 University of Southern California and
 *                            Andrew D. Smith and Benjamin E Decato
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
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeFiles.hpp"
#include "MethpipeSite.hpp"
#include "LevelsCounter.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::to_string;
using std::runtime_error;

// struct LevelsSite {
//   string chrom;
//   size_t pos;
//   string strand;
//   string context;
//   double meth;
//   size_t n_reads;

//   size_t n_meth() const {return std::round(meth*n_reads);}

//   void add(const LevelsSite &other) {
//     if (!is_mutated() && other.is_mutated())
//       context += 'x';
//     // ADS: order matters below as n_reads update invalidates n_meth()
//     // function until meth has been updated
//     const size_t total_c_reads = n_meth() + other.n_meth();
//     n_reads += other.n_reads;
//     meth = static_cast<double>(total_c_reads)/n_reads;
//   }

//   // ADS: function below has redundant check for is_cpg, which is
//   // expensive and might be ok to remove
//   bool is_mate_of(const LevelsSite &first) {
//     return (first.pos + 1 == pos && first.is_cpg() && is_cpg() &&
//             first.strand == "+" && strand == "-");
//   }
//   ////////////////////////////////////////////////////////////////////////
//   /////  Functions below test the type of site. These are CpG, CHH, and
//   /////  CHG divided into two kinds: CCG and CXG, the former including a
//   /////  CpG within. Also included is a function that tests if a site
//   /////  has a mutation.
//   /////  WARNING: None of these functions test for the length of their
//   /////  argument string, which could cause problems.
//   ////////////////////////////////////////////////////////////////////////
//   bool is_cpg() const {
//     return (context[0] == 'C' && context[1] == 'p' && context[2] == 'G');
//   }
//   bool is_chh() const {
//     return (context[0] == 'C' && context[1] == 'H' && context[2] == 'H');
//   }
//   bool is_ccg() const {
//     return (context[0] == 'C' && context[1] == 'C' && context[2] == 'G');
//   }
//   bool is_cxg() const {
//     return (context[0] == 'C' && context[1] == 'X' && context[2] == 'G');
//   }
//   bool is_mutated() const {
//     return context[3] == 'x';
//   }
// };


// struct LevelsCounter {
//   string context;
//   size_t total_sites;
//   size_t sites_covered;
//   size_t max_depth;
//   size_t mutations;
//   size_t total_c, total_t;
//   size_t called_meth, called_unmeth;
//   double mean_agg;
//   LevelsCounter(const string &c) :
//     context(c), total_sites(0), sites_covered(0), max_depth(0),
//     mutations(0), total_c(0), total_t(0),
//     called_meth(0), called_unmeth(0),
//     mean_agg(0.0) {}

//   void update(const LevelsSite &s) {
//     if (s.is_mutated()) {
//       ++mutations;
//     }
//     else if (s.n_reads > 0) {
//       ++sites_covered;
//       max_depth = std::max(max_depth, s.n_reads);
//       total_c += s.n_meth();
//       total_t += s.n_reads - s.n_meth();
//       mean_agg += s.meth;
//       double lower = 0.0, upper = 0.0;
//       wilson_ci_for_binomial(alpha, s.n_reads, s.meth, lower, upper);
//       called_meth += (lower > 0.5);
//       called_unmeth += (upper < 0.5);
//     }
//     ++total_sites;
//   }

//   size_t coverage() const {return total_c + total_t;}
//   size_t total_called() const {return called_meth + called_unmeth;}

//   double mean_meth_weighted() const {
//     return static_cast<double>(total_c)/coverage();
//   }
//   double fractional_meth() const {
//     return static_cast<double>(called_meth)/total_called();
//   }
//   double mean_meth() const {
//     return mean_agg/sites_covered;
//   }

//   string tostring() const {
//     static const string indent = string(2, ' ');
//     const bool good = (sites_covered != 0);
//     std::ostringstream oss;
//     // directly counted values
//     oss << context + ":\n"
//         << indent << "total_sites: " << total_sites << '\n'
//         << indent << "sites_covered: " << sites_covered << '\n'
//         << indent << "total_c: " << total_c << '\n'
//         << indent << "total_t: " << total_t << '\n'
//         << indent << "max_depth: " << max_depth << '\n'
//         << indent << "mutations: " << mutations << '\n'
//         << indent << "called_meth: " << called_meth << '\n'
//         << indent << "called_unmeth: " << called_unmeth << '\n'
//         << indent << "mean_agg: " << mean_agg << '\n';

//     // derived values
//     oss << indent << "coverage: " << coverage() << '\n'
//         << indent << "sites_covered_fraction: "
//         << static_cast<double>(sites_covered)/total_sites << '\n'
//         << indent << "mean_depth: "
//         << static_cast<double>(coverage())/total_sites << '\n'
//         << indent << "mean_depth_covered: "
//         << static_cast<double>(coverage())/sites_covered << '\n'
//         << indent << "mean_meth: "
//         << (good ? to_string(mean_meth()) : "NA")  << '\n'
//         << indent << "mean_meth_weighted: "
//         << (good ? to_string(mean_meth_weighted()) : "NA") << '\n'
//         << indent << "fractional_meth: "
//         << (good ? to_string(fractional_meth()) : "NA");
//     return oss.str();
//   }

//   static double alpha;
// };

// double LevelsCounter::alpha = 0.95;

// std::ostream &
// operator<<(std::ostream &out, const LevelsCounter &cs) {
//   return out << cs.tostring();
// }

// static void
// check_label(const string &observed, const string expected) {
//   if (observed != expected)
//     throw runtime_error("bad levels format [" + observed + "," + expected + "]");
// }

// std::istream &
// operator>>(std::istream &in, LevelsCounter &cs) {
//   in >> cs.context; // get the context
//   cs.context = cs.context.substr(0, cs.context.find_first_of(":"));

//   string label;
//   in >> label >> cs.total_sites; // the total sites
//   check_label(label, "total_sites:");

//   in >> label >> cs.sites_covered; // the sites covered
//   check_label(label, "sites_covered:");

//   in >> label >> cs.total_c; // the total c
//   check_label(label, "total_c:");

//   in >> label >> cs.total_t; // the total t
//   check_label(label, "total_t:");

//   in >> label >> cs.max_depth; // the max depth
//   check_label(label, "max_depth:");

//   in >> label >> cs.mutations; // the number of mutations
//   check_label(label, "mutations:");

//   in >> label >> cs.called_meth; // the number of sites called methylated
//   check_label(label, "called_meth:");

//   in >> label >> cs.called_unmeth; // the number of sites called unmethylated
//   check_label(label, "called_unmeth:");

//   in >> label >> cs.mean_agg; // the mean aggregate
//   check_label(label, "mean_agg:");

//   return in;
// }

// static std::istream &
// get_meth_unmeth(std::istream &in, LevelsSite &site) {
//   return methpipe::read_site(in, site.chrom, site.pos, site.strand,
//                              site.context, site.meth, site.n_reads);
// }

int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "compute methylation levels",
                           "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("alpha", 'a', "alpha for confidence interval",
                      false, LevelsCounter::alpha);
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
      throw std::runtime_error("bad input file: " + meth_file);

    LevelsCounter cpg("cpg");
    LevelsCounter cpg_symmetric("cpg_symmetric");
    LevelsCounter chh("chh");
    LevelsCounter cxg("cxg");
    LevelsCounter ccg("ccg");
    LevelsCounter cytosines("cytosines");

    MSite site, prev_site;
    size_t chrom_count = 0;

    while (in >> site) { //get_meth_unmeth(in, site)) {

      if (site.chrom != prev_site.chrom) {
        ++chrom_count;
        if (VERBOSE)
          cerr << "PROCESSING:\t" << site.chrom << "\n";
      }

      if (site.is_cpg()) {
        cpg.update(site);
        if (site.is_mate_of(prev_site)) {
          site.add(prev_site);
          cpg_symmetric.update(site);
        }
      }
      else if (site.is_chh())
        chh.update(site);
      else if (site.is_ccg())
        ccg.update(site);
      else if (site.is_cxg())
        cxg.update(site);
      else
        throw runtime_error("bad site context: " + site.context);

      cytosines.update(site);

      prev_site = site;
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << "number_of_chromosomes: " << chrom_count << endl
        << cytosines << endl
        << cpg << endl
        << cpg_symmetric << endl
        << chh << endl
        << ccg << endl
        << cxg << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
