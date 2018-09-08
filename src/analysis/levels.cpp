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

    while (in >> site) {

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
      else if (site.is_chh()) chh.update(site);
      else if (site.is_ccg()) ccg.update(site);
      else if (site.is_cxg()) cxg.update(site);
      else throw runtime_error("bad site context: " + site.context);

      cytosines.update(site);

      prev_site = site;
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << cytosines << endl
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
