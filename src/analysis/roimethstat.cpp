/*    roimethstat: compute average methylation in each of a set of regions
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <list>
#include <utility>
#include <tr1/cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::make_pair;
using std::tr1::round;


static pair<bool, bool>
meth_unmeth_calls(const size_t n_meth, const size_t n_unmeth) {
  static const double alpha = 0.95;
  // get info for binomial test
  double lower = 0.0, upper = 0.0;
  const size_t total = n_meth + n_unmeth;
  wilson_ci_for_binomial(alpha, total,
                         static_cast<double>(n_meth)/total, lower, upper);
  return make_pair(lower > 0.5, upper < 0.5);
}


/* Takes the sites and the regions. Creates two identical
   simplegenomicregions, a and b. Sets the end of a to the
   start + 1. Makes an iterator and
*/

static std::pair<size_t, size_t>
region_bounds(const vector<SimpleGenomicRegion> &sites,
              const GenomicRegion &region) {
  SimpleGenomicRegion a(region);
  a.set_end(a.get_start() + 1);

  // a_insert points to the first cpg site inside region
  vector<SimpleGenomicRegion>::const_iterator a_insert =
    lower_bound(sites.begin(), sites.end(), a);

  SimpleGenomicRegion b(region);
  b.set_start(b.get_end());
  b.set_end(b.get_end() + 1);

  // b_insert points to the first cpg site outside the region
  vector<SimpleGenomicRegion>::const_iterator b_insert =
    lower_bound(sites.begin(), sites.end(), b);

  return std::make_pair(a_insert - sites.begin(),
                        b_insert - sites.begin());

}

int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool PRINT_NAN = false;
    bool PRINT_ADDITIONAL_LEVELS = false;

    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("roimethstat", "Compute average CpG "
                           "methylation in each of a set of genomic intervals",
                           "<intervals-bed> <cpgs-bed>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("print-nan", 'P', "print all records (even if NaN score)",
                      false, PRINT_NAN);
    opt_parse.add_opt("more-levels", 'L', "print more meth level information",
                      false, PRINT_ADDITIONAL_LEVELS);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    const string regions_file = leftover_args.front();
    const string cpgs_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "FORMAT = NAME : CPGS : CPGS_WITH_READS : "
        "METH_READS : TOTAL_READS" << endl;

    vector<SimpleGenomicRegion> cpgs;
    vector<pair<double, double> > meths;
    vector<size_t> reads;

    const bool IS_METHPIPE = methpipe::is_methpipe_file_single(cpgs_file);
    if (VERBOSE)
      cerr << "meth file type: " << (IS_METHPIPE ? "methpipe" : "bed") << endl;

    if (IS_METHPIPE)
      methpipe::load_cpgs(cpgs_file, cpgs, meths, reads);
    else
      methpipe::load_cpgs_old(cpgs_file, cpgs, meths, reads);
    //not_methpipe_load_cpgs(cpgs_file, cpgs, meths, reads);

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    assert(check_sorted(regions));
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted in file: " + regions_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    for (size_t i = 0; i < regions.size(); ++i) {

      const std::pair<size_t, size_t> bounds(region_bounds(cpgs, regions[i]));
      //cout << bounds.first << "\t" << bounds.second << "\n";

      size_t meth = 0, read = 0;
      size_t cpgs_with_reads = 0;
      size_t called_total = 0, called_meth = 0;
      double mean_meth = 0.0;

      for (size_t j = bounds.first; j < bounds.second; ++j) {
        //cout << reads.size() << "\n";
        if (reads[j] > 0) {
          meth += static_cast<size_t>(meths[j].first);
          read += reads[j];
          ++cpgs_with_reads;

          const pair<bool, bool> calls =
            meth_unmeth_calls(meths[j].first, meths[j].second);
          called_total += (calls.first || calls.second);
          called_meth += calls.first;

          mean_meth += static_cast<double>(meths[j].first)/reads[j];
        }
      }

      const string name = regions[i].get_name() + ":" +
        toa(bounds.second - bounds.first) + ":" +
        toa(cpgs_with_reads) + ":" + toa(meth) + ":" + toa(read);
      regions[i].set_name(name);
      regions[i].set_score(static_cast<double>(meth)/read);
      if (PRINT_NAN || std::isfinite(regions[i].get_score())) {
        out << regions[i];
        if (PRINT_ADDITIONAL_LEVELS)
          out << '\t'
              << static_cast<double>(called_meth)/called_total << '\t'
              << mean_meth/cpgs_with_reads;
        out << endl;
      }
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
