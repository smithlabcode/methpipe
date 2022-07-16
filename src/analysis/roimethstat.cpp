/*    roimethstat: average methylation in each of a set of regions
 *
 *    Copyright (C) 2014 Andrew D. Smith
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
#include <utility>
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "zlib_wrapper.hpp"

#include "MethpipeSite.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::make_pair;
using std::ios_base;
using std::runtime_error;
using std::ifstream;

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


static std::pair<size_t, size_t>
region_bounds(const vector<MSite> &sites, const GenomicRegion &region) {
  const string chrom(region.get_chrom());
  const char strand(region.get_strand());

  const MSite a(chrom, region.get_start(), strand, "", 0, 0);
  vector<MSite>::const_iterator a_ins(lower_bound(begin(sites), end(sites), a));

  const MSite b(chrom, region.get_end(), strand, "", 0, 0);
  vector<MSite>::const_iterator b_ins(lower_bound(begin(sites), end(sites), b));

  return make_pair(a_ins - begin(sites), b_ins - begin(sites));
}


static void
process_with_cpgs_loaded(const bool VERBOSE,
                         const bool PRINT_NAN,
                         const bool PRINT_ADDITIONAL_LEVELS,
                         const string &cpgs_file,
                         vector<GenomicRegion> &regions,
                         std::ostream &out) {

  igzfstream in(cpgs_file);
  if (!in)
    throw runtime_error("cannot open file: " + cpgs_file);

  vector<MSite> cpgs;
  MSite the_cpg;
  while (in >> the_cpg)
    cpgs.push_back(the_cpg);

  if (VERBOSE)
    cerr << "[n_cpgs=" << cpgs.size() << "]" << endl;

  for (size_t i = 0; i < regions.size(); ++i) {

    const pair<size_t, size_t> bounds(region_bounds(cpgs, regions[i]));

    size_t total_meth = 0, total_reads = 0;
    size_t cpgs_with_reads = 0;
    size_t called_total = 0, called_meth = 0;
    double mean_meth = 0.0;

    for (size_t j = bounds.first; j < bounds.second; ++j) {
      if (cpgs[j].n_reads > 0) {
        total_meth += cpgs[j].n_meth();
        total_reads += cpgs[j].n_reads;
        ++cpgs_with_reads;

        auto calls = meth_unmeth_calls(cpgs[j].n_meth(), cpgs[j].n_unmeth());
        called_total += (calls.first || calls.second);
        called_meth += calls.first;

        mean_meth += cpgs[j].meth;
      }
    }

    const string name = regions[i].get_name() +
      ":" + toa(bounds.second - bounds.first) +
      ":" + toa(cpgs_with_reads) +
      ":" + toa(total_meth) +
      ":" + toa(total_reads);
    regions[i].set_name(name);

    regions[i].set_score(static_cast<double>(total_meth)/total_reads);
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


////////////////////////////////////////////////////////////////////////
///
///  CODE BELOW HERE IS FOR SEARCHING ON DISK
///

static void
move_to_start_of_line(ifstream &in) {
  char next;
  while (in.good() && in.get(next) && next != '\n') {
    in.unget();
    in.unget();
  }
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}

static void
find_start_line(const string &chr, const size_t idx, ifstream &cpg_in) {

  cpg_in.seekg(0, ios_base::beg);
  const size_t begin_pos = cpg_in.tellg();
  cpg_in.seekg(0, ios_base::end);
  const size_t end_pos = cpg_in.tellg();

  if (end_pos - begin_pos < 2)
    throw runtime_error("empty meth file");

  size_t step_size = (end_pos - begin_pos)/2;

  cpg_in.seekg(0, ios_base::beg);
  string low_chr;
  size_t low_idx = 0;
  cpg_in >> low_chr >> low_idx;

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  cpg_in.seekg(-2, ios_base::end);
  move_to_start_of_line(cpg_in);
  string high_chr;
  size_t high_idx;
  cpg_in >> high_chr >> high_idx;

  size_t pos = step_size;
  cpg_in.seekg(pos, ios_base::beg);
  move_to_start_of_line(cpg_in);

  while (step_size > 0) {
    string mid_chr;
    size_t mid_idx = 0;
    cpg_in >> mid_chr >> mid_idx;
    step_size /= 2;
    if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
      std::swap(mid_chr, high_chr);
      std::swap(mid_idx, high_idx);
      pos -= step_size;
    }
    else {
      std::swap(mid_chr, low_chr);
      std::swap(mid_idx, low_idx);
      pos += step_size;
    }
    cpg_in.seekg(pos, ios_base::beg);
    move_to_start_of_line(cpg_in);
  }
}


static bool
cpg_not_past_region(const GenomicRegion &region, const size_t end_pos,
                    const MSite &cpg) {
  return (cpg.chrom == region.get_chrom() && cpg.pos < end_pos) ||
    cpg.chrom < region.get_chrom();
}


static void
get_cpg_stats(ifstream &cpg_in, const GenomicRegion region,
              size_t &total_meth, size_t &total_reads,
              size_t &total_cpgs, size_t &cpgs_with_reads,
              size_t &called_total, size_t &called_meth,
              double &mean_meth) {

  const string chrom(region.get_chrom());
  const size_t start_pos = region.get_start();
  const size_t end_pos = region.get_end();
  find_start_line(chrom, start_pos, cpg_in);

  MSite cpg;
  while (cpg_in >> cpg && cpg_not_past_region(region, end_pos, cpg)) {
    if (start_pos <= cpg.pos && cpg.chrom == chrom) {
      ++total_cpgs;
      if (cpg.n_reads > 0) {

        total_meth += cpg.n_meth();
        total_reads += cpg.n_reads;
        ++cpgs_with_reads;

        auto calls = meth_unmeth_calls(cpg.n_meth(), cpg.n_unmeth());
        called_total += (calls.first || calls.second);
        called_meth += calls.first;

        mean_meth += cpg.meth;
      }
    }
  }
  cpg_in.clear();
}


static void
process_with_cpgs_on_disk(const bool PRINT_NAN,
                          const bool PRINT_ADDITIONAL_LEVELS,
                          const string &cpgs_file,
                          vector<GenomicRegion> &regions,
                          std::ostream &out) {

  ifstream in(cpgs_file);
  for (size_t i = 0; i < regions.size() && in; ++i) {

    size_t meth = 0, read = 0;
    size_t cpgs_with_reads = 0;
    size_t called_total = 0, called_meth = 0;
    size_t total_cpgs = 0;
    double mean_meth = 0.0;

    get_cpg_stats(in, regions[i], meth, read, total_cpgs, cpgs_with_reads,
                  called_total, called_meth, mean_meth);

    const string name = regions[i].get_name() + ":" +
      toa(total_cpgs) + ":" + toa(cpgs_with_reads) + ":" +
      toa(meth) + ":" + toa(read);
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
///
///  END OF CODE FOR SEARCHING ON DISK
///
////////////////////////////////////////////////////////////////////////

int
main_roimethstat(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool PRINT_NAN = false;
    bool LOAD_ENTIRE_FILE = false;
    bool PRINT_ADDITIONAL_LEVELS = false;

    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Compute average CpG "
                           "methylation in each of a set of genomic intervals",
                           "<intervals-bed> <cpgs-bed>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("print-nan", 'P', "print all records (even if NaN score)",
                      false, PRINT_NAN);
    opt_parse.add_opt("preload", 'L', "load all CpG sites",
                      false, LOAD_ENTIRE_FILE);
    opt_parse.add_opt("more-levels", 'M', "print more meth level information",
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

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw runtime_error("regions not sorted in file: " + regions_file);

    if (VERBOSE)
      cerr << "[n_regions=" << regions.size() << "]" << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    if (LOAD_ENTIRE_FILE)
      process_with_cpgs_loaded(VERBOSE, PRINT_NAN, PRINT_ADDITIONAL_LEVELS,
                               cpgs_file, regions, out);
    else
      process_with_cpgs_on_disk(PRINT_NAN, PRINT_ADDITIONAL_LEVELS,
                                cpgs_file, regions, out);
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
