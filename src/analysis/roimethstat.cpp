/*    roimethstat2: average methylation in each of a set of regions
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
#include <list>
#include <utility>

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
using std::ios_base;


static pair<bool, bool>
meth_unmeth_calls(const size_t n_meth, const size_t n_unmeth) {
  static const double alpha = 0.95;
  // get info for binomial test
  double lower = 0.0, upper = 0.0;
  const size_t total = n_meth + n_unmeth;
  wilson_ci_for_binomial(alpha, total,
                         static_cast<double>(n_meth)/total, lower, upper);
  return std::make_pair(lower > 0.5, upper < 0.5);
}



static std::pair<size_t, size_t>
region_bounds(const vector<SimpleGenomicRegion> &sites,
              const GenomicRegion &region) {
  SimpleGenomicRegion a(region);
  a.set_end(a.get_start() + 1);
  vector<SimpleGenomicRegion>::const_iterator a_insert =
    lower_bound(sites.begin(), sites.end(), a);

  SimpleGenomicRegion b(region);
  b.set_start(b.get_end());
  b.set_end(b.get_end() + 1);
  vector<SimpleGenomicRegion>::const_iterator b_insert =
    lower_bound(sites.begin(), sites.end(), b);

  return std::make_pair(a_insert - sites.begin(),
                        b_insert - sites.begin());
}



static void
not_methpipe_load_cpgs(const string &cpgs_file,
                       vector<SimpleGenomicRegion> &cpgs,
                       vector<pair<double, double> > &meths,
                       vector<size_t> &reads) {

  vector<GenomicRegion> cpgs_in;
  ReadBEDFile(cpgs_file, cpgs_in);
  assert(check_sorted(cpgs_in));
  if (!check_sorted(cpgs_in))
    throw SMITHLABException("regions not sorted in file: " + cpgs_file);

  for (size_t i = 0; i < cpgs_in.size(); ++i) {
    cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
    const string name(cpgs_in[i].get_name());
    const size_t total = atoi(smithlab::split(name, ":").back().c_str());
    const double meth_freq = cpgs_in[i].get_score();
    const size_t n_meth = roundf(meth_freq*total);
    const size_t n_unmeth = roundf((1.0 - meth_freq)*total);
    assert(n_meth + n_unmeth == total);
    reads.push_back(total);
    meths.push_back(std::make_pair(n_meth, n_unmeth));
  }
}

static void
process_with_cpgs_loaded(const bool METHPIPE_FORMAT,
                         const bool PRINT_NAN,
                         const bool PRINT_ADDITIONAL_LEVELS,
                         const string &cpgs_file,
                         vector<GenomicRegion> &regions,
                         std::ostream &out) {

  vector<SimpleGenomicRegion> cpgs;
  vector<pair<double, double> > meths;
  vector<size_t> reads;
  if (METHPIPE_FORMAT)
    methpipe::load_cpgs(cpgs_file, cpgs, meths, reads);
  else
    not_methpipe_load_cpgs(cpgs_file, cpgs, meths, reads);

  for (size_t i = 0; i < regions.size(); ++i) {

    const std::pair<size_t, size_t> bounds(region_bounds(cpgs, regions[i]));

    size_t meth = 0, read = 0;
    size_t cpgs_with_reads = 0;
    size_t called_total = 0, called_meth = 0;
    double mean_meth = 0.0;

    for (size_t j = bounds.first; j < bounds.second; ++j) {
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



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///
///  CODE BELOW HERE IS FOR SEARCHING ON DISK
///

static void
move_to_start_of_line(std::ifstream &in) {
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
find_start_line(const string &chr, const size_t idx, std::ifstream &cpg_in) {

  cpg_in.seekg(0, ios_base::beg);
  const size_t begin_pos = cpg_in.tellg();
  cpg_in.seekg(0, ios_base::end);
  const size_t end_pos = cpg_in.tellg();

  if (end_pos - begin_pos < 2)
    throw SMITHLABException("empty meth file");

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



static std::ifstream &
load_cpg(const bool METHPIPE_FORMAT, std::ifstream &cpg_in,
         GenomicRegion &cpg) {

  if (METHPIPE_FORMAT) {
    string chrom, seq, strand;
    size_t pos = 0, coverage = 0;
    double meth = 0.0;
    methpipe::read_site(cpg_in, chrom, pos, strand, seq, meth, coverage);
    cpg = GenomicRegion(chrom, pos, pos + 1, seq + ":" + toa(coverage),
                        meth, strand[0]);
  }
  else {
    cpg_in >> cpg;
  }

  return cpg_in;
}


static void
get_cpg_stats(const bool METHPIPE_FORMAT,
              std::ifstream &cpg_in, const GenomicRegion region,
              size_t &meth, size_t &reads, size_t &total_cpgs,
              size_t &cpgs_with_reads,
              size_t &called_total, size_t &called_meth,
              double &mean_meth) {

  string chrom(region.get_chrom());
  const size_t start_pos = region.get_start();
  const size_t end_pos = region.get_end();
  find_start_line(chrom, start_pos, cpg_in);

  GenomicRegion cpg;
  while (load_cpg(METHPIPE_FORMAT, cpg_in, cpg) &&
         (cpg.same_chrom(region) &&
          cpg.get_end() <= end_pos)) {
    if (start_pos <= cpg.get_start()) {
      ++total_cpgs;
      const size_t n_reads = atoi(smithlab::split(cpg.get_name(), ":").back().c_str());
      if (n_reads > 0) {
        const double meth_freq = cpg.get_score();
        const size_t n_meth = roundf(meth_freq*n_reads);
        const size_t n_unmeth = roundf((1.0 - meth_freq)*n_reads);
        meth += n_meth;
        reads += n_reads;
        ++cpgs_with_reads;

        const pair<bool, bool> calls = meth_unmeth_calls(n_meth, n_unmeth);
        called_total += (calls.first || calls.second);
        called_meth += calls.first;

        mean_meth += static_cast<double>(n_meth)/n_reads;
      }
    }
  }
  cerr << cpg << endl;
  cpg_in.clear();
}


static void
process_with_cpgs_on_disk(const bool METHPIPE_FORMAT,
                          const bool PRINT_NAN,
                          const bool PRINT_ADDITIONAL_LEVELS,
                          const string &cpgs_file,
                          vector<GenomicRegion> &regions,
                          std::ostream &out) {

  std::ifstream in(cpgs_file.c_str());
  for (size_t i = 0; i < regions.size() && in; ++i) {

    size_t meth = 0, read = 0;
    size_t cpgs_with_reads = 0;
    size_t called_total = 0, called_meth = 0;
    size_t total_cpgs = 0;
    double mean_meth = 0.0;

    get_cpg_stats(METHPIPE_FORMAT,
                  in, regions[i], meth, read, total_cpgs, cpgs_with_reads,
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
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



int
main(int argc, const char **argv) {

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
    assert(check_sorted(regions));
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted in file: " + regions_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    const bool METHPIPE_FORMAT =
      methpipe::is_methpipe_file_single(cpgs_file);

    if (VERBOSE)
      cerr << "CPG FILE FORMAT: "
           << (METHPIPE_FORMAT ? "METHPIPE" : "BED") << endl;

    if (LOAD_ENTIRE_FILE)
      process_with_cpgs_loaded(METHPIPE_FORMAT, PRINT_NAN,
                               PRINT_ADDITIONAL_LEVELS,
                               cpgs_file, regions, out);
    else
      process_with_cpgs_on_disk(METHPIPE_FORMAT, PRINT_NAN,
                                PRINT_ADDITIONAL_LEVELS,
                                cpgs_file, regions, out);
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
