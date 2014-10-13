/*    dmr2: computes DMRs based on HMRs and probability of differences
 *    at individual CpGs
 *
 *    Copyright (C) 2012 University of Southern California and
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
#include <fstream>
#include <algorithm>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::max;
using std::ifstream;

static void
read_diffs_file(const string &diffs_file, vector<GenomicRegion> &cpgs, const bool VERBOSE)
{
  std::ifstream in(diffs_file.c_str());
  string chrom, strand, seq;
  double diffscore;
  size_t pos, meth_a, unmeth_a, meth_b, unmeth_b;
  int n = 0;
  while (methpipe::read_methdiff_site(in, chrom, pos, strand, seq,
									  diffscore, meth_a, unmeth_a, meth_b, unmeth_b))
  {
	++n;
	cpgs.push_back(GenomicRegion(chrom, pos, pos + 1, seq, diffscore, strand[0]));
  }

  if (!in.eof() && !in.good())
	throw SMITHLABException("Immature termination when reading " + diffs_file +
							" around line " + smithlab::toa(n));
  if (VERBOSE)
	cerr << "Read " << n << " sites from " + diffs_file << endl;
}

static void
complement_regions(const size_t max_end, const vector<GenomicRegion> &a,
		   const size_t start, const size_t end,
		   vector<GenomicRegion> &cmpl) {
  cmpl.push_back(GenomicRegion(a[start]));
  cmpl.back().set_start(0);
  for (size_t i = start; i < end; ++i) {
    cmpl.back().set_end(a[i].get_start());
    cmpl.push_back(GenomicRegion(a[i]));
    cmpl.back().set_start(a[i].get_end());
  }
  cmpl.back().set_end(max_end);
}


static void
get_chrom_ends(const vector<GenomicRegion> &r, vector<size_t> &ends) {
  for (size_t i = 0; i < r.size() - 1; ++i)
    if (!r[i].same_chrom(r[i+1]))
      ends.push_back(i+1);
  ends.push_back(r.size());
}


static void
complement_regions(const size_t max_end,
      const vector<GenomicRegion> &r, vector<GenomicRegion> &r_cmpl) {

  vector<size_t> r_chroms;
  get_chrom_ends(r, r_chroms);
  size_t t = 0;
  for (size_t i = 0; i < r_chroms.size(); ++i) {
    complement_regions(max_end, r, t, r_chroms[i], r_cmpl);
    t = r_chroms[i];
  }
}


static bool
check_no_overlap(const vector<GenomicRegion> &regions) {
  for (size_t i = 1; i < regions.size(); ++i)
    if (regions[i].same_chrom(regions[i-1]) &&
	regions[i].get_start() < regions[i - 1].get_end())
      return false;
  return true;
}


static void
separate_sites(const vector<GenomicRegion> &dmrs,
        const vector<GenomicRegion> &sites,
        vector<pair<size_t, size_t> > &sep_sites) {
  const size_t n_dmrs = dmrs.size();

  for (size_t i = 0; i < n_dmrs; ++i) {
    GenomicRegion a(dmrs[i]);
    a.set_end(a.get_start() + 1);
    GenomicRegion b(dmrs[i]);
    b.set_start(b.get_end());
    b.set_end(b.get_end() + 1);

    vector<GenomicRegion>::const_iterator a_insert =
      lower_bound(sites.begin(), sites.end(), a);

    vector<GenomicRegion>::const_iterator b_insert =
      lower_bound(sites.begin(), sites.end(), b);

    sep_sites.push_back(std::make_pair(a_insert - sites.begin(),
            b_insert - sites.begin()));
  }
}


template <class T> bool
starts_before(const T &a, const T &b) {
  return (a.get_chrom() < b.get_chrom()) ||
    (a.same_chrom(b) && a.get_start() < b.get_start());
}

template <class T> bool
same_start(const T &a, const T &b) {
  return a.same_chrom(b) && a.get_start() == b.get_start();
}


static void
get_cpg_stats(const bool LOW_CUTOFF, const double sig_cutoff,
        const vector<GenomicRegion> &cpgs,
        const size_t start_idx, const size_t end_idx,
        size_t &total_cpgs, size_t &total_sig) {
  total_cpgs = end_idx - start_idx;
  for (size_t i = start_idx; i < end_idx; ++i) {
    if ((LOW_CUTOFF && (cpgs[i].get_score() < sig_cutoff)) ||
      (!LOW_CUTOFF && (cpgs[i].get_score() > 1.0 - sig_cutoff)))
      ++total_sig;
  }
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    double sig_cutoff = 0.05;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "computes DMRs based on "
        "HMRs and probability of differences at "
        "individual CpGs",
        "<methdiffs_1_gt_2> <hmr_1> <hmr_2> "
        "<dmr_1_lt_2> <dmr_2_lt_1>");
    opt_parse.add_opt("cutoff", 'c', "Significance cutoff (default: 0.05)",
        false, sig_cutoff);
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
    if (leftover_args.size() != 5) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string diffs_file = leftover_args[0];
    const string hmr1_file = leftover_args[1];
    const string hmr2_file = leftover_args[2];
    const string outfile_a = leftover_args[3];
    const string outfile_b = leftover_args[4];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[LOADING HMRS] " << hmr1_file << endl;

    vector<GenomicRegion> regions_a;
    ReadBEDFile(hmr1_file, regions_a);
    assert(check_sorted(regions_a));
    if (!check_sorted(regions_a))
      throw SMITHLABException("regions not sorted in file: " + hmr1_file);
    if (!check_no_overlap(regions_a))
      throw SMITHLABException("regions overlap in file: " + hmr1_file);

    if (VERBOSE)
      cerr << "[LOADING HMRS] " << hmr2_file << endl;

    vector<GenomicRegion> regions_b;
    ReadBEDFile(hmr2_file, regions_b);
    assert(check_sorted(regions_b));
    if (!check_sorted(regions_b))
      throw SMITHLABException("regions not sorted in file: " + hmr2_file);
    if (!check_no_overlap(regions_b))
      throw SMITHLABException("regions overlap in file: " + hmr2_file);

    if (VERBOSE)
      cerr << "[COMPUTING SYMMETRIC DIFFERENCE]" << endl;


    size_t max_end = 0;
    for (size_t i = 0; i < regions_a.size(); ++i)
      max_end = max(max_end, regions_a[i].get_end());
    for (size_t i = 0; i < regions_b.size(); ++i)
      max_end = max(max_end, regions_b[i].get_end());

    vector<GenomicRegion> a_cmpl, b_cmpl;
    complement_regions(max_end, regions_a, a_cmpl);
    complement_regions(max_end, regions_b, b_cmpl);

    vector<GenomicRegion> dmrs_a, dmrs_b;
    genomic_region_intersection_by_base(regions_a, b_cmpl, dmrs_a);
    genomic_region_intersection_by_base(regions_b, a_cmpl, dmrs_b);

    // separate the regions by chrom and by desert
    if (VERBOSE)
      cerr << "[READING CPG METH DIFFS]" << endl;
    vector<GenomicRegion> cpgs;
	read_diffs_file(diffs_file, cpgs, VERBOSE);
    if (!check_sorted(cpgs))
      throw SMITHLABException("CpGs not sorted in: " + diffs_file);
    if (VERBOSE)
      cerr << "[TOTAL CPGS]: " << cpgs.size() << endl;

    vector<pair<size_t, size_t> > sep_sites;
    separate_sites(dmrs_a, cpgs, sep_sites);

    for (size_t i = 0; i < dmrs_a.size(); ++i) {
      size_t total_cpgs = 0, total_sig = 0;
      get_cpg_stats(true, sig_cutoff,
      cpgs, sep_sites[i].first, sep_sites[i].second,
      total_cpgs, total_sig);
      dmrs_a[i].set_name(dmrs_a[i].get_name() + ":" + toa(total_cpgs));
      dmrs_a[i].set_score(total_sig);
    }

    sep_sites.clear();
    separate_sites(dmrs_b, cpgs, sep_sites);

    for (size_t i = 0; i < dmrs_b.size(); ++i) {
      size_t total_cpgs = 0, total_sig = 0;
      get_cpg_stats(false, sig_cutoff,
        cpgs, sep_sites[i].first, sep_sites[i].second,
        total_cpgs, total_sig);
      dmrs_b[i].set_name(dmrs_b[i].get_name() + ":" + toa(total_cpgs));
      dmrs_b[i].set_score(total_sig);
    }

    std::ofstream out_a(outfile_a.c_str());
    copy(dmrs_a.begin(), dmrs_a.end(),
    std::ostream_iterator<GenomicRegion>(out_a, "\n"));

    std::ofstream out_b(outfile_b.c_str());
    copy(dmrs_b.begin(), dmrs_b.end(),
    std::ostream_iterator<GenomicRegion>(out_b, "\n"));

    if (VERBOSE)
      cerr << "[OUTPUT FORMAT] COL4=NAME:N_COVERED_CPGS COL5=N_SIG_CPGS" << endl;
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
