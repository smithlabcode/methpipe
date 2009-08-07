/*    cpgcaller: a program for calling CpG methylation status from
 *    bisulfite capture sequencing with Solexa reads
 *
 *    Copyright (C) 2009 University of Southern California and
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

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "BSUtils.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;


struct meth_base {
  meth_base() : unconv(0.0), conv(0.0) {}
  void inc(char c) {
    if (is_cytosine(c)) unconv++;
    if (is_thymine(c)) conv++;
  }
  double unconv;
  double conv;
};

static size_t
get_min_obs_for_confidence(const double d, const double alpha) {
  size_t i = 0;
  double curr_width = 1.0;
  while (curr_width > d) {
    ++i;
    double lower = 0, upper = 1;
    wilson_ci_for_binomial(alpha, i, 0.5, lower, upper);
    curr_width = upper - lower;
  }
  return i;
}

size_t
get_meth_state(const double alpha, const double critical_value,
	       const double interval_width, 
	       const size_t n, const double p, 
	       bool &confident_call, double &upper, double &lower) {

  static const size_t NO_CALL = 3ul;
  static const size_t PARTIAL_METH = 2ul;
  static const size_t METHYLATED = 1ul;
  static const size_t UNMETHYLATED = 0ul;
  
  size_t meth_state = NO_CALL;

  if (n > 0) {
    wilson_ci_for_binomial(alpha, n, p, lower, upper);
    const double tail_size = (1.0 - critical_value);
    confident_call = ((upper - lower) < interval_width);
    
    assert(upper <= 1.0 && lower >= 0.0);
    
    if (confident_call) meth_state = PARTIAL_METH;
    
    if (upper < tail_size) meth_state = UNMETHYLATED;
    else if (lower > critical_value)
      meth_state = METHYLATED;
  
  }
  return meth_state;
}


string
call_cpg_state(const double critical_value,
	       const double alpha,
	       const size_t min_obs_for_confidence,
	       const double interval_width,
	       
	       const vector<GenomicRegion> &reads, 
	       const vector<string> &rd_seqs,
	       const GenomicRegion &region, 
	       const string &seq) {
  
  // 1) COLLECT INFORMATION ABOUT EACH SITE

  assert(reads.size() == rd_seqs.size());
  assert(region.get_width() == seq.length());

  const size_t region_start = region.get_start();
  const size_t region_size = region.get_width();

  std::vector<meth_base> columns_pos;
  std::vector<meth_base> columns_neg;

  columns_pos.resize(region_size);
  columns_neg.resize(region_size);

  const string seq_rc(revcomp(seq));
  
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t read_width = reads[i].get_width();
    if (reads[i].pos_strand()) {
      assert(reads[i].get_start() >= region_start);
      const size_t offset = reads[i].get_start() - region_start;
      for (size_t j = 0; j < read_width; ++j)
	columns_pos[offset + j].inc(rd_seqs[i][j]);
    }
    else {
      const size_t offset = region_size - read_width -
	(reads[i].get_start() - region_start);
      for (size_t j = 0; j < read_width; ++j)
	columns_neg[offset + j].inc(rd_seqs[i][j]);
    }
  }
  reverse(columns_neg.begin(), columns_neg.end());

  // 2) DO THE STATISTICS AND FORMAT THE OUTPUT
  
  std::ostringstream ss;
  bool first_print = true;

  const string chrom(region.get_chrom());
  const string name(region.get_name());
  const size_t start = region.get_start();
  for (size_t i = 0; i < columns_pos.size() - 1; ++i) {
    if (is_cpg(seq, i)) {
      const size_t meth = (columns_pos[i].unconv + columns_neg[i + 1].unconv);
      const size_t unmeth = (columns_pos[i].conv + columns_neg[i + 1].conv);
      const size_t n = meth + unmeth;
      const double p = meth/std::max(1.0, static_cast<double>(n));
      
      bool confident_call = false;
      double upper = 0, lower = 0;
      const size_t meth_state = 
	get_meth_state(alpha, critical_value, interval_width,
		       n, p, confident_call, upper, lower);
      if (!first_print) 
	ss << endl;
      else first_print = false;
      
      const bool comparable_confident = (n >= min_obs_for_confidence);
      ss << chrom << "\t" << start + i << "\t" << start + i + 1 << "\t"
	 << name << ":" << unmeth << ":" << meth << ":" 
	 << comparable_confident << ":" << confident_call << "\t" 
	 << meth_state << "\t+";
    }
  }
  return ss.str();
}


int 
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    
    string mapped_locations_file;
    string chrom_dir;

    string outfile;
    string regions_file;
    
    double alpha = 0.1;
    double crit = 0.75;
    double interval_width = 0.25;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("cpgcaller", "A program for calling cpg methylation "
			   "status from a Solexa bisulfite capture experiment");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "dir of chrom sequences (single-line format)",
		      true , chrom_dir);
    opt_parse.add_opt("crit", 'C', "critical value", 
		      false , crit);
    opt_parse.add_opt("width", 'w', "interval width for confident call", 
		      false , interval_width);
    opt_parse.add_opt("alpha", 'a', "alpha value", false , alpha);
    opt_parse.add_opt("mapped", 'm', "file of mapped locations", 
		      true , mapped_locations_file);
    opt_parse.add_opt("regions", 'r', "file of target regions", 
		      true , regions_file);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "reading mapped locations" << endl;
    vector<GenomicRegion> mapped_locations;
    ReadBEDFile(mapped_locations_file, mapped_locations);
    assert(check_sorted(mapped_locations));
    if (VERBOSE)
      cerr << "read " << mapped_locations.size() << " total locations" << endl;
    
    if (VERBOSE)
      cerr << "reading regions file" << endl;
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    assert(check_sorted(regions));
    
    if (VERBOSE)
      cerr << "loading read sequences" << endl;
    vector<string> names, sequences;
    read_fasta_file(reads_file.c_str(), names, sequences);
    if (VERBOSE)
      cerr << "read " << names.size() << " sequences" << endl;
    
    if (VERBOSE)
      cerr << "associating read sequences to mapped locations" << endl;
    vector<size_t> lookup;
    relative_sort(mapped_locations, names, lookup);
    vector<string> seq_swapper(sequences.size());
    for (size_t i = 0; i < lookup.size(); ++i)
      seq_swapper[i].swap(sequences[lookup[i]]);
    sequences.swap(seq_swapper);
    
    if (VERBOSE)
      cerr << "associating reads with regions" << endl;
    vector<vector<GenomicRegion> > clusters;
    vector<vector<string> > cluster_seqs;
    separate_regions(regions, mapped_locations, sequences, 
		     clusters, cluster_seqs);
    mapped_locations.clear();
    sequences.clear();
    
    if (VERBOSE)
      cerr << "adjusting regions to cover reads" << endl;
    adjust_region_ends(clusters, regions);
    
    if (VERBOSE)
      cerr << "extracting reference sequence" << endl;
    vector<string> region_sequences;
    extract_regions_fasta(chrom_dir, regions, region_sequences);

    size_t min_obs_for_confidence = 
      get_min_obs_for_confidence(interval_width, alpha);
    
    std::ostream *out = (outfile.empty()) ? &cout : 
      new std::ofstream(outfile.c_str());
    for (size_t i = 0; i < clusters.size(); ++i)
      *out << call_cpg_state(crit, alpha, min_obs_for_confidence, interval_width,
			     clusters[i], cluster_seqs[i], regions[i], 
			     region_sequences[i]) << endl;
    if (out != &cout) delete out;
  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
