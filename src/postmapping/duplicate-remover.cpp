/*    duplicate-remover: a program to select unique reads mapped to
 *    the same location and the same strand
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith, Song Qiang
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
#include <numeric>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "RNG.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;

/**************** FOR CLARITY BELOW WHEN COMPARING MAPPED READS **************/
static inline bool
same_strand(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() == b.get_strand();
}
static inline bool
strand_leq(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() <= a.get_strand();
}
static inline bool
chrom_less(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_chrom() < b.get_chrom();
}
static inline bool
same_start(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() == b.get_start();
}
static inline bool
start_less(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() < b.get_start();
}
static inline bool
same_end(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_end() == b.get_end();
}
static inline bool
end_less(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_end() < b.get_end();
}
/******************************************************************************/


static inline bool
check_sorted(const MappedRead &prev_mr, const MappedRead &mr, 
	     const bool CHECK_SECOND_ENDS) {
  if (!CHECK_SECOND_ENDS)
    if (chrom_less(prev_mr.r, mr.r) || 
	(prev_mr.r.same_chrom(mr.r) &&
	 (start_less(prev_mr.r, mr.r) ||
	  (same_start(prev_mr.r, mr.r) &&
	   (end_less(prev_mr.r, mr.r) ||
	    (same_end(prev_mr.r, mr.r) && 
	     strand_leq(prev_mr.r, mr.r)))))))
      return true;
    else {
      std::ostringstream oss;
      oss << "READS NOT SORTED: " << endl
	  << prev_mr << endl
	  << mr << endl;
      throw RMAPException(oss.str());
    }
  else {
    if (chrom_less(prev_mr.r, mr.r) ||
	(prev_mr.r.same_chrom(mr.r) &&
	 (end_less(prev_mr.r, mr.r) || 
	  (same_end(prev_mr.r, mr.r) &&
	   (start_less(prev_mr.r, mr.r) || 
	    (same_start(prev_mr.r, mr.r) && 
	     strand_leq(prev_mr.r, mr.r)))))))
      return true;
    else {
      std::ostringstream oss;
      oss << "READS NOT SORTED (OPTION -B): " << endl
	  << prev_mr << endl << mr << endl;
      throw RMAPException(oss.str());
    }
  }
  return false;
}


static void
collect_duplicate_stats(const MappedRead &mr, size_t &total_reads, 
			size_t &total_bases, size_t &total_valid_bases,
			size_t &total_mismatches) {
  ++total_reads;
  total_bases += mr.seq.length();
  total_valid_bases += mr.seq.length() - count(mr.seq.begin(), mr.seq.end(), 'N');
  total_mismatches += static_cast<size_t>(mr.r.get_score());
}


static bool
is_complete_fragment(const MappedRead &mr) {
  static const char label[] = "FRAG:";
  static const size_t label_len = 5;
  const string name(mr.r.get_name());
  return lexicographical_compare(name.begin(), name.begin() + label_len,
				 label, label + label_len);
}


struct DuplicateFragmentTester {
  bool CHECK_SECOND_ENDS;
  DuplicateFragmentTester(const bool x) :
    CHECK_SECOND_ENDS(x) {}
  bool operator()(const MappedRead &lhs, const MappedRead &rhs) const {
    
    if (!same_strand(lhs.r, rhs.r) || !lhs.r.same_chrom(rhs.r)) 
      return false;
    
    const bool lhs_is_frag = is_complete_fragment(lhs);
    const bool rhs_is_frag = is_complete_fragment(rhs);
    
    if ((lhs_is_frag && rhs_is_frag) || (!lhs_is_frag && !rhs_is_frag))
      return same_start(lhs.r, rhs.r) && same_end(lhs.r, rhs.r);
    
    if (CHECK_SECOND_ENDS)
      return (lhs.r.pos_strand() && same_end(lhs.r, rhs.r)) ||
	(lhs.r.neg_strand() && same_start(lhs.r, rhs.r));
    
    return (lhs.r.pos_strand() && same_start(lhs.r, rhs.r)) ||
      (lhs.r.neg_strand() && same_end(lhs.r, rhs.r));
  }
};


static size_t
get_representative_read(vector<MappedRead> &mapped_ties) {
  static const Runif rng(time(NULL) + getpid());
  
  vector<MappedRead>::const_iterator iter =
    std::partition(mapped_ties.begin(), mapped_ties.end(), is_complete_fragment);
  
  return rng.runif(0ul, (iter != mapped_ties.begin()) ?
		   iter - mapped_ties.begin() : mapped_ties.size());
}


static void
remove_duplicates(const bool CHECK_SECOND_ENDS,
		  const string &infile, const string &outfile, 
		  size_t &total_reads, size_t &unique_reads,
		  size_t &total_bases, size_t &total_valid_bases, 
		  size_t &total_mismatches) {
  
  const DuplicateFragmentTester dft(CHECK_SECOND_ENDS);
  
  MappedRead mr, representative;
  vector<MappedRead> mapped_ties;
  
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
  
  std::ifstream in(infile.c_str());
  
  if (!(in >> mr)) 
    throw RMAPException("mapped read file seems empty: " + infile);
  
  mapped_ties.push_back(mr);
  
  while (in >> mr) {
    ++total_reads;
    check_sorted(mapped_ties.back(), mr, CHECK_SECOND_ENDS);
    if (!dft(mapped_ties.front(), mr)) {
      const size_t rep_idx = get_representative_read(mapped_ties);
      out << mapped_ties[rep_idx] << endl;
      collect_duplicate_stats(mapped_ties[rep_idx], unique_reads, total_bases,
			      total_valid_bases, total_mismatches);
      mapped_ties.clear();
    }
    mapped_ties.push_back(mr);
  }
  const size_t rep_idx = get_representative_read(mapped_ties);
  out << mapped_ties[rep_idx] << endl;
  collect_duplicate_stats(mapped_ties[rep_idx], unique_reads,
			  total_bases, total_valid_bases,
			  total_mismatches);
}



int main(int argc, const char **argv) {

  try {
    bool VERBOSE = false;
    string outfile;
    string statfile;
    bool CHECK_SECOND_ENDS = false;//remove duplicates by end

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "program to remove duplicate reads from "
			   "a sorted file of mapped reads",
			   "<mapped-reads-file>");
    opt_parse.add_opt("output", 'o', "output file for unique reads",
		      false, outfile);
    opt_parse.add_opt("endtwo", 'B', "[when PE reads are involved] "
		      "similar fragment check on 5' end of second mate",
		      false, CHECK_SECOND_ENDS);
    opt_parse.add_opt("stats", 'S', "statistics output file", false, statfile);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string infile = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    size_t total_reads = 0, unique_reads = 0, total_bases = 0, 
      total_valid_bases = 0, total_mismatches = 0;
    
    remove_duplicates(CHECK_SECOND_ENDS, infile, outfile, total_reads, 
		      unique_reads, total_bases, total_valid_bases, 
		      total_mismatches); 
    
    if (!statfile.empty() || VERBOSE) {
      std::ofstream of;
      if (!statfile.empty()) of.open(statfile.c_str());
      std::ostream out(statfile.empty() ? std::clog.rdbuf() : of.rdbuf());
      out << "TOTAL READS:\t" << total_reads << endl
	  << "UNIQUE READS:\t" << unique_reads << endl
	  << "TOTAL BASES:\t" << total_bases << endl
	  << "VALID BASES:\t" << total_valid_bases << endl
	  << "MISMATCHES:\t" << total_mismatches << endl
	  << "MISMATCH FREQ:\t" << total_mismatches/
	std::max(1.0, static_cast<double>(total_bases)) << endl;
    }
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
