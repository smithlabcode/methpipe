/*    clipmates: a program for merging mates in paired-end BS-seq data
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith, Song Qiang and Elena Harris
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

#include <unistd.h>

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ofstream;
using std::ifstream;
using std::copy;

static bool
check_sorted_by_ID(MappedRead &prev_mr, const MappedRead &mr) {
  if (prev_mr.r.get_name() > mr.r.get_name()) {
    // return false;
    std::ostringstream oss;
    oss << "READS NOT SORTED: " << endl
	<< prev_mr << endl
	<< mr << endl;
    throw SMITHLABException(oss.str());
  }
  prev_mr = mr;
  return true;
}


static void
revcomp(MappedRead &mr) {
  // set the strand to the opposite of the current value
  mr.r.set_strand(mr.r.pos_strand() ? '-' : '+');
  // reverse complement the sequence, and reverse the quality scores
  revcomp_inplace(mr.seq);
  std::reverse(mr.scr.begin(), mr.scr.end());
}


inline static bool
same_read(const MappedRead &a, const MappedRead &b) {
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  return std::equal(sa.begin(), sa.end() - 1, sb.begin());
}


inline static bool
name_smaller(const MappedRead &a, const MappedRead &b) {
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  return std::lexicographical_compare(sa.begin(), sa.end() - 1, 
				      sb.begin(), sb.end() - 1);
}


static void
merge_mates(const size_t MAX_SEGMENT_LENGTH,
	    const MappedRead &one, const MappedRead &two,
            MappedRead &merged, int &len) {
  
  merged = one;
  size_t start_one = 0, end_one = 0, start_two = 0, 
    end_two = 0, start_overlap = 0, end_overlap = 0;
  if (merged.r.pos_strand()) {
    start_overlap = std::max(one.r.get_start(), two.r.get_start());
    end_overlap = std::min(one.r.get_end(), two.r.get_end());
    start_one = one.r.get_start();
    end_one = std::min(start_overlap, one.r.get_end());
    start_two = std::max(end_overlap, two.r.get_start());
    end_two = two.r.get_end();
    len = end_two - start_one;
    merged.r.set_start(start_one);
    merged.r.set_end(end_two);
  }
  else { // if (merged.r.neg_strand())
    start_overlap = std::max(one.r.get_start(), two.r.get_start());
    end_overlap = std::min(one.r.get_end(), two.r.get_end());
    start_one = std::max(end_overlap, one.r.get_start());
    end_one = one.r.get_end();
    start_two = two.r.get_start();
    end_two = std::min(start_overlap, two.r.get_end());
    len = end_one - start_two;
    merged.r.set_start(start_two);
    merged.r.set_end(end_one);
  }
  
  assert(end_one >= start_one && end_two >= start_two);
  assert(start_overlap >= end_overlap || 
	 static_cast<size_t>(len) == end_one - start_one
	 + end_two - start_two + end_overlap - start_overlap);
  
  merged.r.set_score(one.r.get_score() + two.r.get_score());
  
  if (len > 0 && len <= static_cast<int>(MAX_SEGMENT_LENGTH)) {
    merged.seq = string(len, 'N');
    merged.scr = string(len, 'B');
    const string name(one.r.get_name());
    merged.r.set_name("FRAG:" + name.substr(0, name.size() - 1));

    // lim_one ios the offset within the merged sequence where the
    // overlapping portion begins
    const size_t lim_one = end_one - start_one;
    copy(one.seq.begin(), one.seq.begin() + lim_one, merged.seq.begin());
    copy(one.scr.begin(), one.scr.begin() + lim_one, merged.scr.begin());
    const size_t lim_two = end_two - start_two;
    copy(two.seq.end() - lim_two, two.seq.end(), merged.seq.end() - lim_two);
    copy(two.scr.end() - lim_two, two.scr.end(), merged.scr.end() - lim_two);
    
    // deal with overlapping part
    if (start_overlap < end_overlap) {

      const int info_one = one.seq.length() - 
	count(one.seq.begin(), one.seq.end(), 'N') - 
	static_cast<int>(one.r.get_score());
      const int info_two = two.seq.length() - 
	count(two.seq.begin(), two.seq.end(), 'N') - 
	static_cast<int>(two.r.get_score());

      // use the mate with the most information ("info") to fill in
      // the overlapping portion
      if (info_one >= info_two) {
	const size_t source_start = merged.r.pos_strand() ?
	  (start_overlap - one.r.get_start()) : (one.r.get_end() - end_overlap);
	const size_t source_end = merged.r.pos_strand() ?
	  (end_overlap -  one.r.get_start()) : (one.r.get_end() - start_overlap);
	copy(one.seq.begin() + source_start, one.seq.begin() + source_end,
	     merged.seq.begin() + lim_one);
	copy(one.scr.begin() + source_start, one.scr.begin() + source_end,
	     merged.scr.begin() + lim_one);
      }
      else { // if (info_two > info_one)
	const size_t source_start = merged.r.pos_strand() ?
	  (start_overlap - two.r.get_start()) : (two.r.get_end() - end_overlap);
	const size_t source_end = merged.r.pos_strand() ?
	  (end_overlap -  two.r.get_start()) : (two.r.get_end() - start_overlap);
	copy(two.seq.begin() + source_start, two.seq.begin() + source_end,
	     merged.seq.begin() + lim_one);
	copy(two.scr.begin() + source_start, two.scr.begin() + source_end,
	     merged.scr.begin() + lim_one);
      }
    }
  }
}


int 
main(int argc, const char **argv)  {
  try {

    size_t MAX_SEGMENT_LENGTH = 500;
        
    bool VERBOSE = false;
    bool REVCOMP = true;
    string histogram_file;
    string end_one_file, end_two_file; // mapped read input
    string end_one_out, end_two_out;  // mapped read output
    string out_stat;
    string outfile;
      
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program to merge paired-end BS-seq reads");
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size", 
		      false, MAX_SEGMENT_LENGTH); 
    opt_parse.add_opt("out_stat", 'S', "statistics output file", false, out_stat); 
    opt_parse.add_opt("t-rich", 'T', "input file for T-rich mates (s_*_1_...)",
		      true, end_one_file); 
    opt_parse.add_opt("a-rich", 'A', "input file for A-rich mates (s_*_2_...)",
		      true, end_two_file); 
    opt_parse.add_opt("outfile", 'o', "Output file name", false, outfile);
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
    /****************** END COMMAND LINE OPTIONS *****************/
    ifstream in_one(end_one_file.c_str());
    ifstream in_two(end_two_file.c_str());
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    /*************Used for collecting statistics***************/
    vector<size_t> frag_len_distr(MAX_SEGMENT_LENGTH + 1, 0);
    size_t merged_pairs = 0;
    size_t incorrect_chr = 0;
    size_t incorrect_strand = 0;
    size_t problem_x = 0;
    size_t incorrect_orient = 0;
    size_t incorrect_frag_size = 0;
    size_t broken_pairs = 0;
    /*************End for collecting statistics**************/

    MappedRead one, two, prev_one, prev_two;
    prev_one.r.set_name("");
    prev_two.r.set_name("");
    
    bool one_is_good = ((in_one >> one) && check_sorted_by_ID(prev_one, one));
    bool two_is_good = ((in_two >> two) && check_sorted_by_ID(prev_two, two));
    if (REVCOMP) revcomp(two);
    
    while (one_is_good && two_is_good) {
      if (same_read(one, two)) { // one and tow are mates
	
	if (!one.r.same_chrom(two.r)) {
	  incorrect_chr++;
	  out << one << endl << two << endl;
	}
	else if (one.r.get_strand() != two.r.get_strand()) {
	  if (one.r.get_start() == two.r.get_start()) {
	    problem_x++;
	    out << one << endl;
	  }
	  else {
	    incorrect_strand++;
	    out << one << endl << two << endl;
	  }
	}
	else {
	  MappedRead merged;
	  int len = 0;
	  merge_mates(MAX_SEGMENT_LENGTH, one, two, merged, len);
	  if (len > 0 && len <= static_cast<int>(MAX_SEGMENT_LENGTH)) {
	    merged_pairs++;
	    frag_len_distr[len]++;
	    out << merged << endl;
	  }
	  else {
	    out << one << endl << two << endl;
	    if (len < 0) incorrect_orient++;
	    else incorrect_frag_size++;
	  }
	}
	one_is_good = ((in_one >> one) && (check_sorted_by_ID(prev_one, one)));
	two_is_good = ((in_two >> two) && (check_sorted_by_ID(prev_two, two)));
	if (REVCOMP) revcomp(two);
      } 
      else if (name_smaller(one, two)) {
	out << one << endl;
	broken_pairs++;
	one_is_good = ((in_one >> one) && (check_sorted_by_ID(prev_one, one)));
      }
      else { // one comes after two
	out << two << endl;
	broken_pairs++;
	two_is_good = ((in_two >> two) && (check_sorted_by_ID(prev_two, two)));
	if (REVCOMP) revcomp(two);
      }
    }
    while (one_is_good) {
      out << one << endl;
      broken_pairs++;            
      one_is_good = ((in_one >> one) && (check_sorted_by_ID(prev_one, one)));
    }
    while (two_is_good) {
      out << two << endl;
      broken_pairs++;
      two_is_good = ((in_two >> two) && (check_sorted_by_ID(prev_two, two)));
      if (REVCOMP) revcomp(two);
    }
      
    if (!out_stat.empty()) {
      ofstream outst(out_stat.c_str());
      outst << "TOTAL MERGED PAIRS (# OF PAIRS):\t" << merged_pairs << endl
	    << "INCORRECTLY MAPPED TO DIFFERENT CHROM:\t" << incorrect_chr << endl
	    << "UNMERGED BY PROBLEM X:\t" << problem_x << endl
	    << "UNMERGED BY STRAND INCOMPATIBILITY:\t" << incorrect_strand << endl
	    << "UNMERGED BY INCONSISTENT ORIENTATION:\t" << incorrect_orient << endl
	    << "UNMERGED BY INSERT SIZE:\t" << incorrect_frag_size << endl
	    << "TOTAL MAPPED BROKEN PAIRS (MISSING MATES):\t" << broken_pairs << endl
	    << "INSERT_LENGTH" << "\t" << "FREQUENCY" << endl;
      for (size_t i = 0; i < MAX_SEGMENT_LENGTH + 1; i++)
	outst << i << "\t" << frag_len_distr[i] << endl;
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
