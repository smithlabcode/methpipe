/*    duplicate-remover:
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
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "RNG.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

/**************** FOR CLARITY BELOW WHEN COMPARING MAPPED READS **************/
static inline bool
same_strand(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() == b.get_strand();
}
static inline bool
strand_leq(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() <= b.get_strand();
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


class OrderChecker {
public:
  bool operator()(const MappedRead &prev_mr, const MappedRead &mr) const {
    return (CHECK_SECOND_ENDS ? end_two_check(prev_mr.r, mr.r) :
	    end_one_check(prev_mr.r, mr.r));
  }
  static void 
  set_comparison_mode(const bool m) {CHECK_SECOND_ENDS = m;}
  
private:
  static bool CHECK_SECOND_ENDS;

  static bool
  end_one_check(const GenomicRegion &prev_mr, const GenomicRegion &mr) {
    return (chrom_less(prev_mr, mr) || 
	    (prev_mr.same_chrom(mr) &&
	     (start_less(prev_mr, mr) ||
	      (same_start(prev_mr, mr) &&
	       (end_less(prev_mr, mr) ||
		(same_end(prev_mr, mr) && strand_leq(prev_mr, mr)))))));
  }
  
  static bool 
  end_two_check(const GenomicRegion &prev_mr, const GenomicRegion &mr) {
    return (chrom_less(prev_mr, mr) ||
	    (prev_mr.same_chrom(mr) &&
	     (end_less(prev_mr, mr) || 
	      (same_end(prev_mr, mr) &&
	       (start_less(prev_mr, mr) || 
		(same_start(prev_mr, mr) && strand_leq(prev_mr, mr)))))));
  }
};
bool OrderChecker::CHECK_SECOND_ENDS = false;

class DuplicateFragmentTester {
public:

  bool operator()(const MappedRead &lhs, const MappedRead &rhs) const {
    return the_real_test(lhs.r, rhs.r);
  }

  static void 
  set_comparison_mode(const bool m) {CHECK_SECOND_ENDS = m;}

  static bool
  is_complete_fragment(const MappedRead &mr) {
    return internal_is_complete_fragment(mr.r);
  }
private:

  static bool CHECK_SECOND_ENDS;

  static bool
  the_real_test(const GenomicRegion &lhs, const GenomicRegion &rhs) {
    if (!same_strand(lhs, rhs) || !lhs.same_chrom(rhs)) 
      return false;
    
    const bool lhs_is_frag = internal_is_complete_fragment(lhs);
    const bool rhs_is_frag = internal_is_complete_fragment(rhs);
    
    // if ((lhs_is_frag && rhs_is_frag) || (!lhs_is_frag && !rhs_is_frag))
    if (lhs_is_frag == rhs_is_frag)
      return same_start(lhs, rhs) && same_end(lhs, rhs);
    
    if (CHECK_SECOND_ENDS)
      return (lhs.pos_strand() && same_end(lhs, rhs)) ||
	(lhs.neg_strand() && same_start(lhs, rhs));
    
    return (lhs.pos_strand() && same_start(lhs, rhs)) ||
      (lhs.neg_strand() && same_end(lhs, rhs));
  }
  static bool
  internal_is_complete_fragment(const GenomicRegion &mr) {
    static const char label[] = "FRAG:";
    static const size_t label_len = 5;
    const string name(mr.get_name());
    return lexicographical_compare(name.begin(), name.begin() + label_len,
				   label, label + label_len);
  }

};
bool DuplicateFragmentTester::CHECK_SECOND_ENDS = false;


class LessMismatchCmp: public std::binary_function<MappedRead, MappedRead, bool>
{
public:
    bool operator()(const MappedRead &a, const MappedRead &b) const
    {return a.r.get_score() < b.r.get_score();}
};

class SameMismatchCmp: public std::binary_function<MappedRead, MappedRead, bool>
{
public:
    bool operator()(const MappedRead &a, const MappedRead &b) const
    {return a.r.get_score() == b.r.get_score();}
};

static size_t
get_representative_read(vector<MappedRead> &candidates) {
  static const Runif rng(time(NULL) + getpid());
  vector<MappedRead>::iterator iter =
    std::partition(candidates.begin(), candidates.end(), 
		   &DuplicateFragmentTester::is_complete_fragment);
  iter = iter != candidates.begin() ? iter : candidates.end();
  const vector<MappedRead>::const_iterator min_mismatch_iter =
      std::min_element(candidates.begin(), iter, LessMismatchCmp());
  iter = std::partition(candidates.begin(), iter,
                        bind2nd(SameMismatchCmp(), *min_mismatch_iter));
  return rng.runif(static_cast<size_t>(0), iter - candidates.begin());
}


static inline size_t
count_good_bases(const MappedRead &mr) {
  return mr.r.get_width() - count(mr.seq.begin(), mr.seq.end(), 'N');
}


static void
remove_duplicates(const string &infile, const string &outfile, 
		  size_t &total_reads, size_t &unique_reads,
		  size_t &total_bases_in, size_t &good_bases_out) {
  
  DuplicateFragmentTester dft;
  OrderChecker oc;
  
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
  
  std::ifstream ifs;
  if (!infile.empty()) ifs.open(infile.c_str());
  std::istream in(infile.empty() ? cin.rdbuf() : ifs.rdbuf());
  
  MappedRead mr;
  if (!(in >> mr)) 
    throw SMITHLABException("mapped read file seems empty: " + infile);
  
  ++total_reads;
  
  vector<MappedRead> candidates;
  candidates.push_back(mr);
  
  while (in >> mr) {
    ++total_reads;
    total_bases_in += mr.r.get_width();
    
    if (!oc(candidates.back(), mr))
      throw SMITHLABException("NOT SORTED:\n" + candidates.back().r.tostring() + 
			  "\n" + mr.r.tostring());
    
    if (!dft(candidates.front(), mr)) {
      const size_t rep_idx = get_representative_read(candidates);
      out << candidates[rep_idx] << '\n';
      ++unique_reads;
      good_bases_out += count_good_bases(candidates[rep_idx]);
      candidates.clear();
    }
    
    candidates.push_back(mr);
  }
  const size_t rep_idx = get_representative_read(candidates);
  out << candidates[rep_idx] << '\n';
  ++unique_reads;
  good_bases_out += count_good_bases(candidates[rep_idx]);
}


int main(int argc, const char **argv) {

  try {
    bool VERBOSE = false;
    string outfile;
    string statfile;
    
    bool CHECK_SECOND_ENDS = false;
    bool INPUT_FROM_STDIN = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "program to remove duplicate reads from "
			   "sorted mapped reads", "<mapped-reads>");
    opt_parse.add_opt("output", 'o', "output file for unique reads",
		      false, outfile);
    opt_parse.add_opt("stdin", '\0', "take input from stdin",
		      false, INPUT_FROM_STDIN);
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
    if (!INPUT_FROM_STDIN && leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    string infile;
    if (!leftover_args.empty())
      infile = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    size_t total_reads = 0, unique_reads = 0, total_bases_in = 0, good_bases_out = 0;

    OrderChecker::set_comparison_mode(CHECK_SECOND_ENDS);
    DuplicateFragmentTester::set_comparison_mode(CHECK_SECOND_ENDS);
    
    remove_duplicates(infile, outfile, total_reads, 
		      unique_reads, total_bases_in, good_bases_out);
    
    if (!statfile.empty() || VERBOSE) {
      std::ofstream of;
      if (!statfile.empty()) of.open(statfile.c_str());
      std::ostream out(statfile.empty() ? std::clog.rdbuf() : of.rdbuf());
      out << "TOTAL READS:\t" << total_reads << endl
	  << "UNIQUE READS:\t" << unique_reads << endl
	  << "ALL BASES IN:\t" << total_bases_in << endl
	  << "GOOD BASES OUT:\t" << good_bases_out << endl;
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
