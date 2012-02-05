/*    reorder: a program to re-sort genomic intervals based on either
 *    the start or end coordinate, assuming the intervals are already
 *    sorted on the other.
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
#include <numeric>
#include <queue>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::priority_queue;


/**************** FOR CLARITY BELOW WHEN COMPARING MAPPED READS **************/
static inline bool
same_strand(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() == b.get_strand();
}
static inline bool
strand_less(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_strand() <= b.get_strand();
}
static inline bool
start_leq(const GenomicRegion &a, const GenomicRegion &b) {
  return a.get_start() <= b.get_start();
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


struct OrderChecker {
  bool operator()(const MappedRead &prev_mr, const MappedRead &mr) const {
    return end_two_check(mr.r, prev_mr.r);
  }
  static bool 
  is_ready(const priority_queue<MappedRead, vector<MappedRead>, OrderChecker> &pq,
	   const MappedRead &mr) {
    return !pq.top().r.same_chrom(mr.r) || pq.top().r.get_end() < mr.r.get_start();
  }
  static bool 
  end_two_check(const GenomicRegion &prev_mr, const GenomicRegion &mr) {
    return (end_less(prev_mr, mr) || 
	    (same_end(prev_mr, mr) &&
	     (strand_less(prev_mr, mr) ||
	      (same_strand(prev_mr, mr) && start_leq(prev_mr, mr)))));
  }
};

int main(int argc, const char **argv) {
  
  try {
    
    string outfile;
    bool INPUT_FROM_STDIN = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "a program to re-sort genomic "
			   "intervals based on start or end coordinates",
			   "<mapped-reads-file>");
    opt_parse.add_opt("output", 'o', "output file for unique reads",
		      false, outfile);
    opt_parse.add_opt("stdin", '\0', "take input from stdin",
		      false, INPUT_FROM_STDIN);
    
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
    const string infile = INPUT_FROM_STDIN ? "" : leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    std::priority_queue<MappedRead, vector<MappedRead>, OrderChecker> pq;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    std::ifstream ifs;
    if (!INPUT_FROM_STDIN) ifs.open(infile.c_str());
    std::istream in(INPUT_FROM_STDIN ? std::cin.rdbuf() : ifs.rdbuf());
    
    MappedRead mr;
    while (in >> mr) {
      while (!pq.empty() && OrderChecker::is_ready(pq, mr)) {
	out << pq.top() << '\n';
	pq.pop();
      }
      pq.push(mr);
    }
    while (!pq.empty()) {
      out << pq.top() << '\n';
      pq.pop();
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
