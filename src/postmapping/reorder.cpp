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
using std::priority_queue;

struct MappedReadOrder {
  static bool CHECK_SECOND_ENDS;
  static size_t max_frag_len;
  bool operator()(const MappedRead &lhs, const MappedRead &rhs) const {
    return rhs.r <= lhs.r;
  }
  static bool 
  is_ready(const priority_queue<MappedRead, vector<MappedRead>, MappedReadOrder> &pq,
	   const MappedRead &mr) {
    return !pq.top().r.same_chrom(mr.r) || pq.top().r.get_end() < mr.r.get_start();
  }
};

bool MappedReadOrder::CHECK_SECOND_ENDS = false;
size_t MappedReadOrder::max_frag_len = 500;

int main(int argc, const char **argv) {
  
  try {
    
    string outfile;
    size_t max_frag_len = 500;
    bool CHECK_SECOND_ENDS = true; //remove duplicates by end
    bool INPUT_FROM_STDIN = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program to re-sort genomic intervals "
			   "based on start or end coordinates",
			   "<mapped-reads-file>");
    opt_parse.add_opt("output", 'o', "output file for unique reads",
		      false, outfile);
    opt_parse.add_opt("stdin", '\0', "take input from stdin",
		      false, INPUT_FROM_STDIN);
    //     opt_parse.add_opt("endtwo", 'B', "[when PE reads are involved] "
    // 		      "similar fragment check on 5' end of second mate",
    // 		      false, CHECK_SECOND_ENDS);
    opt_parse.add_opt("len", 'L', "max fragment length", false, max_frag_len);
    
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
    
    MappedReadOrder::CHECK_SECOND_ENDS = CHECK_SECOND_ENDS;
    MappedReadOrder::max_frag_len = max_frag_len;
    priority_queue<MappedRead, vector<MappedRead>, MappedReadOrder> pq;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    std::ifstream ifs;
    if (!INPUT_FROM_STDIN) ifs.open(infile.c_str());
    std::istream in(INPUT_FROM_STDIN ? std::cin.rdbuf() : ifs.rdbuf());
    
    MappedRead mr; //, prev_in;
    
    while (in >> mr) {
      // assert(pq.empty() || !mr.r.less1(prev_in.r));
      while (!pq.empty() && MappedReadOrder::is_ready(pq, mr)) {
	// assert(prev_out.r.get_end() == 0 || prev_out.r <= pq.top().r);
	// prev_out = pq.top();
	out << pq.top() << '\n';
	pq.pop();
      }
      pq.push(mr);
      // prev_in = mr;
    }
    
    while (!pq.empty()) {
      // assert(prev_out.r.get_end() == 0 ||prev_out.r <= pq.top().r);
      // prev_out = pq.top();
      out << pq.top() << '\n';
      pq.pop();
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
