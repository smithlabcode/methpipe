/*    checkoverlap: a program for identifying overlapping ends of
 *                  mapped paired end reads.
 *
 *    Copyright (C) 2010 University of Southern California and
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
#include <unistd.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "MappedRead.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::numeric_limits;
using std::ofstream;

inline static size_t 
get_distance(const MappedRead &a, const MappedRead &b) {
  return (a.r.pos_strand()) ? b.r.get_end() - a.r.get_start() :
    a.r.get_end() - b.r.get_start();
}


static void
fix_overlap(const MappedRead &a, MappedRead &b) {
  if (a.r.get_strand() == b.r.get_strand()) {
    if (a.r.get_start() < b.r.get_start()) {
      const size_t overlap = a.r.get_end() - b.r.get_start();
      assert(overlap > 0);
      b.seq = (a.r.pos_strand()) ?
	string(overlap, 'N') + b.seq.substr(overlap) :
	b.seq.substr(0, b.seq.length() - overlap) + string(overlap, 'N');
    }
    else {
      const size_t overlap = b.r.get_end() - a.r.get_start();
      assert(overlap > 0);
      b.seq = (a.r.pos_strand()) ? 
	b.seq.substr(0, b.seq.length() - overlap) + string(overlap, 'N') :
	string(overlap, 'N') + b.seq.substr(overlap);
    }
  }
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
mask_less_informative(MappedRead &one, MappedRead &two) {
  const size_t informative_one = one.seq.length() - 
    count(one.seq.begin(), one.seq.end(), 'N') - one.r.get_score();
  const size_t informative_two = two.seq.length() - 
    count(two.seq.begin(), two.seq.end(), 'N') - two.r.get_score();
  if (informative_one > informative_two) fix_overlap(one, two);
  else fix_overlap(two, one);
}


int 
main(int argc, const char **argv) {

  try {

      cerr << "############################################################" << endl
           << "############################################################" << endl
           << "############################################################" << endl
           << "THIS PROGRAM checkoverlap IS TO BE DEPRECIATED" << endl
           << "Instead use clipmates" << endl 
           << "############################################################" << endl
           << "############################################################" << endl
           << "############################################################" << endl;

    bool VERBOSE = false;
    size_t max_distance = 500;
    string histogram_file;
  
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program for identifying overlapping ends of "
			   "mapped paired end reads.",
			   "<end-1-in> <end-2-in> <end-1-out> <end-2-out>");
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("hist", 'h', "print frag size histogram in this file", 
		      false, histogram_file);
    opt_parse.add_opt("max-dist", 'm', "max distance to print", 
		      false, max_distance);
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
    if (leftover_args.size() != 4) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string end_one_file = leftover_args.front(); 
    const string end_two_file = leftover_args[1];
    const string end_one_out = leftover_args[2];
    const string end_two_out = leftover_args.back(); 
    /****************** END COMMAND LINE OPTIONS *****************/

    const bool HISTOGRAM = !histogram_file.empty();
    
    ifstream in_one(end_one_file.c_str());
    ifstream in_two(end_two_file.c_str());
    
    ofstream out_one(end_one_out.c_str());
    ofstream out_two(end_two_out.c_str());
    
    MappedRead one, two;
    
    unordered_map<size_t, size_t> distance_hist;

    in_two >> two;
    while (!in_one.eof() && !in_two.eof()) {
      in_one >> one;
      while (!in_two.eof() && name_smaller(two, one)) {
	out_two << two << '\n';
	in_two >> two;
      }
      if (same_read(one, two)) {
	if (one.r.overlaps(two.r))
	  mask_less_informative(one, two);
	if (HISTOGRAM) {
	  const size_t distance = get_distance(one, two);
	  if (distance_hist.find(distance) == distance_hist.end())
	    distance_hist[distance] = 0;
	  distance_hist[distance]++;
	}
      }
      out_one << one << '\n';
    }
    out_two << two << '\n';
    while (!in_two.eof()) {
      in_two >> two;
      out_two << two << '\n';
    }
    while (!in_one.eof()) {
      in_one >> one;
      out_one << one << '\n';
    }
    
    if (HISTOGRAM) {
      ofstream out(histogram_file.c_str());
      for (size_t i = 0; i <= max_distance; ++i)
	out << i << "\t" 
	    << ((distance_hist.find(i) != distance_hist.end()) ?
		distance_hist[i] : 0) << endl;
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
