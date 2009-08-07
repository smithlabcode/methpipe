/*    coverage2: a program for determining genomic coverage in a
 *    resequencing project
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
#include <popt.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include "BSUtils.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::pair;
using std::make_pair;


static void
get_coverage(const vector<SimpleGenomicRegion> &mapped, size_t start, size_t end,
	     std::ostream &out) {
  
  const string chrom(mapped[start].get_chrom());
  
  vector<pair<size_t, bool> > ends(2*(end - start));
  for (size_t i = start; i != end; ++i) {
    const size_t id = 2*(i - start);
    ends[id] = make_pair(mapped[i].get_start(), false);
    ends[id + 1] = make_pair(mapped[i].get_end(), true);
  }
  sort(ends.begin(), ends.end());
  
  size_t prev_pos = ends.front().first, prev_depth = 0;
  for (size_t i = 0; i < ends.size(); ) {
    const size_t pos = ends[i].first;
    size_t depth = prev_depth;
    while (i < ends.size() && ends[i].first == pos) {
      depth += (ends[i].second) ? -1 : 1;
      ++i;
    }
    if (depth != prev_depth) {
      if (prev_depth > 0)
	out << chrom << '\t' << prev_pos 
	    << '\t' << pos << '\t' << prev_depth << '\n';
      prev_pos = pos;
      prev_depth = depth;
    }
  }
  if (prev_depth > 0)
    out << chrom << '\t' << prev_pos 
	<< '\t' << ends.back().first << '\t' << prev_depth << '\n';
}

static void
find_chrom_changes(const vector<SimpleGenomicRegion>::const_iterator &a, 
		   const vector<SimpleGenomicRegion>::const_iterator &b, 
		   vector<size_t> &changes) {
  changes.push_back(0);
  string prev_chrom(a->get_chrom());
  for (vector<SimpleGenomicRegion>::const_iterator i(a); i < b; ++i) {
    if (i->get_chrom() != prev_chrom) {
      changes.push_back(i - a);
      prev_chrom = i->get_chrom();
    }
  }
}


static void
get_coverage_wig(const vector<SimpleGenomicRegion> &mapped, size_t start, size_t end,
		 std::ostream &out) {
  
  const string chrom(mapped[start].get_chrom());
  
  vector<pair<size_t, bool> > ends(2*(end - start));
  for (size_t i = start; i != end; ++i) {
    const size_t id = 2*(i - start);
    ends[id] = make_pair(mapped[i].get_start(), false);
    ends[id + 1] = make_pair(mapped[i].get_end(), true);
  }
  sort(ends.begin(), ends.end());
  
  bool print_header = true;
  size_t prev_pos = ends.front().first, prev_depth = 0;
  for (size_t i = 0; i < ends.size(); ) {
    const size_t pos = ends[i].first;
    size_t depth = prev_depth;
    while (i < ends.size() && ends[i].first == pos) {
      depth += (ends[i].second) ? -1 : 1;
      ++i;
    }
    if (depth != prev_depth) {
      if (prev_depth > 0) {
	if (print_header) {
	  out << "fixedStep chrom=" << chrom 
	      << "\tstart=" << prev_pos << "\tstep=1" << '\n';
	  print_header = false;
	}
	for (size_t j = prev_pos; j < pos; ++j)
	  out << prev_depth << '\n';
      }
      else print_header = true; 
      prev_pos = pos;
      prev_depth = depth;
    }
  }
  if (prev_depth > 0) {
    if (print_header)
      out << "fixedStep chrom=" << chrom 
	  << "\tstart=" << prev_pos << "\tstep=1" << '\n';
    for (size_t j = prev_pos; j < ends.back().first; ++j)
      out << prev_depth << '\n';
  }
}


int 
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool PRINT_WIG = false;
    string outfile;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("coverage", "A program for calculating target region "
			   "coverage in a resequencing experiment");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("wig", 'w', "print WIG track", 
		      false, PRINT_WIG);
    opt_parse.add_opt("verbose", 'v', "print run info", 
		      false, VERBOSE);
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
    const string mapped_locations_file = leftover_args.front();
    /**********************************************************************/

    if (VERBOSE)
      cerr << "reading mapped locations" << endl;
    vector<SimpleGenomicRegion> mapped;
    ReadBEDFile(mapped_locations_file, mapped);
    assert(check_sorted(mapped));
    if (VERBOSE)
      cerr << "read " << mapped.size() << " locations" << endl;

    vector<size_t> changes;
    find_chrom_changes(mapped.begin(), mapped.end(), changes);
    changes.push_back(mapped.size());

    size_t prev = 0;
    std::ostream *out = (!outfile.empty()) ? 
      new std::ofstream(outfile.c_str()) : &cout;
    for (size_t i = 1; i < changes.size(); ++i) {
      if (VERBOSE)
	cerr << "processing " << mapped[prev].get_chrom() << " ("
	     << (changes[i] - prev) << " reads)" << endl;
      if (PRINT_WIG)
	get_coverage_wig(mapped, prev, changes[i], *out);
      else get_coverage(mapped, prev, changes[i], *out);
      prev = changes[i];
    }
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
