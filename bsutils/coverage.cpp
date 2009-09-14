/*    coverage: a program for determining target region coverage in a
 *    hybrid-capture experiment
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

static void
get_coverage(const vector<GenomicRegion> &reads, 
	     const GenomicRegion &region, 
	     vector<size_t> &coverage) {
  
  // set individual column data
  const size_t region_start = region.get_start();
  const size_t region_size = region.get_width();
  
  vector<size_t> columns_pos(region_size);
  vector<size_t> columns_neg(region_size);

  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t read_width = reads[i].get_width();
    if (reads[i].pos_strand()) {
      assert(reads[i].get_start() >= region_start);
      const size_t offset = reads[i].get_start() - region_start;
      for (size_t j = 0; j < read_width; ++j)
	columns_pos[offset + j]++;
    }
  }
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t read_width = reads[i].get_width();
    if (reads[i].neg_strand()) {
      const size_t offset = region_size - read_width -
	(reads[i].get_start() - region_start);
      for (size_t j = 0; j < read_width; ++j)
	  columns_neg[offset + j]++;
    }
  }
  reverse(columns_neg.begin(), columns_neg.end());
  
  coverage.resize(region_size);
  for (size_t i = 0; i < columns_pos.size(); ++i)
    coverage[i] = columns_pos[i] + columns_neg[i];
}

static size_t
get_median_from_histogram(vector<size_t> &coverage) {
  size_t total = 0;
  for (size_t i = 0; i < coverage.size(); ++i)
    total += coverage[i];
  size_t counter = 0;
  for (size_t i = 0; i < coverage.size(); ++i) {
    counter += coverage[i];
    if (counter >= total/2) 
      return i;
  }
  return coverage.size();
}

int 
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    
    string outfile;
    string regions_file;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("coverage", "A program for calculating target region "
			   "coverage in a resequencing experiment");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("regions", 'r', "file of relevant regions", 
		      true, regions_file);
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
    vector<GenomicRegion> mapped_locations;
    ReadBEDFile(mapped_locations_file, mapped_locations);
    if (!check_sorted(mapped_locations))
      throw RMAPException("regions in \"" + mapped_locations_file + "\" not sorted");
    
    if (VERBOSE)
      cerr << "read " << mapped_locations.size() << " locations" << endl;
    
    if (VERBOSE)
      cerr << "reading regions file" << endl;
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw RMAPException("regions in \"" + regions_file + "\" not sorted");
    
    if (VERBOSE)
      cerr << "associating reads with regions" << endl;
    vector<vector<GenomicRegion> > clusters;
    separate_regions(regions, mapped_locations, clusters);
    mapped_locations.clear();
    
    if (VERBOSE)
      cerr << "adjusting regions to cover reads" << endl;
    vector<GenomicRegion> adjusted_regions(regions);
    adjust_region_ends(clusters, adjusted_regions);
    
    std::ostream *out = (!outfile.empty()) ? 
      new std::ofstream(outfile.c_str()) : &cout;
    
    vector<size_t> coverage_counts;
    for (size_t i = 0; i < clusters.size(); ++i) {
      vector<size_t> coverage;
      get_coverage(clusters[i], adjusted_regions[i], coverage);
      
      const string chrom(adjusted_regions[i].get_chrom());
      const string name(adjusted_regions[i].get_name());
      const size_t start(adjusted_regions[i].get_start());
      
      size_t j = regions[i].get_start() - adjusted_regions[i].get_start();
      const size_t lim = coverage.size() - 
	(adjusted_regions[i].get_end() - regions[i].get_end());
      for (; j < lim; ++j) {
	*out << chrom << "\t" << start + j << "\t" << start + j + 1 << "\t"
	     << name << "\t" << coverage[j] << "\t+" << endl;
	const size_t c = coverage[j];
	if (c >= coverage_counts.size())
	  coverage_counts.resize(c + 1);
	coverage_counts[coverage[j]]++;
      }
    }
    if (out != &cout) delete out;
    
    if (VERBOSE)
      cerr << "Median read depth:\t" 
	   << get_median_from_histogram(coverage_counts) << endl;
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
