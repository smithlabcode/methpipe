/*    methgram_prep: a program to create the files required for
 *    constructing methgrams
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

/*

This program produces a directory that contains several BED format
files, each containing information about CpG methylation for a
distinct region.

The input is:
(1) BED format file specifying the target regions
(2) BED format file specifying the CpG locations, along with meth and unmeth counts (in the name field)
(4) Value indicating a size of confidence interval
(5) Directory name (must exist) in which to print the results

*/

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include "BSUtils.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

static void
get_meth_unmeth_count(const string name, double &meth, double &unmeth) {
  const vector<string> parts = rmap::split(name, ":");
  unmeth = static_cast<double>(atoi(parts[1].c_str()));
  meth = static_cast<double>(atoi(parts[2].c_str()));
}


static string
reassemble_name(const string &orig_name, const double lower, const double upper) {
  std::ostringstream ss;
  const vector<string> parts = rmap::split(orig_name, ":");
  ss << parts.front() << ":" << parts[1] << ":" << parts[2] 
     << ":" << lower << ":" << upper;
  return ss.str();
}


void
prep_region(const double alpha, const GenomicRegion &region,
	    const vector<GenomicRegion> &cpgs,
	    vector<GenomicRegion> &cpgs_prep) { 
  assert(!cpgs.empty());
  
  const string chrom(cpgs.front().get_chrom());
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const string orig_name(cpgs[i].get_name());
    
    double meth_count = 0, unmeth_count = 0;
    get_meth_unmeth_count(orig_name, meth_count, unmeth_count);
    
    const size_t n_reads = meth_count + unmeth_count;
    const double prop_methylated = meth_count/n_reads;
    
    double lower = 1.0, upper = 0.0;
    wilson_ci_for_binomial(alpha, meth_count + unmeth_count, 
			   prop_methylated, lower, upper);
    
    cpgs_prep.push_back(cpgs[i]);
    cpgs_prep.back().set_name(reassemble_name(orig_name, 
					      lower, upper));
  }
}


int 
main(int argc, const char **argv) {

  try {

    static const string BED_SUFFIX = ".bed";

    bool VERBOSE = false;
    double alpha = 0.1;
    
    string outdir;
    string regions_file;
    
    const string outfiles_suffix = "";
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("methgram_prep", "program to create the files "
			   "required for constructing methgrams");
    opt_parse.add_opt("outdir", 'o', "Name of output dir (default: CWD)", 
		      false, outdir);
    opt_parse.add_opt("alpha", 'a', "alpha value", false , alpha);
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "reading CpG locations" << endl;
    vector<GenomicRegion> cpgs;
    ReadBEDFile(cpgs_file, cpgs);
    assert(check_sorted(cpgs));
    
    if (VERBOSE)
      cerr << "reading regions" << endl;
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    assert(check_sorted(regions));
    
    if (VERBOSE)
      cerr << "separating regions" << endl;
    vector<vector<GenomicRegion> > cpgs_by_region;
    separate_regions(regions, cpgs, cpgs_by_region);
    cpgs.clear();

    for (size_t i = 0; i < regions.size(); ++i) {
      
      vector<GenomicRegion> cpgs_prep;
      prep_region(alpha, regions[i], cpgs_by_region[i], cpgs_prep);
      
      const string region_name(assemble_region_name(regions[i], "_"));
      
      const string filename = 
	path_join(outdir, region_name + 
		  ((outfiles_suffix.empty()) ? 
		   BED_SUFFIX : "_" + outfiles_suffix + BED_SUFFIX));
      std::ofstream out(filename.c_str());
      if (!out) throw RMAPException("could not open file \"" + filename + "\"");
      copy(cpgs_prep.begin(), cpgs_prep.end(), 
	   std::ostream_iterator<GenomicRegion>(out, "\n"));
      out.close();
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
