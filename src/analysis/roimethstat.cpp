/*    roimethstat: compute average methylation in each of a set of regions
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
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <list>
#include <utility>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"


using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;


static void
separate_sites(const vector<GenomicRegion> &regions,
	       const vector<GenomicRegion> &sites, 
	       vector<pair<size_t, size_t> > &sep_sites) {
  const size_t n_regions = regions.size();
  
  for (size_t i = 0; i < n_regions; ++i) {
    GenomicRegion a(regions[i]);
    a.set_end(a.get_start() + 1);
    GenomicRegion b(regions[i]);
    b.set_start(b.get_end());
    b.set_end(b.get_end() + 1);
    
    vector<GenomicRegion>::const_iterator a_insert =
      lower_bound(sites.begin(), sites.end(), a);
    
    vector<GenomicRegion>::const_iterator b_insert =
      lower_bound(sites.begin(), sites.end(), b);
    
    sep_sites.push_back(std::make_pair(a_insert - sites.begin(),
				       b_insert - sites.begin()));
  }
}

static std::pair<size_t, size_t>
region_bounds(const vector<SimpleGenomicRegion> &sites,
              const GenomicRegion &region)
{
  SimpleGenomicRegion a(region);
  a.set_end(a.get_start() + 1);
  vector<SimpleGenomicRegion>::const_iterator a_insert =
    lower_bound(sites.begin(), sites.end(), a);
  
  SimpleGenomicRegion b(region);
  b.set_start(b.get_end());
  b.set_end(b.get_end() + 1);
  vector<SimpleGenomicRegion>::const_iterator b_insert =
    lower_bound(sites.begin(), sites.end(), b);
  
  return std::make_pair(a_insert - sites.begin(),
                        b_insert - sites.begin());
}

static void
get_cpg_stats(const vector<GenomicRegion> &cpgs, 
	      const size_t start_idx, const size_t end_idx,
	      size_t &meth, size_t &reads, size_t &cpgs_with_reads) {
  for (size_t i = start_idx; i < end_idx; ++i) {
    const size_t r = atoi(smithlab::split(cpgs[i].get_name(), 
					  ":").back().c_str());
    meth += static_cast<size_t>(cpgs[i].get_score()*r+0.5);
    // plus 0.5 to make sure the value is rounded correctly
    reads += r;
    cpgs_with_reads += (r > 0);
  }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool PRINT_NAN = false;
    
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("roimethstat", "Compute average CpG "
			   "methylation in each of a set of genomic intervals", 
			   "<intervals-bed> <cpgs-bed>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("print-nan", 'P', "print all records (even if NaN score)", 
		      false, PRINT_NAN);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string regions_file = leftover_args.front();
    const string cpgs_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (methpipe::is_methpipe_file_single(cpgs_file)) {
      if (VERBOSE)
        cerr << "FORMAT = NAME : CPGS : CPGS_WITH_READS : "
          "METH_READS : TOTAL_READS" << endl;
      vector<SimpleGenomicRegion> cpgs;
      vector<pair<double, double> > meths;
      vector<size_t> reads;
      methpipe::load_cpgs(cpgs_file, cpgs, meths, reads);

      vector<GenomicRegion> regions;
      ReadBEDFile(regions_file, regions);
      assert(check_sorted(regions));
      if (!check_sorted(regions))
        throw SMITHLABException("regions not sorted in file: " + regions_file);

      std::ofstream out(outfile.empty() ? "/dev/stdout" : outfile.c_str());

      for (size_t i = 0; i < regions.size(); ++i) {

        const std::pair<size_t, size_t> bounds(region_bounds(cpgs, regions[i]));

        size_t meth = 0, read = 0;
        size_t cpgs_with_reads = 0;
        for (size_t j = bounds.first; j < bounds.second; ++j)
        {
          meth += static_cast<size_t>(meths[j].first);
          read += reads[j];
          cpgs_with_reads += reads[j] > 0;
        }

        const string name = regions[i].get_name() + ":" + 
          toa(bounds.second - bounds.first) + ":" + 
          toa(cpgs_with_reads) + ":" + toa(meth) + ":" + toa(read);
        regions[i].set_name(name);
        regions[i].set_score(static_cast<double>(meth)/read);
        if (PRINT_NAN || std::isfinite(regions[i].get_score()))
          out << regions[i] << endl;
      }
    } else {
    if (VERBOSE)
      cerr << "FORMAT = NAME : CPGS : CPGS_WITH_READS : "
	"METH_READS : TOTAL_READS" << endl;
    
    vector<GenomicRegion> cpgs;
    ReadBEDFile(cpgs_file, cpgs);
    assert(check_sorted(cpgs));
    if (!check_sorted(cpgs))
      throw SMITHLABException("regions not sorted in file: " + cpgs_file);
    
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    assert(check_sorted(regions));
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted in file: " + regions_file);
    
    // separate the CpGs according to the regions in which they occur. 
    vector<pair<size_t, size_t> > roi_cpgs;
    separate_sites(regions, cpgs, roi_cpgs);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    for (size_t i = 0; i < roi_cpgs.size(); ++i) {
      size_t meth = 0, reads = 0;
      size_t cpgs_with_reads = 0;
      get_cpg_stats(cpgs, roi_cpgs[i].first,
		    roi_cpgs[i].second, meth, reads, cpgs_with_reads);
      
      const string name = regions[i].get_name() + ":" + 
	toa(roi_cpgs[i].second - roi_cpgs[i].first) + ":" + 
	toa(cpgs_with_reads) + ":" + toa(meth) + ":" + toa(reads);
      regions[i].set_name(name);
      regions[i].set_score(static_cast<double>(meth)/reads);
      if (PRINT_NAN || std::isfinite(regions[i].get_score()))
	out << regions[i] << endl;
    }
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
