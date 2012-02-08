/*    dmr: 
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
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::numeric_limits;


static void
load_cpgs(const bool VERBOSE, string cpgs_file, 
	  vector<SimpleGenomicRegion> &cpgs, vector<double> &prob) {
  if (VERBOSE)
    cerr << "[READING CPGS AND METH PROPS]" << endl;
  vector<GenomicRegion> cpgs_in;
  ReadBEDFile(cpgs_file, cpgs_in);
  if (!check_sorted(cpgs_in))
    throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
  for (size_t i = 0; i < cpgs_in.size(); ++i) {
    cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
    prob.push_back(cpgs_in[i].get_score());
  }
  if (VERBOSE)
    cerr << "TOTAL CPGS: " << cpgs.size() << endl;
}


template <class T> static void
separate_regions(const bool VERBOSE, const size_t desert_size, 
		 vector<SimpleGenomicRegion> &cpgs,
		 vector<T> &meth, vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    cpgs[j] = cpgs[i];
    meth[j] = meth[i];
    ++j;
  }
  cpgs.erase(cpgs.begin() + j, cpgs.end());
  meth.erase(meth.begin() + j, meth.end());
  
  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ? 
      cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].get_start();
  }
  
  if (VERBOSE)
    cerr << "CPGS RETAINED:    " << cpgs.size() << endl
	 << "MEAN REGION SIZE: " << cpgs.size()/reset_points.size() << endl
	 << "DESERTS REMOVED:  " << reset_points.size() << endl << endl;
}


static inline double
weight(const double dist, double bandwidth) {
  const double u = dist/bandwidth;
  return 0.75*(1.0 - u*u);
}


static void
smooth_diff_region(const vector<SimpleGenomicRegion> &cpgs,
		   vector<double> &diffs, 
		   const size_t context_size, 
		   const size_t start, const size_t end) {
  vector<double> updated_diffs(end - start);
  for (size_t i = start; i < end; ++i) {
    double total_diff = 0.0, total_weight = 0.0;
    for (size_t j = ((i >= start + context_size/2) ? i - context_size/2 : start);
	 j < min(i + context_size/2 + context_size%2, end); ++j) {
      const double dist = std::abs(int(i) - int(j));
      const double w = weight(dist, context_size);
      total_diff += w*diffs[j];
      total_weight += w;
    }
    updated_diffs[i - start] = total_diff/total_weight;
  }
  std::swap_ranges(updated_diffs.begin(), updated_diffs.end(), diffs.begin() + start);
}


static void
smooth_diff_scores(const bool VERBOSE, const size_t context_size, 
		   const vector<size_t> &reset_points, 
		   const vector<SimpleGenomicRegion> &cpgs,
		   vector<double> &diffs) {
  size_t prev_percent = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    if (VERBOSE) {
      const size_t curr_percent = percent(i, reset_points.size());
      if (curr_percent != prev_percent)
	cerr << "\rSMOOTHING " << curr_percent << '%';
      prev_percent = curr_percent;
    }
    smooth_diff_region(cpgs, diffs, context_size, reset_points[i], reset_points[i + 1]);
  }
  if (VERBOSE) 
    cerr << "\rSMOOTHING 100%" << endl;
}


static void
get_dmrs(const bool LOWER,
	 const double cutoff,
	 const vector<SimpleGenomicRegion> &cpgs, 
	 const vector<double> &smooth,  
	 const vector<double> &diffs,  
	 const size_t start, const size_t end,
	 vector<GenomicRegion> &dmrs) {
  static const string DMR_LABEL("DMR:");
  bool inside = false;
  size_t dmr_start = 0, n_cpgs = 0;
  double score = 0.0;
  for (size_t i = start; i < end; ++i) {
    if ((!LOWER && (smooth[i] > cutoff)) ||
	(LOWER && (smooth[i] < (1.0 - cutoff)))) {
      if (!inside) {
	inside = true;
	dmr_start = i;
      }
      score += ((LOWER) ? (1.0 - diffs[i]) : diffs[i]);
      ++n_cpgs;
    }
    else if (inside) {
      inside = false;
      dmrs.push_back(GenomicRegion(cpgs[dmr_start]));
      dmrs.back().set_end(cpgs[i - 1].get_end());
      dmrs.back().set_name(DMR_LABEL + toa(n_cpgs));
      dmrs.back().set_score(score);
      score = 0.0;
      n_cpgs = 0;
    }
  }
}


static void
get_dmrs(const bool LOWER,
	 const double cutoff, 
	 const vector<size_t> &reset_points, 
	 const vector<SimpleGenomicRegion> &cpgs, 
	 const vector<double> &smooth,  
	 const vector<double> &diffs,  
	 vector<GenomicRegion> &dmrs) {
  for (size_t i = 0; i < reset_points.size() - 1; ++i)
    get_dmrs(LOWER, cutoff, cpgs, smooth, diffs, reset_points[i], reset_points[i+1], dmrs);
}


static size_t
get_cpgs(const GenomicRegion &dmr) {
  const string name(dmr.get_name());
  return atoi(name.substr(name.find_first_of(":") + 1).c_str());
}


int
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    
    size_t desert_size = 500;
    size_t bandwidth = 10;
    size_t min_size = 200;
    double cutoff = 0.95;
    size_t min_cpgs = 10;
    string outfile;
    
    bool LOWER = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "",
			   "<diffs>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
    opt_parse.add_opt("bw", 'b', "bandwidth", false, bandwidth);
    opt_parse.add_opt("cut", 'c', "score cutoff", false, cutoff);
    opt_parse.add_opt("min", 'm', "minimum size", false, min_size);
    opt_parse.add_opt("cpgs", 'C', "minimum number of CpGs", false, min_cpgs);
    opt_parse.add_opt("lower", 'l', "use lower cutoff (< 1 - cutoff)", false, LOWER);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string diffs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    // separate the regions by chrom and by desert
    vector<SimpleGenomicRegion> cpgs;
    vector<double> diffs;
    load_cpgs(VERBOSE, diffs_file, cpgs, diffs);

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, diffs, reset_points);
    
    vector<double> smooth(diffs);
    smooth_diff_scores(VERBOSE, bandwidth, reset_points, cpgs, smooth);
    
    vector<GenomicRegion> dmrs;
    get_dmrs(LOWER, cutoff, reset_points, cpgs, smooth, diffs, dmrs);
    
    std::ostream *out = (outfile.empty()) ? &cout : 
      new std::ofstream(outfile.c_str());
    for (size_t i = 0; i < dmrs.size(); ++i)
      if (dmrs[i].get_width() >= min_size && get_cpgs(dmrs[i]) >= min_cpgs)
	*out << dmrs[i] << endl;
    if (out != &cout) delete out;
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
