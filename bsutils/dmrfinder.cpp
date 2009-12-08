/* dmrfinder: A program for identifying DMRs (differentially
 * methylated regions) based on a file showing probability of
 * differential methylation at each CpG or base.
 *
 * Copyright (C) 2009 University of Southern California
 *                    Andrew D Smith
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

#include "FileIterator.hpp"

#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <fstream>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

using std::numeric_limits;
using std::ostream_iterator;
using std::ofstream;

int
main(int argc, const char **argv) {

  try {

    string outfile;
    double crit = 0.05;
    size_t min_cpgs = 20;
    bool UPPER = false;
    bool VERBOSE = false;

    size_t BUFFER_SIZE = 100000;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("dmrfinder", 
			   "A program for identifying DMRs (differentially "
			   "methylated regions) based on a file showing "
			   "probability of differential methylation at "
			   "each CpG or base. ",
			   "<cpg_meth_diffs_file>");
    opt_parse.add_opt("crit", 'c', "critical value (default: 0.05)", 
		      false, crit);
    opt_parse.add_opt("width", 'w', "width in terms of CpGs of min DMR size (default: 20)", 
		      false, min_cpgs);
    opt_parse.add_opt("upper", 'u', "get upper values", false, UPPER);
    opt_parse.add_opt("out", 'o', "output file (BED format)", 
		      false, outfile);
    opt_parse.add_opt("buffer", 'B', "buffer size (in records, not bytes)", 
		      false , BUFFER_SIZE);
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[READING CPGS]";
    vector<GenomicRegion> cpgs;
    ReadBEDFile(cpgs_file, cpgs);
    if (!check_sorted(cpgs))
      throw RMAPException("file not sorted: \"" + cpgs_file + "\"");
    if (VERBOSE)
      cerr << "[DONE]" << endl;
    
    vector<GenomicRegion> regions;
    size_t count = 0, start = 0, end = 0;
    double total_score = 0.0;
    for (size_t i = 0; i < cpgs.size(); ++i) {
      if (i == 0 || !cpgs[i].same_chrom(cpgs[i - 1])) {
	if (VERBOSE)
	  cerr << "processing:\t" << cpgs[i].get_chrom() << endl;
	if (count > min_cpgs) {
	  const string name("CpG:" + toa(total_score/count));
	  regions.push_back(GenomicRegion(cpgs[i - 1].get_chrom(), 
					  start, end, name, count, '+'));
	}
	count = 0;
	total_score = 0.0;
	start = numeric_limits<size_t>::max();
      }
      const double score = (UPPER) ? 1.0 - cpgs[i].get_score() : cpgs[i].get_score();
      if (score < crit) {
	count += 1;
	total_score += score;
	if (start == numeric_limits<size_t>::max())
	  start = cpgs[i].get_start();
	end = cpgs[i].get_end();
      }
      else {
	if (count > min_cpgs) {
	  const string name("CpG:" + toa(total_score/count));
	  regions.push_back(GenomicRegion(cpgs[i - 1].get_chrom(), 
					  start, end, name, count, '+'));
	}
	count = 0;
	total_score = 0.0;
	start = numeric_limits<size_t>::max();
      }
    }
    if (count > min_cpgs) {
      const string name("CpG:" + toa(total_score/count));
      regions.push_back(GenomicRegion(cpgs.back().get_chrom(), 
				      start, end, name, count, '+'));
    }
    
    std::ostream *out = (outfile.empty()) ? &cout : 
      new std::ofstream(outfile.c_str());
    copy(regions.begin(), regions.end(), 
	 ostream_iterator<GenomicRegion>(*out, "\n"));
    if (out != &cout) delete out;
  }
  catch (RMAPException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
