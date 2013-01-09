/*    amrtester: A program for testing whether a genomic region has
 *    allele-specific methylation
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Andrew D. Smith and Fang Fang
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
#include <iomanip>
#include <numeric>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <smithlab_os.hpp>
#include <GenomicRegion.hpp>

#include "Epiread.hpp"
#include "EpireadStats.hpp"
#include "EpireadIO.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


static void
clip_reads(const size_t start_pos, const size_t end_pos, vector<epiread> &r) {
  for (size_t i = 0; i < r.size(); ++i) {
    if (r[i].pos < start_pos) {
      r[i].seq = r[i].seq.substr(start_pos - r[i].pos);
      r[i].pos = start_pos;
    }
    if (r[i].end() > end_pos)
      r[i].seq = r[i].seq.substr(0, end_pos - r[i].pos);
  }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool EPIREAD_FORMAT = false;
    bool PROGRESS = false;
    bool USE_BIC = false;

    string outfile, chroms_dir;
    size_t max_itr = 10;
    double high_prob = 0.75, low_prob = 0.25;

    bool IGNORE_BALANCED_PARTITION_INFO = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "resolve epi-alleles", 
			   "<bed-regions> <mapped-reads>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("no-bal", 'g', "ignore balanced partition info", 
		      false, IGNORE_BALANCED_PARTITION_INFO);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("progress", 'P', "print progress info", false, PROGRESS);
    opt_parse.add_opt("chrom", 'c', "dir of chroms (.fa extn)", false, chroms_dir);
    opt_parse.add_opt("epiread", 'E', "reads in epiread format", false, EPIREAD_FORMAT);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, USE_BIC);
    
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
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string regions_file(leftover_args.front());
    const string reads_file_name(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);

    EpireadIO eio(reads_file_name, VERBOSE, EPIREAD_FORMAT, chroms_dir);
    
    size_t n_regions  = regions.size();
    if (VERBOSE)
      cerr << "NUMBER OF REGIONS: " << n_regions << endl;

    vector<GenomicRegion> amrs;
    for (size_t i = 0; i < regions.size(); ++i) {
      if (PROGRESS) 
	cerr << '\r' << percent(i, n_regions) << "%\r";
      
      vector<epiread> reads;
      eio.load_reads(regions[i], reads);
      
      if (!reads.empty()) {
	clip_reads(regions[i].get_start(), regions[i].get_end(), reads);
	regions[i].set_score((USE_BIC) ?
			     test_asm_bic(max_itr, low_prob, high_prob, reads) :
			     ((IGNORE_BALANCED_PARTITION_INFO) ?
			      test_asm_lrt(max_itr, low_prob, high_prob, reads) :
			      test_asm_lrt2(max_itr, low_prob, high_prob, reads)));
      }
      else regions[i].set_score(1.0);
    }
    if (PROGRESS) cerr << "\r100%" << endl;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    copy(regions.begin(), regions.end(), 
	 std::ostream_iterator<GenomicRegion>(out, "\n"));
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
