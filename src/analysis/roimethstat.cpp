/*    roimethstat: a program for obtaining methylation statistics
 *    about each of a set of regions of interest (ROIs)
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

// #include <string>
// #include <vector>
// // #include <popt.h>
#include <iostream>
// #include <iterator>
// #include <fstream>
// #include <algorithm>
// #include <numeric>
// #include <list>

// // #include <tr1/unordered_map>

// // #include "OptionParser.hpp"
// // #include "rmap_utils.hpp"
// // #include "rmap_os.hpp"
// // #include "GenomicRegion.hpp"

// // #include <gsl/gsl_statistics_double.h>

// using std::tr1::unordered_map;

// using std::string;
// using std::vector;
// using std::pair;
// using std::make_pair;
// using std::cout;
// using std::cerr;
// using std::endl;
// using std::ifstream;
// using std::istream;
// using std::max;

// void
// separate_regions(const std::vector<GenomicRegion> &big_regions,
// 		 const std::vector<GenomicRegion> &regions, 
// 		 std::vector<std::vector<GenomicRegion> > &sep_regions) {
//   size_t rr_id = 0;
//   const size_t n_regions = regions.size();
//   const size_t n_big_regions = big_regions.size();
//   sep_regions.resize(n_big_regions);
//   for (size_t i = 0; i < n_big_regions; ++i) {
//     const std::string current_chrom(big_regions[i].get_chrom());
//     const size_t current_start = big_regions[i].get_start();
//     const size_t current_end = big_regions[i].get_end();
//     while (rr_id < n_regions &&
// 	   (regions[rr_id].get_chrom() < current_chrom ||
// 	    (regions[rr_id].get_chrom() == current_chrom &&
// 	     regions[rr_id].get_end() < current_start)))
//       ++rr_id;
//     while (rr_id < n_regions &&
// 	   (regions[rr_id].get_chrom() == current_chrom &&
// 	    regions[rr_id].get_start() < current_end)) {
//       sep_regions[i].push_back(regions[rr_id]);
//       ++rr_id;
//     }
//   }
// }



// static void
// get_cpg_stats(const vector<GenomicRegion> &cpgs, 
// 	      size_t &meth,
// 	      size_t &reads,
// 	      size_t &cpgs_with_reads) {
//   for (size_t i = 0; i < cpgs.size(); ++i) {
//     const size_t r = atoi(rmap::split(cpgs[i].get_name(), ":").back().c_str());
//     meth += cpgs[i].get_score()*r;
//     reads += r;
//     cpgs_with_reads += (r > 0);
//   }
// }


int 
main(int argc, const char **argv) {

    std::cerr << "############################################################" << std::endl
              << "#############   roimethstat     ############################" << std::endl
              << "############################################################" << std::endl
              << "THIS PROGRAM roimethstat is deprecated" << std::endl
              << "############################################################" << std::endl
              << "############################################################" << std::endl
              << "############################################################" << std::endl;
    
  // try {
  //   // bool VERBOSE = false;
  //   // bool PRINT_NAN = false;
    
  //   // string outfile;
  //   // string regions_file;
    
  //   // /****************** COMMAND LINE OPTIONS ********************/
  //   // OptionParser opt_parse("roimethstat", "", "<cpgs-bed>");
  //   // opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
  //   // 	      false, outfile);
  //   // opt_parse.add_opt("regions", 'r', "file of regions of interest", 
  //   // 	      false, regions_file);
  //   // opt_parse.add_opt("print-nan", 'P', "print all records (even if NaN score)", 
  //   // 	      false, PRINT_NAN);
  //   // opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  //   // vector<string> leftover_args;
  //   // opt_parse.parse(argc, argv, leftover_args);
  //   // if (argc == 1 || opt_parse.help_requested()) {
  //   //   cerr << opt_parse.help_message() << endl;
  //   //   return EXIT_SUCCESS;
  //   // }
  //   // if (opt_parse.about_requested()) {
  //   //   cerr << opt_parse.about_message() << endl;
  //   //   return EXIT_SUCCESS;
  //   // }
  //   // if (opt_parse.option_missing()) {
  //   //   cerr << opt_parse.option_missing_message() << endl;
  //   //   return EXIT_SUCCESS;
  //   // }
  //   // if (leftover_args.empty()) {
  //   //   cerr << opt_parse.help_message() << endl;
  //   //   return EXIT_SUCCESS;
  //   // }
  //   // const string cpgs_file = leftover_args.front();
  //   // /****************** END COMMAND LINE OPTIONS *****************/
    
  //   // if (VERBOSE)
  //   //   cerr << "format = name:cpgs:cpgs_with_reads:meth:reads" << endl;
    
  //   // vector<GenomicRegion> cpgs;
  //   // ReadBEDFile(cpgs_file, cpgs);
  //   // assert(check_sorted(cpgs));
    
  //   // vector<GenomicRegion> regions;
  //   // ReadBEDFile(regions_file, regions);
  //   // assert(check_sorted(regions));
    
  //   // vector<vector<GenomicRegion> > roi_cpgs;
  //   // separate_regions(regions, cpgs, roi_cpgs);
      
  //   // std::ostream *out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
  //   // for (size_t i = 0; i < roi_cpgs.size(); ++i) {
  //   //   size_t meth = 0;
  //   //   size_t reads = 0;
  //   //   size_t cpgs_with_reads = 0;
  //   //   get_cpg_stats(roi_cpgs[i], meth, reads, cpgs_with_reads);
      
  //   //   const string name = regions[i].get_name() + ":" + toa(roi_cpgs[i].size()) 
  //   // + ":" + toa(cpgs_with_reads) + ":" + 
  //   // toa(meth) + ":" + toa(reads);
  //   //   regions[i].set_name(name);
  //   //   regions[i].set_score(static_cast<double>(meth)/reads);
  //   //   if (PRINT_NAN || std::isfinite(regions[i].get_score()))
  //   // *out << regions[i] << endl;
  //   // }
  //   // if (out != &cout) delete out;
  // }
  // catch (const RMAPException &e) {
  //   cerr << e.what() << endl;
  //   return EXIT_FAILURE;
  // }
  // catch (std::bad_alloc &ba) {
  //   cerr << "ERROR: could not allocate memory" << endl;
  //   return EXIT_FAILURE;
  // }
  // return EXIT_SUCCESS;
}
