/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
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
 */
 
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "OptionParser.hpp"

#include "locus.hpp"
#include "bin_for_distance.hpp"
#include "combine_pvals.hpp"
#include "fdr.hpp"

using std::cout;  using std::endl;
using std::cerr;  using std::vector;
using std::string;  using std::pair;

int main(int argc, const char **argv) {
  try {
    /* FILES */
    string outfile;
    string bin_spec = "1:200:25";
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("adjust_pval", "a program for computing "
                           "adjust p values using autocorrelation",
                           "<bed-p-values>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("bins", 'b', "corrlation bin specification", 
		      false , bin_spec);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string bed_filename = leftover_args.front();
    /*************************************************************************/
    
    std::ifstream bed_file(bed_filename.c_str());
    
    if (!bed_file)
      throw "could not open file: " + bed_filename;
    
    vector<Locus> loci;
    read_loci(bed_file, loci);
    
    vector<LocusIterator> good_loci_iterators;
    
    cerr << "Getting iterators to loci associated with valid p-values." << endl;
    get_iterators_to_good_loci(loci, good_loci_iterators);
    cerr << "[done]" << endl;
    
    cerr << "Combining p-values." << endl;
    BinForDistance bin_for_dist(bin_spec);
    combine_pvals(good_loci_iterators, bin_for_dist);
    cerr << "[done]" << endl;
    
    cerr << "Running multiple test adjustment." << endl;
    fdr(good_loci_iterators);
    cerr << "[done]" << endl;
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
        
    for (vector<Locus>::const_iterator it = loci.begin(); 
          it != loci.end(); ++it)
      out << *it << std::endl;
    
    //TODO: Check that the regions do not overlap & sorted
        
  } catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
